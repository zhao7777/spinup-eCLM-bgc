from pathlib import Path
import subprocess
import time
import os
from datetime import datetime


def calculate_endtime(sdate, edate, baseunit=0.25):
    """
    Calculates the end time in hours between two dates.
    """
    sdate = str(sdate)
    edate = str(edate)
    start = datetime.strptime(sdate, "%Y%m%d")
    end = datetime.strptime(edate, "%Y%m%d")
    duration = end - start
    return duration.total_seconds() / 3600.0 + baseunit


# === Function: Check if a restart was successful ===
def restart_file_exists(rundir):
    return Path(rundir, "ready.txt").exists()


# === Function: generate SLURM script ===
def generate_slurm_script(name, nprocs, nreal, ntasks_per_node, stime, spart, sacc):
    ntasks = nprocs * nreal
    nnodes = ntasks // ntasks_per_node

    slurm_template = """#!/bin/bash

#SBATCH --job-name="{jobname}"
#SBATCH --nodes={nnodes}
#SBATCH --ntasks={ntasks}
#SBATCH --ntasks-per-node={ntasks_per_node}
#SBATCH --output=mpiMPMD-out.%j
#SBATCH --error=mpiMPMD-err.%j
#SBATCH --time={stime}
#SBATCH --partition={spart}
#SBATCH --mail-type=NONE
#SBATCH --account={sacc}

export PSP_RENDEZVOUS_OPENIB=-1

echo "PWD: $(pwd)"
source ./loadenvs

date
echo "started" > started.txt
srun ./eclm.exe
date
echo "ready" > ready.txt
exit 0
"""
    filled_script = slurm_template.format(
        jobname=name,
        nnodes=nnodes,
        ntasks=ntasks,
        ntasks_per_node=ntasks_per_node,
        stime=stime,
        spart=spart,
        sacc=sacc
    )

    return filled_script


# === Helper: Check for spawn failure in slurm log ===
def check_spawn_failure(rundir, job_id):
    logfile = Path(rundir) / f"mpiMPMD-out.{job_id}"
    if not logfile.exists():
        return False
    with open(logfile, "r", errors="ignore") as f:
        text = f.read()
    spawn_signatures = [
        "PSI: doSpawn: spawn to node",
        "Could not spawn './eclm.exe'",
        "PSI_spawnRsrvtn() failed",
        "Unable to start all processes. Aborting"
    ]
    return any(sig in text for sig in spawn_signatures)


# === Function: Submit and monitor job ===
def submit_and_monitor_job(rundir, jobscript, check_interval=1800, max_retries=2):
    """
    Submit a SLURM job and monitor its status, and retry if spawn failure occurs.

    Parameters:
        rundir (str|Path): Directory where job runs
        jobscript (str|Path): Path to the SLURM script
        check_interval (int): Interval in seconds between job status checks
        max_retries (int): Maximum number of retries if spawn failure occurs
    """

    rundir = Path(rundir).resolve()
    os.chdir(rundir)
    print(f"[INFO] Entered directory: {rundir}")

    retries = 0
    while retries <= max_retries:
        result = subprocess.run(["sbatch", jobscript], capture_output=True, text=True)
        if result.returncode != 0:
            print(f"[ERROR] Job submission failed:\n{result.stderr}")
            return False

        job_id = result.stdout.strip().split()[-1]
        print(f"[INFO] Job submitted: {job_id}")
        print(f"[INFO] Monitoring job every {check_interval} seconds...\n")

        while True:
            status = subprocess.run(["squeue", "-j", job_id], capture_output=True, text=True)
            if job_id in status.stdout:
                print(f"[{time.strftime('%H:%M:%S')}] Job {job_id} still running...")
            else:
                print(f"[{time.strftime('%H:%M:%S')}] Job {job_id} has finished.")
                break
            time.sleep(check_interval)

        # --- Check for spawn failure ---
        if check_spawn_failure(rundir, job_id):
            print(f"[{datetime.now()}] Spawn failure detected in job {job_id}. Retrying...")
            retries += 1
            continue
        else:
            # If job finished without spawn error → return True
            return True

    print(f"[{datetime.now()}] Exceeded maximum retries ({max_retries}) due to spawn failures.")
    return False