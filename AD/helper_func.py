import os
from datetime import datetime

def calculate_endtime(sdate, edate, baseunit = 0.25):
    """
    Calculates the end time in hours between two dates.
    """
    sdate = str(sdate)
    edate = str(edate)
    start = datetime.strptime(sdate, "%Y%m%d")
    end = datetime.strptime(edate, "%Y%m%d")
    duration = end - start
    return duration.total_seconds() / 3600.0 + baseunit

import subprocess
import time
from pathlib import Path

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

source loadenvs

date
echo "started" > started.txt
srun eclm.exe
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


# === Function: Submit and monitor job ===
def submit_and_monitor_job(rundir, jobscript, check_interval=1800):
    """
    Submit a SLURM job and monitor its status.
    
    Parameters:
        script_path (str): Path to the SLURM script (e.g., run_job.slurm)
        check_interval (int): Interval in seconds between job status checks
    """
    rundir = Path(rundir).resolve()
    os.chdir(rundir)
    print(f"[INFO] Entered directory: {rundir}") 
    
    result = subprocess.run(["sbatch", jobscript], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] Job submission failed:\n{result.stderr}")
        return
    
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

    return True

