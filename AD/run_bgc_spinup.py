#!/usr/bin/env python3

import os
import sys
import shutil
from datetime import datetime
import importlib

sys.path.append("/p/project1/cjibg36/jibg3674/BGC_spinup/")

import create_ensemble_namelists
import helper_func

# importlib.reload(create_ensemble_namelists)
# importlib.reload(helper_func)

from create_ensemble_namelists import *
from helper_func import *


def restart_successful(rundir):
    restart_file = os.path.join(rundir, "clmoas.clm2.r.1970-01-01-00000.nc")
    return os.path.isfile(restart_file)

# === Function: Setup rundir for each cycle ===
def setup_spinup_rundir(cycle, rundir, prev_rundir, sdate, edate):
    os.makedirs(rundir, exist_ok=True)
    os.makedirs(os.path.join(rundir, 'logs'), exist_ok=True)
    os.makedirs(os.path.join(rundir, 'timing', 'checkpoints'), exist_ok=True)

    # Link build and env
    file_build = '/p/project/cjibg36/jibg3674/eCLM_PyDA/TSMP2/bin/JUWELS_eCLM/bin/eclm.exe'
    link_build = os.path.join(rundir, "eclm.exe")
    if os.path.islink(link_build) or os.path.exists(link_build):
        os.remove(link_build)
    os.symlink(file_build, link_build)
    file_env = '/p/project1/cjibg36/jibg3674/eCLM_PyDA/TSMP2/env/jsc.2024_Intel.sh'
    link_env = os.path.join(rundir, "loadenvs")
    if os.path.islink(link_env) or os.path.exists(link_env):
        os.remove(link_env)
    os.symlink(file_env, link_env)

    # Write input files
    write_cpl_modelio(rundir)
    write_drv_flds_in(rundir)
    write_drv_in(rundir, sdate, edate, nprocs, prefix=None)
    write_presaero_stream_files(rundir)
    write_topo_stream_files(rundir)

    components = ["atm", "esp", "glc", "ice", "lnd", "ocn", "rof", "wav"]
    for m in components:
        write_modelio(m, rundir, None)

    write_mosart_in(None, rundir)
    write_datm_in(int(str(sdate)[:4]), int(str(edate)[:4]), None, rundir)

    write_lnd_in(rundir, None, restart_from=prev_rundir)  # restart mode

    # Streams
    streams_ctrl = ["Precip", "Solar", "TPQW"]
    for stream in streams_ctrl:
        write_stream_files(rundir, None, None, sdate, edate, stream)


# Simulation time range
sdate = 19600101
edate = 19700101
syear = int(str(sdate)[:4])
eyear = int(str(edate)[:4])

# Spinup config
total_cycles = 20
prefix = 'eur11_bgc_spinup'
name = "BGC_SPINUP"
nprocs = 384
nreal = 1
ntasks = nprocs * nreal
ntasks_per_node = 48
stime = "05:00:00"
spart = "batch"
sacc = "jibg36"
check_interval = 1800  # 30 minutes

# Directories
base_rundir = "/p/scratch/cjibg36/jibg3674/eCLM_BGC_SPINUP/AD/"

# Detect the last completed cycle ---
last_completed_cycle = 0
for cycle in range(1, total_cycles + 1):
    rundir_name = f"{sdate}_{edate}_cycle{cycle:02d}"
    rundir = os.path.join(base_rundir, rundir_name)
    if restart_file_exists(rundir):
        last_completed_cycle = cycle
    else:
        break

# Decide where to start ---
start_cycle = last_completed_cycle + 1
if start_cycle > total_cycles:
    print(f"[{datetime.now()}] All {total_cycles} cycles already completed. Exiting.")
    exit()
print(f"[{datetime.now()}] Resuming from cycle {start_cycle} of {total_cycles}...")


# Main spinup loop
for cycle in range(start_cycle, total_cycles + 1):
    rundir_name = f"{sdate}_{edate}_cycle{cycle:02d}"
    rundir = os.path.join(base_rundir, rundir_name)
    os.makedirs(os.path.join(rundir, 'logs'), exist_ok=True)
    os.makedirs(os.path.join(rundir, 'timing', 'checkpoints'), exist_ok=True)

    print(f"\n[{datetime.now()}] Starting cycle {cycle}: {rundir_name}")

    # Set previous rundir
    prev_rundir = os.path.join(base_rundir, f"{sdate}_{edate}_cycle{cycle - 1:02d}")
    if cycle > 1:
        while not restart_file_exists(prev_rundir):
        	print(f"[{datetime.now()}] Restart file is missing. Awaiting completion of cycle {cycle - 1}...")

    # Setup rundir
    setup_spinup_rundir(cycle, rundir, prev_rundir, sdate, edate)

    # Create SLURM job script
    jobscript = os.path.join(rundir, "jobscript.slurm")
    with open(jobscript, "w") as f:
        f.write(generate_slurm_script(name, nprocs, nreal, ntasks_per_node, stime, spart, sacc))

    # Submit job and monitor
    success = submit_and_monitor_job(rundir, jobscript, check_interval=check_interval)
    if not success or not restart_successful(rundir):
        print(f"[{datetime.now()}] ERROR: Cycle {cycle} failed.")
        break
    else:
        print(f"[{datetime.now()}] Cycle {cycle} completed successfully.")