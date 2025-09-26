#!/usr/bin/env python3

import os
import sys
import shutil
from datetime import datetime
import importlib
from datetime import datetime, timedelta

sys.path.append("/p/project1/cjibg36/jibg3674/BGC_spinup/postAD")

import create_ensemble_namelists
import helper_func

# importlib.reload(create_ensemble_namelists)
# importlib.reload(helper_func)

from create_ensemble_namelists import *
from helper_func import *

def get_cycle_dates(cycle):
    start_year = base_year + (cycle - 1) * cycle_years
    end_year = start_year + cycle_years - 1
    sdate = int(datetime(start_year, 1, 1).strftime('%Y%m%d'))
    edate = int(datetime(end_year, 12, 31).strftime('%Y%m%d'))
    return sdate, edate

def restart_successful(rundir, edate_str):
    eyear_str = edate_str[:4]
    restart_file = os.path.join(rundir, f"clmoas.clm2.r.{eyear_str}-01-01-00000.nc")
    return os.path.isfile(restart_file)

# === Function: Setup rundir for each cycle ===
def setup_spinup_rundir(cycle, rundir, restart_file, sdate_str, edate_str, atmf_sdate, atmf_edate):
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
    write_drv_in(rundir, sdate_str, edate_str, nprocs, prefix=None)
    write_presaero_stream_files(rundir)
    write_topo_stream_files(rundir)

    components = ["atm", "esp", "glc", "ice", "lnd", "ocn", "rof", "wav"]
    for m in components:
        write_modelio(m, rundir, None)

    write_mosart_in(None, rundir)
    write_datm_in(int(str(atmf_sdate)[:4]), int(str(atmf_edate)[:4]), None, rundir)

    write_lnd_in(rundir, None, restart_from=restart_file)  # restart mode

    # Streams
    streams_ctrl = ["Precip", "Solar", "TPQW"]
    for stream in streams_ctrl:
        write_stream_files(rundir, None, None, str(atmf_sdate), str(atmf_edate), stream)


# Spinup config
total_cycles = 80
prefix = 'eur11_bgc_spinup'
name = "BGC_SPINUP"
nprocs = 384
nreal = 1
ntasks = nprocs * nreal
ntasks_per_node = 48
stime = "06:30:00"
spart = "batch"
sacc = "jibg36"
check_interval = 1800  # 30 minutes

initial_restart = "/p/scratch/cjibg36/jibg3674/eCLM_BGC_SPINUP/AD/19600101_19700101_cycle39/clmoas.clm2.r.1970-01-01-00000.nc "

# Directories
base_rundir = "/p/scratch/cjibg36/jibg3674/eCLM_BGC_SPINUP/postAD/"

base_year = 400  
cycle_years = 10

atmf_sdate = 19600101
atmf_edate = 19691231
atmf_syear = int(str(atmf_sdate)[:4])
atmf_eyear = int(str(atmf_edate)[:4])

# Detect the last completed cycle ---
last_completed_cycle = 0
for cycle in range(1, total_cycles + 1):
    sdate, edate = get_cycle_dates(cycle)
    sdate_str = f"{sdate:08d}"
    edate_str = f"{edate:08d}"
    rundir_name = f"run_{sdate_str}_{edate_str}"
    rundir = os.path.join(base_rundir, rundir_name)
    if restart_file_exists(rundir) and restart_successful(rundir, edate_str):
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
    sdate, edate = get_cycle_dates(cycle)
    sdate_str = f"{sdate:08d}"
    edate_str = f"{edate:08d}"
    rundir_name = f"run_{sdate_str}_{edate_str}"
    rundir = os.path.join(base_rundir, rundir_name)
    os.makedirs(os.path.join(rundir, 'logs'), exist_ok=True)
    os.makedirs(os.path.join(rundir, 'timing', 'checkpoints'), exist_ok=True)

    print(f"\n[{datetime.now()}] Starting cycle {cycle}: {rundir_name}")

    # Set previous rundir
    if cycle == 1:
        restart_file = initial_restart
    if cycle > 1:
        prev_sdate, prev_edate = get_cycle_dates(cycle - 1)
        prev_sdate_str = f"{prev_sdate:08d}"
        prev_edate_str = f"{prev_edate:08d}"
        prev_rundir_name = f"run_{prev_sdate_str}_{prev_edate_str}"
        prev_rundir = os.path.join(base_rundir, prev_rundir_name)
        restart_file = os.path.join(prev_rundir, f"clmoas.clm2.r.{sdate_str[:4]}-01-01-00000.nc")        
        while not restart_file_exists(prev_rundir):
        	print(f"[{datetime.now()}] Restart file is missing. Awaiting completion of cycle {cycle - 1} in {prev_rundir}...")

    # Setup rundir
    setup_spinup_rundir(cycle, rundir, restart_file, sdate_str, edate_str, atmf_sdate, atmf_edate)

    # Create SLURM job script
    jobscript = os.path.join(rundir, "jobscript.slurm")
    with open(jobscript, "w") as f:
        f.write(generate_slurm_script(name, nprocs, nreal, ntasks_per_node, stime, spart, sacc))

    # Submit job and monitor
    success = submit_and_monitor_job(rundir, jobscript, check_interval=check_interval)
    if not success or not restart_successful(rundir, edate_str):
        print(f"[{datetime.now()}] ERROR: Cycle {cycle} failed.")
        break
    else:
        print(f"[{datetime.now()}] Cycle {cycle} completed successfully.")
