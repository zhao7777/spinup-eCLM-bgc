import os
import re
import shutil
import argparse

# === ARGPARSE ===
parser = argparse.ArgumentParser(description="Clean restart files and optionally remove trash folders.")
parser.add_argument(
    "--clean-trash", action="store_true",
    help="If set, also delete all trash/ folders after moving files."
)
args = parser.parse_args()

# === CONFIGURATION ===
start_cycle = 20
total_cycles = 40
base_rundir = "/p/scratch/cjibg36/jibg3674/eCLM_BGC_SPINUP/AD"  
sdate = "19600101"
edate = "19700101"
last_year = "1970"
DELETE_TRASH = False

# Regex patterns
h0_pattern = re.compile(r'.*\.h0\..*\.nc$')

# Restart file prefixes to clean
restart_prefixes = [
    "clmoas.clm2.r.",
    "clmoas.clm2.rh0.",
    "clmoas.cpl.r.",
    "clmoas.datm.rs1."
]

# === MAIN LOOP OVER RUN FOLDERS ===
for cycle in range(start_cycle, total_cycles):
    rundir_name = f"{sdate}_{edate}_cycle{cycle:02d}"
    rundir = os.path.join(base_rundir, rundir_name)

    if not os.path.exists(rundir):
        print(f"Skipping missing folder: {rundir}")
        continue
    
    print(f"\nProcessing: {rundir}")
    trash_dir = os.path.join(rundir, "trash")
    os.makedirs(trash_dir, exist_ok=True)

    # Move unwanted files to trash
    for filename in os.listdir(rundir):
        filepath = os.path.join(rundir, filename)

        if os.path.isdir(filepath):
            continue  
            
        if h0_pattern.match(filename):
            # print(f"  Keeping h0 file:        {filename}")
            continue
            
        for prefix in restart_prefixes:
            if filename.startswith(prefix) and last_year not in filename:
                print(f"  Moving old restart:     {filename}")
                shutil.move(filepath, os.path.join(trash_dir, filename))
                break

# === OPTIONAL CLEANUP ===
if args.clean_trash:
    print("\nCleaning all trash folders...")
    for cycle in range(start_cycle, total_cycles):
        rundir_name = f"{sdate}_{edate}_cycle{cycle:02d}"
        rundir = os.path.join(base_rundir, rundir_name)
        trash_dir = os.path.join(rundir, "trash")
        if os.path.exists(trash_dir):
            print(f"  Deleting: {trash_dir}")
            shutil.rmtree(trash_dir)