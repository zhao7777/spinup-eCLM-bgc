import os
import re
import shutil
import argparse

# === ARGPARSE ===
parser = argparse.ArgumentParser(
	description="Clean restart files and optionally remove trash folders.")
parser.add_argument(
    "--clean-trash", action="store_true",
    help="If set, also delete all trash/ folders after moving files."
)
args = parser.parse_args()

# === CONFIGURATION ===
base_rundir = "/p/scratch/cjibg36/jibg3674/eCLM_BGC_SPINUP/postAD"  
DELETE_TRASH = False

# Regex patterns
h0_pattern = re.compile(r'.*\.h0\..*\.nc$')
run_pattern = re.compile(r"^run_(\d{8})_(\d{8})$")

# Restart file prefixes to clean
restart_prefixes = [
    "clmoas.clm2.r.",
    "clmoas.clm2.rh0.",
    "clmoas.cpl.r.",
    "clmoas.datm.rs1."
]

# === MAIN LOOP OVER RUN FOLDERS ===
for folder in sorted(os.listdir(base_rundir)):
    m = run_pattern.match(folder)
    if not m:
        continue
    
    rundir = os.path.join(base_rundir, folder)
    if not os.path.isdir(rundir):
        continue
        
    sdate, edate = m.groups()
    print(f"\nProcessing: {rundir} (start={sdate}, end={edate})")

    # Extract only the year portion of edate (first 4 digits)
    end_year = int(edate[:4])
    last_year = str(end_year + 1)

    trash_dir = os.path.join(rundir, "trash")
    os.makedirs(trash_dir, exist_ok=True)
    
    # Move unwanted files to trash
    for filename in os.listdir(rundir):
        filepath = os.path.join(rundir, filename)

        if os.path.isdir(filepath):
            continue

        if h0_pattern.match(filename):
            # Keep h0 files
            continue

        for prefix in restart_prefixes:
            if filename.startswith(prefix) and last_year not in filename:
                print(f"  Moving old restart: {filename}")
                shutil.move(filepath, os.path.join(trash_dir, filename))
                break
    

# === OPTIONAL CLEANUP ===
if args.clean_trash:
    print("\nCleaning all trash folders...")
    for folder in sorted(os.listdir(base_rundir)):
        if not run_pattern.match(folder):
            continue
            
        rundir = os.path.join(base_rundir, folder)
        trash_dir = os.path.join(rundir, "trash")
        if os.path.exists(trash_dir):
            print(f"  Deleting: {trash_dir}")
            shutil.rmtree(trash_dir)