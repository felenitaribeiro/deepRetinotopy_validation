import os
import neuropythy
from argparse import ArgumentParser
import subprocess
import shutil
import sys
sys.path.append("/home/cbuerger/deepRetinotopy_validation/functions/project")
import methods
import time
import numpy as np


os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

parser = ArgumentParser(description="Run neuropythy retinotopy registration")
parser.add_argument("--execute", action="store_true", help="Execute the commands instead of just printing them")
args = parser.parse_args()
execute = args.execute


subjects_dir="/BULK/LABDATA/openneuro/nyu4christian/HCP/freesurfer"

params = ["polarAngle_neuropythy", "eccentricity", "pRFsize", "ones"]
groups = ["predicted"]
# , "empirical"

subjects_done = []
for sub in os.listdir("/BULK/LABDATA/openneuro/nyu4christian/HCP/subjects"):
    subjects_done.append(sub)

subjects = []
for sub in os.listdir(subjects_dir):
    subjects.append(sub)
subjects = np.setdiff1d(subjects, subjects_done)
total = len(subjects)
#subjects = ["100610"]
#, "102311"


hemispheres=("lh", "rh")
counter = 0
for group in groups:
    for sub in subjects:
        ## Run the Benson14 retinotopy before register_retinotopy
        #benson_cmd = "python -m neuropythy benson14_retinotopy --verbose --surf-format=mgz {subjects_dir}/{sub}/deepRetinotopy".format(subjects_dir=subjects_dir, sub=sub)
        #print("# Benson14 Retinotopy for {sub}".format(sub=sub))
        #print(benson_cmd)
        #print()



        mgz_results = {}
        
        for hemi in hemispheres:
            # Create dictionary of file paths for each hemisphere and parameter
            mgz_results[hemi] = {}
            for param in params:
                mgz_results[hemi][param] = "{subjects_dir}/{sub}/deepRetinotopy/{sub}.{group}_{param}.{hemi}.native.func.mgz".format(sub=sub, subjects_dir=subjects_dir, hemi=hemi, group=group, param=param)


        if group == "predicted":
            out_dir = "{subjects_dir}/{sub}/deepRetinotopy/inferred_deepRetinotopy_ones".format(subjects_dir=subjects_dir, sub=sub)
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
        else:
            out_dir = "{subjects_dir}/{sub}/deepRetinotopy/inferred_empirical_ones".format(subjects_dir=subjects_dir, sub=sub)
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)


        # Generate the neuropythy command for the current subject
        register_cmd = """python -m neuropythy \\
        register_retinotopy {sub} \\
        --subjects-dir={subjects_dir} \\
        --verbose \\
        --surf-outdir={out_dir} \\
        --surf-format="mgz" \\
        --no-volume-export \\
        --max-input-eccen=8 \\
        --lh-angle={lh_angle} \\
        --lh-eccen={lh_eccen} \\
        --lh-radius={lh_pRFsize} \\
        --lh-weight={lh_weight} \\
        --rh-angle={rh_angle} \\
        --rh-eccen={rh_eccen} \\
        --rh-radius={rh_pRFsize} \\
        --rh-weight={rh_weight}""".format(
            sub=sub,
            subjects_dir=subjects_dir,
            out_dir=out_dir,
            lh_angle=mgz_results['lh']['polarAngle_neuropythy'],
            lh_eccen=mgz_results['lh']['eccentricity'],
            lh_pRFsize=mgz_results['lh']['pRFsize'],
            lh_weight=mgz_results['lh']['ones'],
            rh_angle=mgz_results['rh']['polarAngle_neuropythy'],
            rh_eccen=mgz_results['rh']['eccentricity'],
            rh_pRFsize=mgz_results['rh']['pRFsize'],
            rh_weight=mgz_results['rh']['ones']
        )

        print("# Command for {sub}".format(sub=sub))
        print(register_cmd)
        print()
        counter = counter + 1 
        print("# {counter} of {total}".format(counter=counter, total=total))

        # Execute commands if execute is True
        if execute:
            # print("Executing Benson14 retinotopy for {sub}...".format(sub=sub))
            # subprocess.run(benson_cmd, shell=True, check=True)

            print("Executing register_retinotopy for {sub}".format(sub=sub))
            result = subprocess.run(register_cmd, shell=True, capture_output=True, text=True)

            print("---- STDOUT ----")
            print(result.stdout)

            print("---- STDERR ----")
            print(result.stderr)
