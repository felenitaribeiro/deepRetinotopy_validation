import os
#import neuropythy
from argparse import ArgumentParser
import subprocess
import numpy as np
from pathlib import Path


os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

parser = ArgumentParser(description="Run neuropythy retinotopy registration")
parser.add_argument("--execute", action="store_true", help="Execute the commands instead of just printing them")
parser.add_argument("--subjects_dir", required=True, help="Path to subjects directory")
args = parser.parse_args()
execute = args.execute

subjects_dir = Path(args.subjects_dir)
# in my case the subjects_dir was "/BULK/LABDATA/openneuro/nyu4christian/HCP/freesurfer"

#
#!!! be aware of the weights. I named them exactly like the other parameter files to make it easier
#
params = ["polarAngle_neuropythy", "eccentricity", "pRFsize", "weight"]
groups = ["predicted", "empirical"]


# list all subjects
subjects = []
for sub in os.listdir(subjects_dir):
    subjects.append(sub)
total = len(subjects)


hemispheres=("lh", "rh")
counter = 0
for group in groups:
    for sub in subjects:
        mgz_results = {}
        
        for hemi in hemispheres:
            # Create dictionary of file paths for each hemisphere and parameter
            mgz_results[hemi] = {}
            for param in params:
                mgz_results[hemi][param] = "{subjects_dir}/{sub}/deepRetinotopy/{sub}.{group}_{param}.{hemi}.native.func.mgz".format(sub=sub, subjects_dir=subjects_dir, hemi=hemi, group=group, param=param)


        if group == "predicted":
            out_dir = "{subjects_dir}/{sub}/deepRetinotopy/inferred_deepRetinotopy".format(subjects_dir=subjects_dir, sub=sub)
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
        else:
            out_dir = "{subjects_dir}/{sub}/deepRetinotopy/inferred_empirical".format(subjects_dir=subjects_dir, sub=sub)
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
            lh_weight=mgz_results['lh']['weight'],
            rh_angle=mgz_results['rh']['polarAngle_neuropythy'],
            rh_eccen=mgz_results['rh']['eccentricity'],
            rh_pRFsize=mgz_results['rh']['pRFsize'],
            rh_weight=mgz_results['rh']['weight']
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
