import os
import subprocess
import numpy as np
from pathlib import Path
from argparse import ArgumentParser

def process_empirical_data(args):
    group = 'empirical'
    counter = 0  # Local counter
    
    for sub in args.subjects:
        mgz_results = {}
        
        for hemi in args.hemispheres:
            # Create dictionary of file paths for each hemisphere and parameter
            mgz_results[hemi] = {}
            for param in args.params:
                mgz_results[hemi][param] = "{subjects_dir}/{sub}/deepRetinotopy/{sub}.{group}_{param}.{hemi}.native.func.mgz".format(sub=sub, subjects_dir=args.subjects_dir, hemi=hemi, group=group, param=param)

        out_dir = "{subjects_dir}/{sub}/deepRetinotopy/inferred_empirical".format(subjects_dir=args.subjects_dir, sub=sub)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # Generate the neuropythy command for the current subject (single line)
        register_cmd = "python -m neuropythy register_retinotopy {sub} --subjects-dir={subjects_dir} --verbose --surf-outdir={out_dir} --surf-format=mgz --no-volume-export --max-input-eccen=8 --lh-angle={lh_angle} --lh-eccen={lh_eccen} --lh-radius={lh_pRFsize} --lh-weight={lh_weight} --rh-angle={rh_angle} --rh-eccen={rh_eccen} --rh-radius={rh_pRFsize} --rh-weight={rh_weight}".format(
            sub=sub,
            subjects_dir=args.subjects_dir,
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

        counter += 1
        print("# Command for {sub} ({counter}/{total})".format(sub=sub, counter=counter, total=args.total))
        print(register_cmd)
        print()

        # Execute commands if execute is True
        if args.execute:
            print("Executing register_retinotopy for {sub}".format(sub=sub))
            result = subprocess.run(register_cmd, shell=True, capture_output=True, text=True)

            print("---- STDOUT ----")
            print(result.stdout)

            print("---- STDERR ----")
            print(result.stderr)

def process_deepRetinotopy_data(dummy_weight, args):
    group = "predicted"
    counter = 0  # Local counter
    
    for sub in args.subjects:
        mgz_results = {}
        
        for hemi in args.hemispheres:
            # Create dictionary of file paths for each hemisphere and parameter
            mgz_results[hemi] = {}
            for param in args.params:
                if param == 'weight':
                    param_weight_type = f'{dummy_weight}_{param}'
                    mgz_results[hemi][param] = "{subjects_dir}/{sub}/deepRetinotopy/{sub}.{group}_{param}.{hemi}.native.func.mgz".format(sub=sub, subjects_dir=args.subjects_dir, hemi=hemi, group=group, param=param_weight_type)
                    print(f"Weight file: {mgz_results[hemi]['weight']}")
                else:
                    mgz_results[hemi][param] = "{subjects_dir}/{sub}/deepRetinotopy/{sub}.{group}_{param}.{hemi}.native.func.mgz".format(sub=sub, subjects_dir=args.subjects_dir, hemi=hemi, group=group, param=param)
        
        if dummy_weight == 'mean':
            out_dir = "{subjects_dir}/{sub}/deepRetinotopy/inferred_deepRetinotopy".format(subjects_dir=args.subjects_dir, sub=sub)
        else:
            out_dir = "{subjects_dir}/{sub}/deepRetinotopy/inferred_deepRetinotopy_{dummy_weight}".format(subjects_dir=args.subjects_dir, sub=sub, dummy_weight=dummy_weight)
        
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # Generate the neuropythy command for the current subject (single line)
        register_cmd = "python -m neuropythy register_retinotopy {sub} --subjects-dir={subjects_dir} --verbose --surf-outdir={out_dir} --surf-format=mgz --no-volume-export --max-input-eccen=8 --lh-angle={lh_angle} --lh-eccen={lh_eccen} --lh-radius={lh_pRFsize} --lh-weight={lh_weight} --rh-angle={rh_angle} --rh-eccen={rh_eccen} --rh-radius={rh_pRFsize} --rh-weight={rh_weight}".format(
            sub=sub,
            subjects_dir=args.subjects_dir,
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

        counter += 1
        print("# Command for {sub} with {dummy_weight} weights ({counter}/{total})".format(sub=sub, dummy_weight=dummy_weight, counter=counter, total=args.total))
        print(register_cmd)
        print()

        # Execute commands if execute is True
        if args.execute:
            print("Executing register_retinotopy for {sub} with {dummy_weight} weights".format(sub=sub, dummy_weight=dummy_weight))
            result = subprocess.run(register_cmd, shell=True, capture_output=True, text=True)

            print("---- STDOUT ----")
            print(result.stdout)

            print("---- STDERR ----")
            print(result.stderr)

# Main execution logic
if __name__ == "__main__":

    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'

    parser = ArgumentParser(description="Run neuropythy retinotopy registration")
    parser.add_argument("--execute", action="store_true", help="Execute the commands instead of just printing them")
    parser.add_argument("--subjects_dir", required=True, help="Path to subjects directory")
    parser.add_argument("--subject_id", default=None, help="Subject ID to be processed")
    args = parser.parse_args()

    print(f"Processing subject: {args.subject_id}")
    args.subjects_dir = Path(args.subjects_dir)
    args.params = ["polarAngle_neuropythy", "eccentricity", "pRFsize", "weight"]

    # list all subjects
    if args.subject_id is None:
        subjects = []
        for sub in os.listdir(args.subjects_dir):
            subjects.append(sub)
        args.total = len(subjects)
    else:
        subjects = [args.subject_id]
        args.total = len(subjects)
    args.subjects = subjects
    args.hemispheres=("lh", "rh")

    process_empirical_data(args)
    for dummy_weight in ['ones', 'mean']: 
        process_deepRetinotopy_data(dummy_weight, args)