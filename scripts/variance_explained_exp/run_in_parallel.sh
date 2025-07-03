#!/bin/bash

# List of subject IDs
subjects=("sub-01" "sub-02" "sub-03" "sub-04" "sub-05" "sub-06" "sub-07" "sub-08" "sub-09" "sub-10" "sub-11" "sub-12")


# Experiment type
experiment="logbar"
export experiment

# Run the Python script for each subject in parallel using xargs
echo 'python logbar_vexp_exp.py --sub_id "{}" --experiment "$experiment"'
printf "%s\n" "${subjects[@]}" | xargs -I {} -P 5 bash -c 'python logbar_vexp_exp.py --sub_id "{}" --experiment "$experiment"'

# Wait for all background processes to finish
wait

echo "All processes completed."