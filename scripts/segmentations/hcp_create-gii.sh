#!/bin/bash

# this script converts all created .mgz files to .gii files

# list of all subjects
file_names=()
for d in /BULK/LABDATA/openneuro/nyu4christian/HCP/freesurfer/* ; do
    if [ -d "$d" ]; then
    file_names+=("$(basename "$d")")
  fi
done

subjects_dir=""
while getopts s: flag
do
    case "${flag}" in
        s) subjects_dir=${OPTARG};;
        ?)
           echo "script usage: $(basename "$0") [-s path to subs]" >&2
           exit 1;;
    esac
done
# in my case the subjects_dir was "/BULK/LABDATA/openneuro/nyu4christian/HCP/freesurfer"



# define hemispheres
hemisphere=("lh" "rh")

for sub in "${file_names[@]}"; do
    
    path_emp=""$subjects_dir"/$sub/deepRetinotopy/inferred_empirical"
    path_deep=""$subjects_dir"/$sub/deepRetinotopy/inferred_deepRetinotopy"

    
    for hemi in "${hemisphere[@]}"; do

        # save all files as .gii files
        params=("angle" "eccen" "sigma" "varea")
        for param in "${params[@]}"; do
            mris_convert -c "$path_emp"/"$hemi".inferred_"$param".mgz "$subjects_dir"/$sub/surf/"$hemi".white "$path_emp"/"$hemi".inferred_"$param".gii
            wait
            mris_convert -c "$path_deep"/"$hemi".inferred_"$param".mgz "$subjects_dir"/$sub/surf/"$hemi".white "$path_deep"/"$hemi".inferred_"$param".gii
            wait
        done
    done
done