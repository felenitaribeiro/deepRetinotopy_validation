#!/bin/bash
ml freesurfer/7.3.2
# this script converts all created .mgz files to .gii files

subjects_dir=""
subject_id=""
while getopts s:i: flag
do
    case "${flag}" in
        s) subjects_dir=${OPTARG};;
        i) subject_id=${OPTARG};;
        ?)
           echo "script usage: $(basename "$0") [-s path to subs]" >&2
           exit 1;;
    esac
done

# list of all subjects
if [ -n "$subject_id" ]; then
    echo "Converting files from subject $subject_id"
    file_names=($subject_id)
else
    echo "Converting files from all subjects in the path"
    file_names=()
    for d in $subjects_dir/* ; do
        if [ -d "$d" ]; then
        file_names+=("$(basename "$d")")
    fi
    done
fi

# define hemispheres
hemisphere=("lh" "rh")

for sub in "${file_names[@]}"; do
    
    # path_emp=""$subjects_dir"/$sub/deepRetinotopy/inferred_empirical"
    # path_deep=""$subjects_dir"/$sub/deepRetinotopy/inferred_deepRetinotopy"
    for path in ""$subjects_dir"/$sub/deepRetinotopy/inferred_empirical" ""$subjects_dir"/$sub/deepRetinotopy/inferred_deepRetinotopy" \
            ""$subjects_dir"/$sub/deepRetinotopy/inferred_deepRetinotopy_ones"; do 
        for hemi in "${hemisphere[@]}"; do
            # save all files as .gii files
            params=("angle" "eccen" "sigma" "varea")
            for param in "${params[@]}"; do
                mris_convert -c "$path"/"$hemi".inferred_"$param".mgz "$subjects_dir"/$sub/surf/"$hemi".white "$path"/"$hemi".inferred_"$param".gii
            done
        done
    done
done