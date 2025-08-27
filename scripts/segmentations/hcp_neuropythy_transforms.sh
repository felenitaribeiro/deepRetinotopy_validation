#!/bin/bash
ml freesurfer/7.3.2
# this script transforms the anlges from deepRetinotopy to the required angles for neuropythy and save them as ..._neuropythy

subjects_dir=""
while getopts s:r: flag
do
    case "${flag}" in
        s) subjects_dir=${OPTARG};;
        r) path_to_validation_repo=${OPTARG};;
        ?)
           echo "script usage: $(basename "$0") [-s path to subs] [-r path to validation repo]" >&2
           exit 1;;
    esac
done

# list of all subjects
file_names=()
for d in "$subjects_dir"/* ; do
    if [ -d "$d" ]; then
    file_names+=("$(basename "$d")")
  fi
done

# define hemispheres
hemisphere=("lh" "rh")

#define groups
groups=("empirical" "predicted")


for sub in "${file_names[@]}"; do
    path_sub=""$subjects_dir"/$sub/deepRetinotopy"

    for hem in "${hemisphere[@]}"; do
        for group in "${groups[@]}"; do

            # prework to create mgz files with corrected axis origin
            echo Creating corrected angles for "$group" "$hem" hemisphere in "$sub"...

            # transform the angles
            path_to_use="$path_sub"/"$sub"."$group"_polarAngle."$hem".native.func.gii
            path_to_save="$path_sub"/"$sub"."$group"_polarAngle_neuropythy."$hem".native.func.gii
            python $path_to_validation_repo/functions/preprocess.py transform_polarangle_to_benson14 --path_to_use "$path_to_use" --path_to_save "$path_to_save" --hemisphere $hem
            # wait for all background processes to finish
            wait

            # save all files as .mgz files
            params=("eccentricity" "polarAngle_neuropythy" "pRFsize")
            cd $path_sub
            for param in "${params[@]}"; do
                mri_convert "$sub"."$group"_"$param"."$hem".native.func.gii "$sub"."$group"_"$param"."$hem".native.func.mgz
                # wait for all background processes to finish
                wait
            done
            if [ "$group" = "predicted" ]; then
                mri_convert "$sub".dummy_variance_explained."$hem".native.func.gii "$sub".predicted_weight."$hem".native.func.mgz
            else
                mri_convert "$sub".empirical_variance_explained."$hem".native.func.gii "$sub".empirical_weight."$hem".native.func.mgz
            fi
        done
    done
done