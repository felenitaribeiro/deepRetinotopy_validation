#!/bin/bash

# list of all subjects
file_names=()
for d in /BULK/LABDATA/openneuro/nyu4christian/HCP/subjects/* ; do
    if [ -d "$d" ]; then
    file_names+=("$(basename "$d")")
  fi
done

## Ausgabe zur Kontrolle
#echo "Gefundene Dateien:"
#for name in "${file_names[@]}"; do
#  echo "$name"
#done
#exit
file_names=("100610")

# path to the methods
path_meth="/home/cbuerger/deepRetinotopy_validation/functions/project"

# define hemispheres
hemisphere=("lh" "rh")

#define groups
groups=("empirical" "predicted")


for sub in "${file_names[@]}"; do
    path_sub="/BULK/LABDATA/openneuro/nyu4christian/HCP/freesurfer/$sub/deepRetinotopy"

    for hem in "${hemisphere[@]}"; do
        for group in "${groups[@]}"; do

            # prework to create mgz files with corrected axis origin
            echo Creating corrected angles for "$group" "$hem" hemisphere in "$sub"...

            # transform the angles
            path_to_use="$path_sub"/"$sub"."$group"_polarAngle."$hem".native.func.gii
            path_to_save="$path_sub"/"$sub"."$group"_polarAngle_neuropythy."$hem".native.func.gii
            cd $path_meth
            python methods.py transform_angle_native_neuropythy --path_to_use "$path_to_use" --path_to_save "$path_to_save" --hemisphere "$hem"
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