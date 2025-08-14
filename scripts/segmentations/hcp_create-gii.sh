#!/bin/bash

# list of all subjects
file_names=()
for d in /BULK/LABDATA/openneuro/nyu4christian/HCP/freesurfer/* ; do
    if [ -d "$d" ]; then
    file_names+=("$(basename "$d")")
  fi
done
#file_names=("100610")

# Ausgabe zur Kontrolle
#echo "Gefundene Dateien:"
#for name in "${file_names[@]}"; do
#  echo "$name"
#done
#exit



# define hemispheres
hemisphere=("lh" "rh")

for sub in "${file_names[@]}"; do
    
    path_emp="/BULK/LABDATA/openneuro/nyu4christian/HCP/freesurfer/$sub/deepRetinotopy/inferred_empirical"
    path_deep="/BULK/LABDATA/openneuro/nyu4christian/HCP/freesurfer/$sub/deepRetinotopy/inferred_deepRetinotopy"
    path_emp_ones="/BULK/LABDATA/openneuro/nyu4christian/HCP/freesurfer/$sub/deepRetinotopy/inferred_empirical_ones"
    path_deep_ones="/BULK/LABDATA/openneuro/nyu4christian/HCP/freesurfer/$sub/deepRetinotopy/inferred_deepRetinotopy_ones"

    
    for hemi in "${hemisphere[@]}"; do

        # save all files as .gii files
        params=("angle" "eccen" "sigma" "varea")
        for param in "${params[@]}"; do
            mris_convert -c "$path_emp"/"$hemi".inferred_"$param".mgz /BULK/LABDATA/openneuro/nyu4christian/HCP/freesurfer/$sub/surf/"$hemi".white "$path_emp"/"$hemi".inferred_"$param".gii
            wait
            mris_convert -c "$path_deep"/"$hemi".inferred_"$param".mgz /BULK/LABDATA/openneuro/nyu4christian/HCP/freesurfer/$sub/surf/"$hemi".white "$path_deep"/"$hemi".inferred_"$param".gii
            wait
            mris_convert -c "$path_emp_ones"/"$hemi".inferred_"$param".mgz /BULK/LABDATA/openneuro/nyu4christian/HCP/freesurfer/$sub/surf/"$hemi".white "$path_emp_ones"/"$hemi".inferred_"$param".gii
            wait
            mris_convert -c "$path_deep_ones"/"$hemi".inferred_"$param".mgz /BULK/LABDATA/openneuro/nyu4christian/HCP/freesurfer/$sub/surf/"$hemi".white "$path_deep_ones"/"$hemi".inferred_"$param".gii
            wait
        done


    done

    ## do register_retinotopy for empiricals
    ## --prior=benson17 --prior-dir=/home/cbuerger/deepRetinotopy_validation/functions/project/priors/ 
    #echo Starting with bensons register_retinotopy for "$sub"...
    #conda activate /home/cbuerger/miniconda3/envs/neuropythy
    #path_free="/BULK/LABDATA/openneuro/nyu4christian/to_use/ds003787/derivatives/freesurfer/"
    #cd $path_free
    #python -m neuropythy register_retinotopy "$sub" --verbose --surf-outdir=. --surf-format=mgz --no-volume-export --lh-angle="$sub"/surf/empirical_lh.angle_neuropythy.mgz --lh-eccen="$sub"/surf/empirical_lh.eccen.mgz --lh-weight="$sub"/surf/empirical_lh.vexpl.mgz --lh-radius="$sub"/surf/empirical_lh.sigma.mgz --rh-angle="$sub"/surf/empirical_rh.angle_neuropythy.mgz --rh-eccen="$sub"/surf/empirical_rh.eccen.mgz --rh-weight="$sub"/surf/empirical_rh.vexpl.mgz --rh-radius="$sub"/surf/empirical_rh.sigma.mgz 
#
    ## wait for all background processes to finish
    #wait
#
    #echo Done for "$sub"
done