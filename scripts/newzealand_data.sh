#!/usr/bin/env bash
ml connectomeworkbench/1.5.0
ml freesurfer/7.3.2
ml deepretinotopy/1.0.19

source ~/miniforge3/etc/profile.d/conda.sh
conda activate deepretinotopy_validation

while getopts d:t:r:o: flag
do
    case "${flag}" in
        d) dataDir=$(realpath -m "${OPTARG}");;   # resolved, since the script changes directory before using it
        t) dirHCP=$(realpath "${OPTARG}");;
        r) validationRepo=$(realpath "${OPTARG}");;
        o) outputDir=$(realpath -m "${OPTARG}");;
	?)
		echo "script usage: $(basename "$0") [-d path to datasets directory] [-t path to directory with HCP template surfaces] [-r path to deepRetinotopy_validation repo] [-o path to output directory]" >&2
		exit 1
	esac
done

# Tests for paths arguments
if [ -z "$dataDir" ] || [ -z "$dirHCP" ] || [ -z "$validationRepo" ] || [ -z "$outputDir" ]; then
    echo "Usage: $(basename "$0") [-d path to datasets directory] [-t path to directory with HCP template surfaces] [-r path to deepRetinotopy_validation repo] [-o path to output directory]"
    exit 1
fi

projectDir="RetinotopyKiwi"
cd "$dataDir"/"$projectDir"/

# Run deepRetinotopy
echo "--------------------------------------------------------------------------------"
echo "[Step 1] Run deepRetinotopy..."
echo "--------------------------------------------------------------------------------"
deepRetinotopy -s "$dataDir"/"$projectDir"/ -t $dirHCP -d kiwi -m "polarAngle,eccentricity,pRFsize" -j 64 -o "$outputDir"

# Convert the retinotopy data to .gii format in the fs_32k space
echo "--------------------------------------------------------------------------------"
echo "[Step 2] Register data from native space to fs_average space..."
echo "--------------------------------------------------------------------------------"
cd "$dataDir"/"$projectDir"/
for subject in `ls .`;
do
    # Skip anything that is not a subject directory (for example logs/)
    if [ ! -d "$subject"/surf ]; then
        continue
    fi
    # Surfaces written by deepRetinotopy in Step 1 (same path as the -o passed above), so that the
    # empirical and predicted maps are resampled with the same template sphere and area surfaces
    surfDir="$outputDir"/"$subject"/surf
    for hemisphere in lh rh;
    do
        if [ $hemisphere == "lh" ]; then
            hemi="L"
        else
            hemi="R"
        fi

        # Polar angle and eccentricity are reconstructed from the resampled x0/y0 maps below (metric y0),
        # so that the 0/360 wrap of the angle is never interpolated
        for metric in x0 y0 sigma r^2;
        do  
            if [ $metric == "sigma" ]; then
                metric_new="pRFsize"
            elif [ $metric == "r^2" ]; then
                metric_new="variance_explained"
            elif [ $metric == "x0" ]; then
                metric_new="x0"
            elif [ $metric == "y0" ]; then
                metric_new="y0"
            fi

            for experiment in CanHrf FitHrf; do
                echo "Resampling native data to fsaverage space..."
                echo "Resampling $metric data..."
                
                echo "Resampling $metric data..."
                wb_command -metric-resample "$subject"/"$hemisphere"_"$experiment"_"$metric".gii \
                        "$surfDir"/"$hemisphere".sphere.reg.surf.gii "$dirHCP"/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                        ADAP_BARY_AREA $surfDir/"$subject".fs_empirical_"$metric_new"_"$experiment"_"$hemisphere".func.gii \
                        -area-surfs "$surfDir"/"$hemisphere".midthickness.surf.gii "$surfDir"/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii   
            

                if [ $metric == "y0" ]; then
                    echo "Reconstructing polar angle and eccentricity from the resampled x0/y0 maps..."
                    reconstruct_coords_native.py \
                                --x $surfDir/"$subject".fs_empirical_x0_"$experiment"_"$hemisphere".func.gii \
                                --y $surfDir/"$subject".fs_empirical_y0_"$experiment"_"$hemisphere".func.gii \
                                --polarangle $surfDir/"$subject".fs_empirical_polarAngle_"$experiment"_"$hemisphere".func.gii \
                                --eccentricity $surfDir/"$subject".fs_empirical_eccentricity_"$experiment"_"$hemisphere".func.gii
                fi
            done
        done
    done
done
