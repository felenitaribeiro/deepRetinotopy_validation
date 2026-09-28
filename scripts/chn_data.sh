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


projectURL=https://github.com/OpenNeuroDatasets/ds004698.git
projectDir=${projectURL:37:-4} # after slash before .git
cd $dataDir

# git refuses to operate on a repository that appears to be owned by another user, which is the
# case on network shares that squash ownership, and then every datalad call fails silently. Mark
# the dataset location as safe once (the entry goes to the global git config)
if ! git config --global --get-all safe.directory | grep -qxF "$dataDir"/"$projectDir"; then
    git config --global --add safe.directory "$dataDir"/"$projectDir"
fi

datalad install $projectURL

# On the first "datalad get", git-annex probes the GitHub remote for file content, which fails
# and costs that first file (it only marks the remote as annex-ignore afterwards). Record it
# upfront so that every file is fetched from the S3 remote, including the first one
git -C "$dataDir"/"$projectDir" config remote.origin.annex-ignore true

# Dataset download
echo "--------------------------------------------------------------------------------"
echo "[Step 1] Data download..."
echo "--------------------------------------------------------------------------------"
cd "$dataDir"/"$projectDir"/derivatives/freesurfer
for subject in `ls .`; 
do
    if [ ${subject:0:3} != "sub" ]; then
        continue
    else
        for hemisphere in lh rh; 
        do
            if [ $hemisphere == "lh" ]; then
                hemi="L"
            else
                hemi="R"
            fi
            # freesurfer data
            datalad get $subject/surf/"$hemisphere".white
            datalad get $subject/surf/"$hemisphere".pial
            datalad get $subject/surf/"$hemisphere".sphere
            datalad get $subject/surf/"$hemisphere".sphere.reg
            datalad get $subject/surf/"$hemisphere".thickness
            # prf estimates
            datalad get "$dataDir"/"$projectDir"/derivatives/prf-estimation/"$subject"/prfs/"$subject"_ses-all_task-all_hemi-"$hemi"_space-fsnative_prf/*
            datalad get "$dataDir"/"$projectDir"/derivatives/prf-estimation/"$subject"/prfs/"$subject"_ses-all_task-fixedbar_hemi-"$hemi"_space-fsnative_prf/*
            datalad get "$dataDir"/"$projectDir"/derivatives/prf-estimation/"$subject"/prfs/"$subject"_ses-all_task-logbar_hemi-"$hemi"_space-fsnative_prf/*
        done
    fi
done

# Download aperture data
datalad get "$dataDir"/"$projectDir"/derivatives/prf-estimation/stimuli/*

# Run deepRetinotopy
echo "--------------------------------------------------------------------------------"
echo "[Step 2] Run deepRetinotopy..."
echo "--------------------------------------------------------------------------------"
deepRetinotopy -s "$dataDir"/"$projectDir"/derivatives/freesurfer/ -t $dirHCP -d chn -m "polarAngle,eccentricity,pRFsize" -j 64 -o "$outputDir"


# Data processing
echo "--------------------------------------------------------------------------------"
echo "[Step 3] Register data from native space to fs_average space..."
echo "--------------------------------------------------------------------------------"
cd "$dataDir"/"$projectDir"/derivatives/freesurfer
for subject in `ls .`;
do
    if [ ${subject:0:3} != "sub" ]; then
        continue
    else
        # Surfaces written by deepRetinotopy in Step 2 (same path as the -o passed above), so that the
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
            for metric in x0 y0 sigma vexpl;
            do  
                if [ $metric == "sigma" ]; then
                    metric_new="pRFsize"
                elif [ $metric == "vexpl" ]; then
                    metric_new="variance_explained"
                elif [ $metric == "x0" ]; then
                    metric_new="x0"
                elif [ $metric == "y0" ]; then
                    metric_new="y0"
                fi

                for experiment in fixedbar logbar all; do
                    echo "Resampling native data to fsaverage space..."
                    echo "Resampling $metric data..."
                    
                    echo "Convert $metric data to gii format..."
                    mris_convert -c "$dataDir"/"$projectDir"/derivatives/prf-estimation/"$subject"/prfs/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_prf/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_$metric.mgz $subject/surf/"$hemisphere".white \
                    "$dataDir"/"$projectDir"/derivatives/prf-estimation/"$subject"/prfs/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_prf/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_$metric.gii \

                    echo "Resampling $metric data..."
                    wb_command -metric-resample "$dataDir"/"$projectDir"/derivatives/prf-estimation/"$subject"/prfs/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_prf/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_"$metric".gii \
                            "$surfDir"/"$hemisphere".sphere.reg.surf.gii "$dirHCP"/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                            ADAP_BARY_AREA "$surfDir"/"$subject".fs_empirical_"$metric_new"_"$experiment"_"$hemisphere".func.gii \
                            -area-surfs "$surfDir"/"$hemisphere".midthickness.surf.gii "$surfDir"/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii   
              
                    if [ $metric == "y0" ]; then
                        echo "Reconstructing polar angle and eccentricity from the resampled x0/y0 maps..."
                        reconstruct_coords_native.py \
                                --x "$surfDir"/"$subject".fs_empirical_x0_"$experiment"_"$hemisphere".func.gii \
                                --y "$surfDir"/"$subject".fs_empirical_"$metric_new"_"$experiment"_"$hemisphere".func.gii \
                                --polarangle "$surfDir"/"$subject".fs_empirical_polarAngle_"$experiment"_"$hemisphere".func.gii \
                                --eccentricity "$surfDir"/"$subject".fs_empirical_eccentricity_"$experiment"_"$hemisphere".func.gii
                    fi
                done
            done
        done
    fi
done