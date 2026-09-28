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

projectURL=https://github.com/OpenNeuroDatasets/ds004440.git
projectDir=${projectURL:37:-4} # after slash before .git
cd $dataDir
echo `pwd $dataDir`

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
echo `pwd .`
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
            datalad get "$subject"/surf/"$hemisphere".white
            if [ -L "$subject"/surf/"$hemisphere".pial.T1 ]; then
                datalad get "$subject"/surf/"$hemisphere".pial.T1
                cp "$subject"/surf/"$hemisphere".pial.T1 "$subject"/surf/"$hemisphere".pial
            else
                datalad get "$subject"/surf/"$hemisphere".pial
            fi
            datalad get "$subject"/surf/"$hemisphere".sphere
            datalad get "$subject"/surf/"$hemisphere".sphere.reg
            datalad get "$subject"/surf/"$hemisphere".thickness
        done
    fi
done
# prf estimates
datalad get "$dataDir"/"$projectDir"/derivatives/prfanalyze-vista/children/*
datalad get "$dataDir"/"$projectDir"/derivatives/prfanalyze-vista/adults/*

cd "$dataDir"/"$projectDir"/
datalad unlock .

# Run deepRetinotopy
echo "--------------------------------------------------------------------------------"
echo "[Step 2] Run deepRetinotopy..."
echo "--------------------------------------------------------------------------------"
deepRetinotopy -s "$dataDir"/"$projectDir"/derivatives/freesurfer -t $dirHCP -d stanford -m "polarAngle,eccentricity,pRFsize" -j 64 -o "$outputDir"

# Convert the retinotopy data to .gii format in the fs_32k space
echo "--------------------------------------------------------------------------------"
echo "[Step 3] Register data from native space to fs_average space..."
echo "--------------------------------------------------------------------------------"

for hemisphere in lh rh;
do
    if [ $hemisphere == "lh" ]; then
        hemi="L"
    else
        hemi="R"
    fi
    # Polar angle and eccentricity are not resampled themselves: they are reconstructed from the
    # resampled x/y maps below, so that the 0/360 wrap of the angle is never interpolated
    for metric in sigma vexpl x y;
    do
        if [ $metric == "sigma" ]; then
            metric_new="pRFsize"
        elif [ $metric == "vexpl" ]; then
            metric_new="variance_explained"
        elif [ $metric == "x" ]; then
            metric_new="x0"
        elif [ $metric == "y" ]; then
            metric_new="y0"
        fi

        echo "Converting $metric data to .gii format..."
        for data_folder in adults children; do
            cd "$dataDir"/"$projectDir"/derivatives/prfanalyze-vista/"$data_folder"/
            for subject in `ls .`; do
                # Surfaces written by deepRetinotopy in Step 2 (same path as the -o passed above), so that the
                # empirical and predicted maps are resampled with the same template sphere and area surfaces
                surfDir="$outputDir"/"$subject"/surf
                mris_convert -c "$dataDir"/"$projectDir"/derivatives/prfanalyze-vista/"$data_folder"/"$subject"/"$hemisphere"."$metric".mgz "$dataDir"/"$projectDir"/derivatives/freesurfer/"$subject"/surf/"$hemisphere".white \
                    "$dataDir"/"$projectDir"/derivatives/prfanalyze-vista/"$data_folder"/"$subject"/"$hemisphere"."$metric".gii \
                
                echo "Resampling native data to fsaverage space..."
                wb_command -metric-resample "$dataDir"/"$projectDir"/derivatives/prfanalyze-vista/"$data_folder"/"$subject"/"$hemisphere"."$metric".gii \
                        "$surfDir"/"$hemisphere".sphere.reg.surf.gii "$dirHCP"/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                        ADAP_BARY_AREA "$surfDir"/"$subject".fs_empirical_"$metric_new"_"$hemisphere".func.gii \
                        -area-surfs "$surfDir"/"$hemisphere".midthickness.surf.gii "$surfDir"/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii
                if [ $metric == "y" ]; then
                    echo "Reconstructing polar angle and eccentricity from the resampled x/y maps..."
                    reconstruct_coords_native.py \
                        --x "$surfDir"/"$subject".fs_empirical_x0_"$hemisphere".func.gii \
                        --y "$surfDir"/"$subject".fs_empirical_"$metric_new"_"$hemisphere".func.gii \
                        --polarangle "$surfDir"/"$subject".fs_empirical_polarAngle_"$hemisphere".func.gii \
                        --eccentricity "$surfDir"/"$subject".fs_empirical_eccentricity_"$hemisphere".func.gii
                    echo "Done!"
                fi
            done
        done
    done
done
