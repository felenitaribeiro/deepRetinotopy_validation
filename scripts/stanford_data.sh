#!/usr/bin/env bash
ml connectomeworkbench/1.5.0
ml freesurfer/7.3.2
ml deepretinotopy/1.0.11

while getopts d:t:r: flag
do
    case "${flag}" in
        d) dataDir=${OPTARG};;
        t) dirHCP=$(realpath "${OPTARG}");;
        r) validationRepo=${OPTARG};;
	?)
		echo "script usage: $(basename "$0") [-d path to datasets directory] [-t path to directory with HCP template surfaces] [-r path to deepRetinotopy_validation repo]" >&2
		exit 1
	esac
done

# Tests for paths arguments
if [ -z "$dataDir" ] || [ -z "$dirHCP" ] || [ -z "$validationRepo" ]; then
    echo "Usage: $(basename "$0") [-d path to datasets directory] [-t path to directory with HCP template surfaces] [-r path to deepRetinotopy_validation repo]"
    exit 1
fi

projectURL=https://github.com/OpenNeuroDatasets/ds004440.git
projectDir=${projectURL:37:-4} # after slash before .git
cd $dataDir
echo `pwd $dataDir`

datalad install $projectURL

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
datalad get "$dataDir"/"$projectDir"/derivatives/prfanalyze-vista/"$data_folder"/children/*
datalad get "$dataDir"/"$projectDir"/derivatives/prfanalyze-vista/"$data_folder"/adults/*

cd "$dataDir"/"$projectDir"/
datalad unlock .

# Run deepRetinotopy
echo "--------------------------------------------------------------------------------"
echo "[Step 2] Run deepRetinotopy..."
echo "--------------------------------------------------------------------------------"
deepRetinotopy -s "$dataDir"/"$projectDir"/derivatives/freesurfer -t $dirHCP -d stanford -m "polarAngle,eccentricity,pRFsize"


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
    for metric in eccen sigma vexpl angle;
    do
        if [ $metric == "angle" ]; then
            metric_new="polarAngle"
        elif [ $metric == "eccen" ]; then
            metric_new="eccentricity"
        elif [ $metric == "sigma" ]; then
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
                mris_convert -c "$dataDir"/"$projectDir"/derivatives/prfanalyze-vista/"$data_folder"/"$subject"/"$hemisphere"."$metric".mgz "$dataDir"/"$projectDir"/derivatives/freesurfer/"$subject"/surf/"$hemisphere".white \
                    "$dataDir"/"$projectDir"/derivatives/prfanalyze-vista/"$data_folder"/"$subject"/"$hemisphere"."$metric".gii \
                
                echo "Resampling native data to fsaverage space..."
                # Transform polar angle data before resampling
                if [ $metric == "angle" ]; then
                    python -c "import sys; sys.path.append('"$validationRepo"/'); from functions.preprocess import transform_angle; transform_angle('"$dataDir"/"$projectDir"/derivatives/prfanalyze-vista/"$data_folder"/"$subject"/"$hemisphere"."$metric".gii', '$hemisphere', radians = True)"
                    wb_command -metric-resample "$dataDir"/"$projectDir"/derivatives/prfanalyze-vista/"$data_folder"/"$subject"/"$hemisphere"."$metric"_transformed.gii \
                        "$dataDir"/"$projectDir"/derivatives/freesurfer/"$subject"/surf/"$hemisphere".sphere.reg.surf.gii "$dirHCP"/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                        ADAP_BARY_AREA "$dataDir"/"$projectDir"/derivatives/freesurfer/"$subject"/surf/"$subject".fs_empirical_"$metric_new"_"$hemisphere".func.gii \
                        -area-surfs "$dataDir"/"$projectDir"/derivatives/freesurfer/"$subject"/surf/"$hemisphere".midthickness.surf.gii "$dataDir"/"$projectDir"/derivatives/freesurfer/"$subject"/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii      
                    echo "Done!"
               else
                    wb_command -metric-resample "$dataDir"/"$projectDir"/derivatives/prfanalyze-vista/"$data_folder"/"$subject"/"$hemisphere"."$metric".gii \
                        "$dataDir"/"$projectDir"/derivatives/freesurfer/"$subject"/surf/"$hemisphere".sphere.reg.surf.gii "$dirHCP"/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                        ADAP_BARY_AREA "$dataDir"/"$projectDir"/derivatives/freesurfer/"$subject"/surf/"$subject".fs_empirical_"$metric_new"_"$hemisphere".func.gii \
                        -area-surfs "$dataDir"/"$projectDir"/derivatives/freesurfer/"$subject"/surf/"$hemisphere".midthickness.surf.gii "$dataDir"/"$projectDir"/derivatives/freesurfer/"$subject"/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii      
                fi
            done
        done
    done
done
