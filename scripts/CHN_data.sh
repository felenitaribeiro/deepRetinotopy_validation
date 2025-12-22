#!/usr/bin/env bash
ml connectomeworkbench/1.5.0
ml freesurfer/7.3.2
ml deepretinotopy/1.0.11

source ~/miniforge3/etc/profile.d/conda.sh
conda activate deepretinotopy_validation

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


projectURL=https://github.com/OpenNeuroDatasets/ds004698.git
projectDir=${projectURL:37:-4} # after slash before .git
cd $dataDir

datalad install $projectURL

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
deepRetinotopy -s "$dataDir"/"$projectDir"/derivatives/freesurfer/ -t $dirHCP -d chn -m "polarAngle,eccentricity,pRFsize" 


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
        for hemisphere in lh rh;
        do
            if [ $hemisphere == "lh" ]; then
                hemi="L"
            else
                hemi="R"
            fi

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
                            $subject/surf/"$hemisphere".sphere.reg.surf.gii "$dirHCP"/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                            ADAP_BARY_AREA $subject/surf/"$subject".fs_empirical_"$metric_new"_"$experiment"_"$hemisphere".func.gii \
                            -area-surfs $subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii   
              

                    if [ $metric == "y0" ] && [ $hemisphere == "lh" ]; then
                        echo "Generate polar angle and eccentricity data..."
                        python -c "import sys; sys.path.append('"$validationRepo"'); \
                                from functions.preprocess import polarcoord; \
                                polarcoord('"$dataDir"/"$projectDir"/derivatives/prf-estimation/"$subject"/prfs/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_prf/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_x0.gii', '"$dataDir"/"$projectDir"/derivatives/prf-estimation/"$subject"/prfs/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_prf/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_y0.gii')"
                        
                        echo "Tranforming polar angle data from lh..."
                        python -c "import sys; sys.path.append('"$validationRepo"'); \
                            from functions.preprocess import transform_angle; \
                            transform_angle('"$dataDir"/"$projectDir"/derivatives/prf-estimation/"$subject"/prfs/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_prf/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_angle_new.gii', '$hemisphere')"

                        echo "Resampling polar angle data..."
                        wb_command -metric-resample "$dataDir"/"$projectDir"/derivatives/prf-estimation/"$subject"/prfs/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_prf/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_angle_new_transformed.gii \
                            $subject/surf/"$hemisphere".sphere.reg.surf.gii "$dirHCP"/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                            ADAP_BARY_AREA $subject/surf/"$subject".fs_empirical_polarAngle_"$experiment"_"$hemisphere".func.gii \
                            -area-surfs $subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii   

                        echo "Resampling eccentricity data..."
                        wb_command -metric-resample "$dataDir"/"$projectDir"/derivatives/prf-estimation/"$subject"/prfs/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_prf/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_eccen_new.gii \
                            $subject/surf/"$hemisphere".sphere.reg.surf.gii "$dirHCP"/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                            ADAP_BARY_AREA $subject/surf/"$subject".fs_empirical_eccentricity_"$experiment"_"$hemisphere".func.gii \
                            -area-surfs $subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii
                            
                    elif [ $metric == "y0" ] && [ $hemisphere == "rh" ]; then
                        echo "Generate polar angle and eccentricity data..."
                        python -c "import sys; sys.path.append('"$validationRepo"'); \
                                from functions.preprocess import polarcoord; \
                                polarcoord('"$dataDir"/"$projectDir"/derivatives/prf-estimation/"$subject"/prfs/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_prf/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_x0.gii',
                                    '"$dataDir"/"$projectDir"/derivatives/prf-estimation/"$subject"/prfs/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_prf/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_y0.gii')"

                        echo "Resampling polar angle data..."
                        wb_command -metric-resample "$dataDir"/"$projectDir"/derivatives/prf-estimation/"$subject"/prfs/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_prf/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_angle_new.gii \
                            $subject/surf/"$hemisphere".sphere.reg.surf.gii "$dirHCP"/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                            ADAP_BARY_AREA $subject/surf/"$subject".fs_empirical_polarAngle_"$experiment"_"$hemisphere".func.gii \
                            -area-surfs $subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii   

                        echo "Resampling eccentricity data..."
                        wb_command -metric-resample "$dataDir"/"$projectDir"/derivatives/prf-estimation/"$subject"/prfs/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_prf/"$subject"_ses-all_task-"$experiment"_hemi-"$hemi"_space-fsnative_eccen_new.gii \
                            $subject/surf/"$hemisphere".sphere.reg.surf.gii "$dirHCP"/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                            ADAP_BARY_AREA $subject/surf/"$subject".fs_empirical_eccentricity_"$experiment"_"$hemisphere".func.gii \
                            -area-surfs $subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii
                    fi
                done
            done
        done
    fi
done