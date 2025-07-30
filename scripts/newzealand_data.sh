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

projectDir="RetinotopyKiwi"
cd "$dataDir"/"$projectDir"/

# Run deepRetinotopy
echo "--------------------------------------------------------------------------------"
echo "[Step 1] Run deepRetinotopy..."
echo "--------------------------------------------------------------------------------"
deepRetinotopy -s "$dataDir"/"$projectDir"/ -t $dirHCP -d kiwi -m "polarAngle,eccentricity,pRFsize"

# Convert the retinotopy data to .gii format in the fs_32k space
echo "--------------------------------------------------------------------------------"
echo "[Step 2] Register data from native space to fs_average space..."
echo "--------------------------------------------------------------------------------"
cd "$dataDir"/"$projectDir"/
for subject in `ls .`;
do
    for hemisphere in lh rh;
    do
        if [ $hemisphere == "lh" ]; then
            hemi="L"
        else
            hemi="R"
        fi

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
                        $subject/surf/"$hemisphere".sphere.reg.surf.gii "$dirHCP"/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                        ADAP_BARY_AREA $subject/surf/"$subject".fs_empirical_"$metric_new"_"$experiment"_"$hemisphere".func.gii \
                        -area-surfs $subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii   
            

                if [ $metric == "y0" ] && [ $hemisphere == "lh" ]; then
                    echo "Generate polar angle and eccentricity data..."
                    python -c "import sys; sys.path.append('"$validationRepo"'); \
                            from functions.preprocess import polarcoord; \
                            polarcoord('"$subject"/"$hemisphere"_"$experiment"_x0.gii', '"$subject"/"$hemisphere"_"$experiment"_y0.gii')"
                    
                    echo "Tranforming polar angle data from lh..."
                    python -c "import sys; sys.path.append('"$validationRepo"'); \
                        from functions.preprocess import transform_angle; \
                        transform_angle('"$subject"/"$hemisphere"_"$experiment"_angle_new.gii', '$hemisphere')"

                    echo "Resampling polar angle data..."
                    wb_command -metric-resample "$subject"/"$hemisphere"_"$experiment"_angle_new_transformed.gii \
                        $subject/surf/"$hemisphere".sphere.reg.surf.gii "$dirHCP"/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                        ADAP_BARY_AREA $subject/surf/"$subject".fs_empirical_polarAngle_"$experiment"_"$hemisphere".func.gii \
                        -area-surfs $subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii   

                    echo "Resampling eccentricity data..."
                    wb_command -metric-resample "$subject"/"$hemisphere"_"$experiment"_eccen_new.gii \
                        $subject/surf/"$hemisphere".sphere.reg.surf.gii "$dirHCP"/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                        ADAP_BARY_AREA $subject/surf/"$subject".fs_empirical_eccentricity_"$experiment"_"$hemisphere".func.gii \
                        -area-surfs $subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii
                        
                elif [ $metric == "y0" ] && [ $hemisphere == "rh" ]; then
                    echo "Generate polar angle and eccentricity data..."
                    python -c "import sys; sys.path.append('"$validationRepo"'); \
                            from functions.preprocess import polarcoord; \
                            polarcoord('"$subject"/"$hemisphere"_"$experiment"_x0.gii',
                                '"$subject"/"$hemisphere"_"$experiment"_y0.gii')"

                    echo "Resampling polar angle data..."
                    wb_command -metric-resample "$subject"/"$hemisphere"_"$experiment"_angle_new.gii \
                        $subject/surf/"$hemisphere".sphere.reg.surf.gii "$dirHCP"/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                        ADAP_BARY_AREA $subject/surf/"$subject".fs_empirical_polarAngle_"$experiment"_"$hemisphere".func.gii \
                        -area-surfs $subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii   

                    echo "Resampling eccentricity data..."
                    wb_command -metric-resample "$subject"/"$hemisphere"_"$experiment"_eccen_new.gii \
                        $subject/surf/"$hemisphere".sphere.reg.surf.gii "$dirHCP"/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                        ADAP_BARY_AREA $subject/surf/"$subject".fs_empirical_eccentricity_"$experiment"_"$hemisphere".func.gii \
                        -area-surfs $subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii
                fi
            done
        done
    done
done
