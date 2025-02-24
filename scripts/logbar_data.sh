ml connectomeworkbench/1.5.0
ml freesurfer/7.3.2
ml deepretinotopy/1.0.8

dataDir=/BULK/ribeiro/datasets
projectURL=https://github.com/OpenNeuroDatasets/ds004698.git
projectDir=${projectURL:37:-4} # after slash before .git
cd $dataDir

datalad install $projectURL

# Dataset download
# Download the required freesurfer data
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
            datalad get "$dataDir"/"$projectDir"/derivatives/prf-estimation/"$subject"/prfs/"$subject"_ses-all_task-logbar_hemi-"$hemi"_space-fsnative_prf/*
            # functional data
            datalad get "$dataDir"/"$projectDir"/derivatives/prf-estimation/"$subject"/func/*
        done
    fi
done

# Download aperture data
datalad get "$dataDir"/"$projectDir"/derivatives/prf-estimation/stimuli/*

# Data processing
# Convert data to .gii format
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

            echo "Generate midthickness surface..."

            mris_convert $subject/surf/"$hemisphere".white $subject/surf/"$hemisphere".white.gii
            mris_convert $subject/surf/"$hemisphere".pial $subject/surf/"$hemisphere".pial.gii
            /home/ribeiro/Projects/deepRetinotopy_validation/deepRetinotopy_TheToolbox/utils/midthickness_surf.py --path $subject/surf/ --hemisphere $hemisphere 
            mris_convert $subject/surf/"$hemisphere".graymid.gii $subject/surf/"$hemisphere".graymid
            mris_curvature -w $subject/surf/"$hemisphere".graymid

            echo "Resampling native surface to fs_LR space..."
            wb_shortcuts -freesurfer-resample-prep $subject/surf/"$hemisphere".white $subject/surf/"$hemisphere".pial \
            $subject/surf/"$hemisphere".sphere.reg /home/ribeiro/Projects/deepRetinotopy_validation/templates/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
            $subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii \
            $subject/surf/"$hemisphere".sphere.reg.surf.gii

            for metric in angle eccen sigma vexpl x0 y0;
            do
                if [ $metric == "angle" ]; then
                    metric_new="polarAngle"
                elif [ $metric == "eccen" ]; then
                    metric_new="eccentricity"
                elif [ $metric == "sigma" ]; then
                    metric_new="pRFsize"
                elif [ $metric == "vexpl" ]; then
                    metric_new="explainedVariance"
                elif [ $metric == "x0" ]; then
                    metric_new="x0"
                elif [ $metric == "y0" ]; then
                    metric_new="y0"
                fi

                echo "Converting $metric data to .gii format..."
                echo `ls "$dataDir"/"$projectDir"/derivatives/prf-estimation/$subject/prfs/"$subject"_ses-all_task-logbar_hemi-"$hemi"_space-fsnative_prf/`
                mris_convert -c "$dataDir"/"$projectDir"/derivatives/prf-estimation/$subject/prfs/"$subject"_ses-all_task-logbar_hemi-"$hemi"_space-fsnative_prf/"$subject"_ses-all_task-logbar_hemi-"$hemi"_space-fsnative_"$metric".mgz $subject/surf/"$hemisphere".white \
                    "$dataDir"/"$projectDir"/derivatives/prf-estimation/$subject/prfs/"$subject"_ses-all_task-logbar_hemi-"$hemi"_space-fsnative_prf/"$subject"_ses-all_task-logbar_hemi-"$hemi"_space-fsnative_"$metric".gii \

                echo "Resampling native data to fsaverage space..."
                wb_command -metric-resample "$dataDir"/"$projectDir"/derivatives/prf-estimation/$subject/prfs/"$subject"_ses-all_task-logbar_hemi-"$hemi"_space-fsnative_prf/"$subject"_ses-all_task-logbar_hemi-"$hemi"_space-fsnative_"$metric".gii \
                    $subject/surf/"$hemisphere".sphere.reg.surf.gii /home/ribeiro/Projects/deepRetinotopy_validation/templates/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                    ADAP_BARY_AREA $subject/surf/"$subject".fs_empirical_"$metric_new"_"$hemisphere".func.gii \
                    -area-surfs $subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii            
            done
        done
    fi
done



# Run deepRetinotopy
deepRetinotopy -s /BULK/ribeiro/datasets/ds004698/derivatives/freesurfer/ -t ~/Projects/deepRetinotopy_validation/templates/ -d logbar -m "polarAngle,eccentricity,pRFsize" 
# signMaps -s /BULK/ribeiro/datasets/ds004698/derivatives/freesurfer/ -t ~/Projects/deepRetinotopy_validation/ -d logbar 
