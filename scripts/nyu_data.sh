ml connectomeworkbench/1.5.0
ml freesurfer/7.3.2
ml deepretinotopy/1.0.8

dataDir=/BULK/LABDATA/openneuro/
projectURL=https://github.com/OpenNeuroDatasets/ds003787.git
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
            datalad get $subject/surf/"$hemisphere".white
            datalad get $subject/surf/"$hemisphere".pial
            datalad get $subject/surf/"$hemisphere".sphere
            datalad get $subject/surf/"$hemisphere".sphere.reg
            datalad get $subject/surf/"$hemisphere".thickness
            # prf estimates
            datalad get "$dataDir"/"$projectDir"/derivatives/prfanalyze-vista/$subject/ses-nyu3t01/*
            # functional data
            datalad get $dataDir/"$projectDir"/derivatives/fmriprep/$subject/ses-nyu3t01/func/*
            # stimulus apertures
            datalad get $dataDir/"$projectDir"/derivatives/stimulus_apertures/"$subject"/ses-nyu3t01/*.mat
        done
    fi
done

cd "$dataDir"/"$projectDir"/
datalad unlock .

# Convert the retinotopy data to .gii format in the fs_32k space
echo "--------------------------------------------------------------------------------"
echo "[Step 2] Register data from native space to fs_average space..."
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

            echo "Generate midthickness surface..."
            mris_convert $subject/surf/"$hemisphere".white $subject/surf/"$hemisphere".white.gii
            mris_convert $subject/surf/"$hemisphere".pial $subject/surf/"$hemisphere".pial.gii
            /home/ribeiro/Projects/deepRetinotopy_validation/deepRetinotopy_TheToolbox/utils/midthickness_surf.py --path $subject/surf/ --hemisphere $hemisphere # to edit after deepRetinotopy container is updated
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
                    metric_new="variance_explained"
                elif [ $metric == "x" ]; then
                    metric_new="x0"
                elif [ $metric == "y" ]; then
                    metric_new="y0"
                fi

                echo "Converting $metric data to .gii format..."

                # echo "Converting $metric data to .gii format..."
                mris_convert -c "$dataDir"/"$projectDir"/derivatives/prfanalyze-vista/$subject/ses-nyu3t01/"$hemisphere"."$metric".mgz $subject/surf/"$hemisphere".white \
                    "$dataDir"/"$projectDir"/derivatives/prfanalyze-vista/$subject/ses-nyu3t01/"$hemisphere"."$metric".gii \
                
                # Transform polar angle data before resampling
                if [ $metric == "angle" ]; then
                    python -c "import sys; sys.path.append('/home/ribeiro/Projects/deepRetinotopy_validation/'); from functions.preprocess import transform_angle_nyu; transform_angle_nyu('"$dataDir"/"$projectDir"/derivatives/prfanalyze-vista/$subject/ses-nyu3t01/"$hemisphere"."$metric".gii', '$hemisphere')"
                fi
                
                echo "Resampling native data to fsaverage space..."
                wb_command -metric-resample "$dataDir"/"$projectDir"/derivatives/prfanalyze-vista/$subject/ses-nyu3t01/"$hemisphere"."$metric"_transformed.gii \
                    "$dataDir"/"$projectDir"/derivatives/freesurfer/$subject/surf/"$hemisphere".sphere.reg.surf.gii /home/ribeiro/Projects/deepRetinotopy_validation/templates/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                    ADAP_BARY_AREA "$dataDir"/"$projectDir"/derivatives/freesurfer/$subject/surf/"$subject".fs_empirical_"$metric_new"_"$hemisphere".func.gii \
                    -area-surfs "$dataDir"/"$projectDir"/derivatives/freesurfer/$subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii            
            done
        done
    fi
done

# Run deepRetinotopy
echo "--------------------------------------------------------------------------------"
echo "[Step 3] Run deepRetinotopy..."
echo "--------------------------------------------------------------------------------"
deepRetinotopy -s "$dataDir"/"$projectDir"/derivatives/freesurfer -t ~/Projects/deepRetinotopy_validation/templates/ -d lnyu -m "polarAngle,eccentricity,pRFsize"  -g 'yes'