ml connectomeworkbench/1.5.0
ml freesurfer/7.3.2
ml deepretinotopy/1.0.8

dataDir=/BULK/ribeiro
projectURL=https://github.com/OpenNeuroDatasets/ds003787.git
projectDir=${projectURL:37:-4} # after slash before .git
cd $dataDir
echo `pwd $dataDir`

datalad install $projectURL

# # Download the required freesurfer data
# cd "$dataDir"/"$projectDir"/derivatives/freesurfer
# echo `pwd $dataDir/$projectDir/derivatives/freesurfer`
# for subject in `ls .`; do
#     if [ ${subject:0:3} != "sub" ]; then
#         continue
#     else
#         for hemisphere in lh rh; do
#             datalad get $subject/surf/${hemisphere}.white
#             datalad get $subject/surf/${hemisphere}.pial
#             datalad get $subject/surf/${hemisphere}.sphere
#             datalad get $subject/surf/${hemisphere}.sphere.reg
#             datalad get $subject/surf/${hemisphere}.thickness
#         done
#     fi
# done
# echo "Required freesurfer data downloaded"

# # Download the required retinotopy data
# cd "$dataDir"/"$projectDir"/derivatives/prfanalyze-vista
# for subject in `ls .`; 
# do
#     if [ ${subject:0:3} != "sub" ]; then
#         continue
#     else
#         for hemisphere in lh rh; 
#         do 
#             datalad get $subject/ses-nyu3t01/"$hemisphere".angle.mgz
#             datalad get $subject/ses-nyu3t01/"$hemisphere".eccen.mgz
#             datalad get $subject/ses-nyu3t01/"$hemisphere".sigma.mgz
#             datalad get $subject/ses-nyu3t01/"$hemisphere".vexpl.mgz
#             datalad get $subject/ses-nyu3t01/"$hemisphere".x.mgz
#             datalad get $subject/ses-nyu3t01/"$hemisphere".y.mgz
#         done
#     fi
# done
# echo "Required retinotopy data downloaded"

# Download the required retinotopy data
cd $dataDir/"$projectDir"/derivatives/stimulus_apertures/
for subject in `ls .`; 
do  
    datalad get $subject/* 
    datalad get $dataDir/"$projectDir"/derivatives/fmriprep/$subject/ses-nyu3t01/func/*
done
echo "Required stimulus apertures downloaded and fMRI data downloaded"


# Convert the retinotopy data to .gii format in the fs_32k space
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
            ../../../../functions/midthickness_surf.py --path $subject/surf/ --hemisphere $hemisphere # to edit after deepRetinotopy container is updated
            mris_convert $subject/surf/"$hemisphere".graymid.gii $subject/surf/"$hemisphere".graymid
            mris_curvature -w $subject/surf/"$hemisphere".graymid

            echo "Resampling native surface to fs_LR space..."
            wb_shortcuts -freesurfer-resample-prep $subject/surf/"$hemisphere".white $subject/surf/"$hemisphere".pial \
            $subject/surf/"$hemisphere".sphere.reg ../../../../templates/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
            $subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii \
            $subject/surf/"$hemisphere".sphere.reg.surf.gii

            for metric in angle eccen sigma vexpl x y;
            do
                if [ $metric == "angle" ]; then
                    metric_new="polarAngle"
                elif [ $metric == "eccen" ]; then
                    metric_new="eccentricity"
                elif [ $metric == "sigma" ]; then
                    metric_new="pRFsize"
                elif [ $metric == "vexpl" ]; then
                    metric_new="vexpl"
                elif [ $metric == "x" ]; then
                    metric_new="x0"
                elif [ $metric == "y" ]; then
                    metric_new="y0"
                fi

                echo "Converting $metric data to .gii format..."
                if [ $metric == "angle " ] && [ $hemisphere == "lh" ]; then
                    # transform angle range
                fi
                mris_convert -c $dataDir/ds003787/derivatives/prfanalyze-vista/$subject/ses-nyu3t01/"$hemisphere"."$metric".mgz $subject/surf/"$hemisphere".white \
                $dataDir/ds003787/derivatives/prfanalyze-vista/$subject/ses-nyu3t01/"$hemisphere"."$metric".gii \

                echo "Resampling native data to fsaverage space..."
                wb_command -metric-resample $dataDir/ds003787/derivatives/prfanalyze-vista/$subject/ses-nyu3t01/"$hemisphere"."$metric".gii \
                    $dataDir/ds003787/derivatives/freesurfer/$subject/surf/"$hemisphere".sphere.reg.surf.gii /home/ribeiro/Projects/deepRetinotopy_validation/templates/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                    ADAP_BARY_AREA $dataDir/ds003787/derivatives/freesurfer/$subject/surf/"$subject".empirical_"$metric_new"."$hemisphere".32k_fs_LR.func.gii \
                    -area-surfs $dataDir/ds003787/derivatives/freesurfer/$subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii            
            done

            for run in run-01 run-02 run-03 run-04 run-05 run-06 run-07 run-08 run-09 run-10 run-11 run-12;
            do 
                if [ ! -f '$dataDir/ds003787/derivatives/fmriprep/$subject/ses-nyu3t01/func/"$subject"_ses-nyu3t01_task-prf_acq-PA_"$run"_space-fsnative_hemi-"$hemi".func.mgz' ]; then
                    continue
                else
                    mris_convert -c $dataDir/ds003787/derivatives/fmriprep/$subject/ses-nyu3t01/func/"$subject"_ses-nyu3t01_task-prf_acq-PA_"$run"_space-fsnative_hemi-"$hemi".func.mgz $subject/surf/"$hemisphere".white \
                    $dataDir/ds003787/derivatives/fmriprep/$subject/ses-nyu3t01/func/"$subject"_ses-nyu3t01_task-prf_acq-PA_"$run"_space-fsnative_hemi-"$hemi".func.gii \

                    wb_command -metric-resample $dataDir/ds003787/derivatives/fmriprep/$subject/ses-nyu3t01/func/"$subject"_ses-nyu3t01_task-prf_acq-PA_"$run"_space-fsnative_hemi-"$hemi".func.gii \
                        $dataDir/ds003787/derivatives/freesurfer/$subject/surf/"$hemisphere".sphere.reg.surf.gii /home/ribeiro/Projects/deepRetinotopy_validation/templates/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                        ADAP_BARY_AREA $dataDir/ds003787/derivatives/freesurfer/$subject/surf/"$subject"_ses-nyu3t01_task-prf_acq-PA_"$run"_space-fsnative_hemi-"$hemi".32k_fs_LR.func.gii \
                        -area-surfs $dataDir/ds003787/derivatives/freesurfer/$subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii
                fi
            done
        done
    fi
done

# # # Run deepRetinotopy
# # deepRetinotopy -s ./ -t ../../../../templates/ -d logbar -m "polarAngle" -g 'yes'

# # # # Convert polar angle range 

# # # # TODO