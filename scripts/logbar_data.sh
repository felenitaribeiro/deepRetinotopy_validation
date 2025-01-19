ml connectomeworkbench/1.5.0
ml freesurfer/7.3.2
ml deepretinotopy/1.0.6

cd ../datasets
# datalad install https://github.com/OpenNeuroDatasets/ds004698.git

# Download the required freesurfer data
cd ds004698/derivatives/freesurfer
for subject in `ls .`; do
    for hemisphere in lh rh; do
        datalad get $subject/surf/${hemisphere}.white
        datalad get $subject/surf/${hemisphere}.pial
        datalad get $subject/surf/${hemisphere}.sphere
        datalad get $subject/surf/${hemisphere}.sphere.reg
    done
done

# Download the required retinotopy data
cd ../prf-estimation
for subject in `ls .`; 
do
    if [ ${subject:0:3} != "sub" ]; then
        continue
    else
        for hemisphere in L R; 
        do
            datalad get $subject/prfs/"$subject"_ses-all_task-logbar_hemi-"$hemisphere"_space-fsnative_prf/"$subject"_ses-all_task-logbar_hemi-"$hemisphere"_space-fsnative_angle.mgz
            datalad get $subject/prfs/"$subject"_ses-all_task-logbar_hemi-"$hemisphere"_space-fsnative_prf/"$subject"_ses-all_task-logbar_hemi-"$hemisphere"_space-fsnative_eccen.mgz
            datalad get $subject/prfs/"$subject"_ses-all_task-logbar_hemi-"$hemisphere"_space-fsnative_prf/"$subject"_ses-all_task-logbar_hemi-"$hemisphere"_space-fsnative_sigma.mgz
        done
    fi
done

# Convert the retinotopy data to .gii format in the fs_32k space
cd ../freesurfer
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

            for metric in angle eccen sigma;
            do
                if [ $metric == "angle" ]; then
                    metric_new="polarAngle"
                elif [ $metric == "eccen" ]; then
                    metric_new="eccentricity"
                elif [ $metric == "sigma" ]; then
                    metric_new="pRFsize"
                fi

                echo "Converting $metric data to .gii format..."
                mris_convert -c ../prf-estimation/$subject/prfs/"$subject"_ses-all_task-logbar_hemi-"$hemi"_space-fsnative_prf/"$subject"_ses-all_task-logbar_hemi-"$hemi"_space-fsnative_"$metric".mgz $subject/surf/"$hemisphere".white \
                ../prf-estimation/$subject/prfs/"$subject"_ses-all_task-logbar_hemi-"$hemi"_space-fsnative_prf/"$subject"_ses-all_task-logbar_hemi-"$hemi"_space-fsnative_"$metric".gii \

                echo "Resampling native data to fsaverage space..."
                wb_command -metric-resample ../prf-estimation/$subject/prfs/"$subject"_ses-all_task-logbar_hemi-"$hemi"_space-fsnative_prf/"$subject"_ses-all_task-logbar_hemi-"$hemi"_space-fsnative_"$metric".gii \
                    $subject/surf/"$hemisphere".sphere.reg.surf.gii ../../../../templates/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                    ADAP_BARY_AREA $subject/surf/"$subject".empirical_"$metric_new"."$hemisphere".32k_fs_LR.func.gii \
                    -area-surfs $subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii            
            done
        done
    fi
done

# # Run deepRetinotopy
# deepRetinotopy -s ./ -t ../../../../templates/ -d logbar -m "polarAngle,eccentricity,pRFsize" -g 'yes'
# signMaps -s ./ -t ../../../../templates/ -d logbar 
# # Convert polar angle range 

# # TODO