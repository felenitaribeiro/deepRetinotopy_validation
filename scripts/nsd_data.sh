ml connectomeworkbench/1.5.0
ml freesurfer/7.3.2
ml deepretinotopy/1.0.8

dataDir=/BULK/LABDATA/NSD
cd $dataDir


# Run deepRetinotopy
echo "--------------------------------------------------------------------------------"
echo "[Step 1] Run deepRetinotopy..."
echo "--------------------------------------------------------------------------------"
deepRetinotopy -s "$dataDir"/freesurfer/ -t /home/ribeiro/Projects/deepRetinotopy_validation/templates/ -d nsd -m 'polarAngle,prfeccentricitytricity,pRFsize' 

# Data processing
echo "--------------------------------------------------------------------------------"
echo "[Step 2] Register data from native space to fs_average space..."
echo "--------------------------------------------------------------------------------"
cd "$dataDir"/freesurfer/
echo `pwd`
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

            for metric in prfangle prfeccentricity prfsize prfR2 prfexponent;
            do  
                if [ $metric == "prfangle" ]; then
                    metric_new="polarAngle"
                elif [ $metric == "prfeccentricity" ]; then
                    metric_new="eccentricity"
                elif [ $metric == "prfsize" ]; then
                    metric_new="pRFsize"
                elif [ $metric == "prfR2" ]; then
                    metric_new="variance_explained"
                elif [ $metric == "prfexponent" ]; then
                    metric_new="prf_exponent"
                fi
                # echo "mris_convert "$dataDir"/freesurfer/$subject/label/"$hemisphere"."$metric".mgz "$dataDir"/freesurfer/$subject/label/"$hemisphere"."$metric".gii"
                mris_convert -c "$dataDir"/freesurfer/"$subject"/label/"$hemisphere"."$metric".mgz "$dataDir"/freesurfer/"$subject"/surf/"$hemisphere".white \
                    "$dataDir"/freesurfer/"$subject"/label/"$hemisphere"."$metric".gii
                cp "$dataDir"/freesurfer/"$subject"/label/"$hemisphere"."$metric".gii "$dataDir"/freesurfer/"$subject"/surf/"$hemisphere"."$metric".gii

                if [ $metric == "prfangle" ] && [ $hemisphere == "lh" ]; then
                    python -c "import sys; sys.path.append('/home/ribeiro/Projects/deepRetinotopy_validation/'); \
                            from functions.preprocess import transform_angle_lh_nsd; \
                            transform_angle_lh_nsd('"$dataDir"/freesurfer/"$subject"/label/"$hemisphere"."$metric".gii')"

                    echo "Resampling $metric data..."
                    wb_command -metric-resample "$dataDir"/freesurfer/$subject/label/"$hemisphere"."$metric"_transformed.gii \
                        $subject/surf/"$hemisphere".sphere.reg.surf.gii /home/ribeiro/Projects/deepRetinotopy_validation/templates/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                        ADAP_BARY_AREA $subject/deepRetinotopy/"$subject".fs_empirical_"$metric_new"_"$hemisphere".func.gii \
                        -area-surfs $subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii   
                        
                else
                    echo "Resampling $metric data..."
                    wb_command -metric-resample "$dataDir"/freesurfer/$subject/label/"$hemisphere"."$metric".gii \
                        $subject/surf/"$hemisphere".sphere.reg.surf.gii /home/ribeiro/Projects/deepRetinotopy_validation/templates/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                        ADAP_BARY_AREA $subject/deepRetinotopy/"$subject".fs_empirical_"$metric_new"_"$hemisphere".func.gii \
                        -area-surfs $subject/surf/"$hemisphere".midthickness.surf.gii $subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii   
                fi
            done
        done
    fi
done
