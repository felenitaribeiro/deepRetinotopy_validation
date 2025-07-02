#!/bin/bash
module use /sw/local/rocky8/noarch/neuro/software/neurocommand/local/containers/modules/
export APPTAINER_BINDPATH=/scratch,/QRISdata

ml connectomeworkbench/1.5.0


dataDir=../samsrf_all/

# Data processing
echo "--------------------------------------------------------------------------------"
echo "[Step 1] Register data from native space to fs_average space..."
echo "--------------------------------------------------------------------------------"
cd "$dataDir"
echo `pwd`
for subject in `ls .`;
do
    for hemisphere in lh rh;
    do
        if [ $hemisphere == "lh" ]; then
            hemi="L"
        else
            hemi="R"
        fi

        for metric in pol ecc sigma r^2;
        do  
            if [ $metric == "pol" ]; then
                metric_new="polarAngle"
            elif [ $metric == "ecc" ]; then
                metric_new="eccentricity"
            elif [ $metric == "sigma" ]; then
                metric_new="pRFsize"
            elif [ $metric == "r^2" ]; then
                metric_new="variance_explained"
            fi

            if [ $metric == "pol" ]; then
                if [ $hemisphere == 'lh' ]; then
                    echo "Resampling $metric data..."
                    wb_command -metric-resample "$dataDir"/$subject/"$hemisphere"_"$subject"_"$metric".gii \
                        /scratch/project/recyle_dl/deepRetinotopy_TheToolbox/HCP/freesurfer/$subject/surf/"$hemisphere".sphere.reg.surf.gii /scratch/project/recyle_dl/deepRetinotopy_TheToolbox/templates/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                        ADAP_BARY_AREA "$dataDir"/$subject/"$subject".fs_empirical_samsrf_"$metric_new"_"$hemisphere".func.gii \
                        -area-surfs /scratch/project/recyle_dl/deepRetinotopy_TheToolbox/HCP/freesurfer/$subject/surf/"$hemisphere".midthickness.surf.gii /scratch/project/recyle_dl/deepRetinotopy_TheToolbox/HCP/freesurfer/$subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii   
                    python -c "import sys; sys.path.append('/scratch/project/recyle_dl/deepRetinotopy_validation'); \
                            from functions.preprocess import transform_angle; \
                            transform_angle('"$dataDir"/$subject/"$subject".fs_empirical_samsrf_"$metric_new"_"$hemisphere".func.gii', '$hemisphere', left_hemi_shift=False)"
                                    echo "Resampling $metric data..." 
                    mv "$dataDir"/$subject/"$subject".fs_empirical_samsrf_"$metric_new"_"$hemisphere".func_transformed.gii "$dataDir"/$subject/"$subject".fs_empirical_samsrf_"$metric_new"_"$hemisphere".func.gii 
                else
                    python -c "import sys; sys.path.append('/scratch/project/recyle_dl/deepRetinotopy_validation'); \
                        from functions.preprocess import transform_angle; \
                        transform_angle('"$dataDir"/"$subject"/"$hemisphere"_"$subject"_"$metric".gii', '$hemisphere', left_hemi_shift=False)"
                                echo "Resampling $metric data..."
                    echo "Resampling $metric data..."
                    wb_command -metric-resample "$dataDir"/$subject/"$hemisphere"_"$subject"_"$metric"_transformed.gii \
                        /scratch/project/recyle_dl/deepRetinotopy_TheToolbox/HCP/freesurfer/$subject/surf/"$hemisphere".sphere.reg.surf.gii /scratch/project/recyle_dl/deepRetinotopy_TheToolbox/templates/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                        ADAP_BARY_AREA "$dataDir"/$subject/"$subject".fs_empirical_samsrf_"$metric_new"_"$hemisphere".func.gii \
                        -area-surfs /scratch/project/recyle_dl/deepRetinotopy_TheToolbox/HCP/freesurfer/$subject/surf/"$hemisphere".midthickness.surf.gii /scratch/project/recyle_dl/deepRetinotopy_TheToolbox/HCP/freesurfer/$subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii   
                fi
            else
                echo "Resampling $metric data..."
                wb_command -metric-resample "$dataDir"/$subject/"$hemisphere"_"$subject"_"$metric".gii \
                    /scratch/project/recyle_dl/deepRetinotopy_TheToolbox/HCP/freesurfer/$subject/surf/"$hemisphere".sphere.reg.surf.gii /scratch/project/recyle_dl/deepRetinotopy_TheToolbox/templates/fs_LR-deformed_to-fsaverage."$hemi".sphere.32k_fs_LR.surf.gii \
                    ADAP_BARY_AREA "$dataDir"/$subject/"$subject".fs_empirical_samsrf_"$metric_new"_"$hemisphere".func.gii \
                    -area-surfs /scratch/project/recyle_dl/deepRetinotopy_TheToolbox/HCP/freesurfer/$subject/surf/"$hemisphere".midthickness.surf.gii /scratch/project/recyle_dl/deepRetinotopy_TheToolbox/HCP/freesurfer/$subject/surf/"$subject"."$hemisphere".midthickness.32k_fs_LR.surf.gii   
            fi
        done
    done
done
