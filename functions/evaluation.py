import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
import os.path as osp
import seaborn as sns
import pandas as pd
import os
import sys
import scipy

sys.path.append('..')

from functions.datasets import *
from deepRetinotopy_TheToolbox.utils.metrics import smallest_angle
from deepRetinotopy_TheToolbox.utils.rois import ROI_WangParcelsPlusFovea as roi
from deepRetinotopy_TheToolbox.utils.rois import ROIs_WangParcels as roi_parcel
from sklearn.metrics import jaccard_score
from nilearn.glm.first_level.hemodynamic_models import _gamma_difference_hrf
from scipy.stats import multivariate_normal
from astropy.stats import circcorrcoef
from astropy import units as u

def transform_polarangle(data):
        """Transform polar angle values from the left hemisphere to the range from 0-90 degrees 
        (UVF) and from 360-270 degrees (LVF).

        Args:
            data (numpy array): Polar angle values

        Returns:
            data (numpy array): Transformed polar angle values
        """
        mask = data < 0
        subtract = data > 180
        add = data < 180
        data[subtract] = data[subtract] - 180
        data[add] = data[add] + 180
        data[mask] = -1
        return data

def roi_earlyvisualcortex(list_of_labels):
    """Mask for the selection of the region of interest in the surface
    template.

    Args:
        list_of_labels (list): list with the file name (.mat) containing the
            region of interest (from both L and R hemispheres)

    Returns:
        final_mask_L (numpy array): Mask of the region of interest from left
            hemisphere (32492,)

        final_mask_R (numpy array): Mask of the region of interest from
            right hemisphere (32492,)

        index_L_mask (list): Indices of the non-zero elements from
            final_mask_L (number of nonzero elements,)

        index_R_mask (list): Indices of the non-zero elements from
            final_mask_R (number of nonzero elements,)

    """

    # Defining number of nodes
    number_cortical_nodes = int(64984)
    number_hemi_nodes = int(number_cortical_nodes / 2)

    list_primary_visual_areas = np.zeros([len(list_of_labels), 64984])
    for i in range(len(list_of_labels)):
        list_primary_visual_areas[i] = np.reshape(scipy.io.loadmat(
            osp.join(osp.dirname(osp.realpath(__file__)),
                     'rois/EarlyVisualCortex',
                     list_of_labels[i] + '_V1-3.mat'))[list_of_labels[i]][0:64984],
                                                  (-1))

    final_mask_L = np.sum(list_primary_visual_areas, axis=0)[
                   0:number_hemi_nodes]
    final_mask_R = np.sum(list_primary_visual_areas, axis=0)[
                   number_hemi_nodes:number_cortical_nodes]

    index_L_mask = [i for i, j in enumerate(final_mask_L) if j == 1]
    index_R_mask = [i for i, j in enumerate(final_mask_R) if j == 1]

    return final_mask_L, final_mask_R, index_L_mask, index_R_mask

def discretize(ind_map, retinotopic_map, hemisphere = 'lh'):
    '''
    Binarize the retinotopic map according to the type of map.
    
    Args:
        ind_map (numpy array): Retinotopic map
        retinotopic_map (str): Type of retinotopic map ('polarAngle', 'eccentricity' or 'pRFsize')
        hemisphere (str): Hemisphere ('lh' or 'rh')

    Returns:
        ind_map (numpy array): Binarized retinotopic map
    '''
    if retinotopic_map == 'polarAngle':
        if hemisphere == 'lh':
            ind_map[(ind_map >= 0) & (ind_map <= 45)] = 0
            ind_map[(ind_map > 45) & (ind_map <= 180)] = 90
            ind_map[(ind_map >= 315) & (ind_map <= 360)] = 360
            ind_map[(ind_map > 180) & (ind_map < 315)] = 270
        else:
            ind_map[(ind_map >= 0) & (ind_map <= 135)] = 90
            ind_map[(ind_map > 135) & (ind_map <= 180)] = 130
            ind_map[(ind_map > 180) & (ind_map <= 225)] = 220
            ind_map[(ind_map > 225) & (ind_map <= 360)] = 270
    elif retinotopic_map == 'eccentricity':
        ind_map[(ind_map >= 0) & (ind_map <= 2)] = 0
        ind_map[(ind_map > 2) & (ind_map <= 4)] = 2
        ind_map[(ind_map > 4) & (ind_map <= 6)] = 4
        ind_map[(ind_map > 6)] = 6
    elif retinotopic_map == 'pRFsize':
        ind_map[(ind_map > 4)] = 9
        ind_map[(ind_map >= 3.5) & (ind_map <= 4)] = 8
        ind_map[(ind_map >= 3) & (ind_map < 3.5)] = 7
        ind_map[(ind_map >= 2.5) & (ind_map < 3)] = 6
        ind_map[(ind_map >= 2) & (ind_map < 2.5)] = 5
        ind_map[(ind_map >= 1.5) & (ind_map < 2)] = 4
        ind_map[(ind_map >= 1) & (ind_map < 1.5)] = 3
        ind_map[(ind_map >= .5) & (ind_map < 1)] = 2
        ind_map[(ind_map >= .25) & (ind_map < .5)] = 1
        ind_map[(ind_map >= 0) & (ind_map < .25)] = 0
    return ind_map

def calculate_overlap(map_1, map_2, retinotopic_map, angle = 'rad', hemisphere = 'lh'):
    """ Calculate the overlap (Jaccard index) between two retinotopic maps.

    Args:
        map_1 (numpy array): Predicted (or empirical) retinotopic map from subject 1
        map_2 (numpy array): Predicted (or empirical) retinotopic map from subject 2 (or 1)
        retinotopic_map (str): Type of retinotopic map ('polarAngle', 'eccentricity' or 'pRFsize')
        mask (numpy array): Mask of the region of interest
        angle (str): Units of the angle (rad or original)
        hemisphere (str): Hemisphere ('lh' or 'rh')
    
    Returns:
        score (float): Overlap between the two retinotopic maps
    """
    if angle == 'rad':
        map_1 = map_1 * (180/np.pi)
        map_2 = map_2 * (180/np.pi)
    
    map_1 = discretize(map_1, retinotopic_map, hemisphere)
    map_2 = discretize(map_2, retinotopic_map, hemisphere)

    jaccard_1 = jaccard_score(map_1, map_2, average='weighted')
    jaccard_2 = jaccard_score(map_2, map_1, average='weighted')
    score = (jaccard_1 + jaccard_2) / 2
    return score

def create_mask(final_mask_L_ROI, final_mask_R_ROI, final_mask_L, final_mask_R, hemisphere):
    """Create ROI and visual area masks based on hemisphere."""
    ROI_masked = np.zeros((32492, 1))
    visualarea = np.zeros((32492, 1))
    if hemisphere == 'lh':
        ROI_masked[final_mask_L_ROI == 1] = 1
        visualarea[final_mask_L == 1] = 1
    else:
        ROI_masked[final_mask_R_ROI == 1] = 1
        visualarea[final_mask_R == 1] = 1
    mask = ROI_masked + visualarea
    return ROI_masked, mask == 2

def return_list_of_subs(dataset_name, pool_of_participants = 'all'):
    """Return the list of subject IDs for the dataset.

    Args:
        dataset_name (str): Name of the dataset ('hcp', 'logbar' or 'nyu')

    Returns:
        list_of_subs (list): List of subject IDs
    """

    if dataset_name == 'hcp':
        list_of_subs = ['680957', '191841', '617748', '725751', '198653',
                         '191336', '572045', '601127', '644246', '157336']
    else:
        if dataset_name == 'logbar':
            path = '/BULK/LABDATA/openneuro/ds004698/derivatives/freesurfer/'
            list_of_subs = os.listdir(path)
        elif dataset_name == 'nyu':
            path = '/BULK/LABDATA/openneuro/ds003787/derivatives/freesurfer/'
            list_of_subs = os.listdir(path)
        elif dataset_name == 'nsd':
            path = '/BULK/LABDATA/NSD/freesurfer/'
            list_of_subs = os.listdir(path)
        elif dataset_name == 'HCP':
            path = '/home/ribeiro/Projects/deepRetinotopy_validation/HCP/freesurfer/'
            list_of_subs = os.listdir(path)
        elif dataset_name == 'stanford':
            if pool_of_participants == 'all':
                path = '/BULK/LABDATA/openneuro/stanford-data/ds004440/derivatives/freesurfer/'
            elif pool_of_participants == 'children':
                path = '/BULK/LABDATA/openneuro/stanford-data/ds004440/derivatives/prfanalyze-vista/children/'
            elif pool_of_participants == 'adults':
                path = '/BULK/LABDATA/openneuro/stanford-data/ds004440/derivatives/prfanalyze-vista/adults/'
            list_of_subs = os.listdir(path)
        
    return list_of_subs

def apply_threshold(data, threshold):
    """Apply threshold to predicted and empirical maps."""
    if threshold is not None:
        data.predicted_map = data.predicted_map[data.variance_explained > threshold]
        data.empirical_map = data.empirical_map[data.variance_explained > threshold]

def process_subjects(list_of_sub_ids, path, dataset_name, hemisphere, mask_final, retinotopic_map, experiment, threshold):
    """Process subjects and return predicted and empirical maps."""
    predicted_map_hemi = []
    empirical_map_hemi = []
    for sub_id in list_of_sub_ids:
        if dataset_name == 'logbar':
            print(retinotopic_map)
            data = RetinotopyData_logbar(path, sub_id, hemisphere, retinotopic_map, experiment=experiment)
        else:
            data = RetinotopyData(path, sub_id, hemisphere, retinotopic_map)
        data.apply_mask_to_maps(mask_final)
        apply_threshold(data, threshold)
        if retinotopic_map == 'polarAngle' and hemisphere == 'lh':
            data.apply_transform_polarangle()
        if dataset_name == 'logbar':
            second_mask = data.empirical_map > 0
            data.predicted_map = data.predicted_map[second_mask]
            data.empirical_map = data.empirical_map[second_mask]
        predicted_map_hemi += list(data.predicted_map.flatten())
        empirical_map_hemi += list(data.empirical_map.flatten())
    return predicted_map_hemi, empirical_map_hemi

def remove_outliers(predicted_map, empirical_map, retinotopic_map, dataset_name):
    """Remove outliers for eccentricity and pRFsize maps."""
    if retinotopic_map == 'eccentricity' and dataset_name in ['logbar', 'stanford']:
        predicted_map = np.clip(predicted_map, 0, 10.5)
        empirical_map = np.clip(empirical_map, 0, 10.5)
    elif retinotopic_map == 'pRFsize' and dataset_name == 'logbar':
        predicted_map = np.clip(predicted_map, 0, 3.5)
        empirical_map = np.clip(empirical_map, 0, 3.5)
    else:
        predicted_map = np.clip(predicted_map, 0, None)
    return predicted_map, empirical_map

def metric_model_selection(path, retinotopic_map, hemisphere, retinotopic_mapping = 'continuous', threshold = 10):
    """Calculate the inter-individual variability in predicted maps and the error
    between predicted and empirical maps for different models.
    
    Args:
        path (str): Path to the Freesurfer directory (both empirical and predicted maps should be 
            located at the deepRetinotopy folder within the individual's directory)
        retinotopic_map (str): Type of retinotopic map ('polarAngle', 'eccentricity' or 'pRFsize')
        hemisphere (str): Hemisphere ('lh' or 'rh')
        retinotopic_mapping (str): Type of retinotopic mapping ('continuous' or 'discrete')
        threshold (float): Threshold for the explained variance

    Returns:
        df (pandas dataframe): Dataframe containing the inter-individual variability in predicted maps 
            and the error.
        plt (matplotlib.pyplot): Plot of the inter-individual variability in predicted maps and the error. 
            PNG file saved at the output folder.
    """
    mean_ind_variability = []
    error = []

    # Development dataset
    dev_set = ['186949', '169747', '826353', '825048', '671855',
                    '751550', '318637', '131722', '137128', '706040'] 
    
    # Region of interest 
    ## Region of interest used for training
    final_mask_L_ROI, final_mask_R_ROI, index_L_mask, index_R_mask = roi(['ROI'])
    ## Early visual cortex
    final_mask_L, final_mask_R, index_L_mask, index_R_mask = roi_earlyvisualcortex(['ROI'])
    
    ROI_masked, mask = create_mask(final_mask_L_ROI, final_mask_R_ROI, final_mask_L, final_mask_R, hemisphere)
    mask = mask[ROI_masked == 1]
    for model in range(5):
        theta_withinsubj = []
        theta_acrosssubj_pred = []
        print('Model ' + str(model + 1) + ' of 5')
        for j in range(len(dev_set)):
            theta_pred_across_temp = []
            for i in range(len(dev_set)):
                # Error
                if i == j:
                    # Loading maps
                    data = RetinotopyData(path, dev_set[i], hemisphere, retinotopic_map, model_index=str(model + 1))
                    data.apply_mask_to_maps(ROI_masked)
                    data.apply_mask_to_maps(mask)

                    # # Transform polar angle values from lh
                    if retinotopic_map == 'polarAngle' and hemisphere == 'lh':
                        data.apply_transform_polarangle()

                    assert np.min(data.predicted_map) >= 0.
                    assert np.min(data.empirical_map) >= 0.
                    
                    # Transforming to radians
                    if retinotopic_map == 'polarAngle' or retinotopic_map == 'eccentricity':
                        data.convert_to_radian()

                    # Calculating error
                    if retinotopic_mapping == 'continuous':
                        if retinotopic_map == 'polarAngle' or retinotopic_map == 'eccentricity':
                            theta = smallest_angle(data.predicted_map, data.empirical_map)
                            if threshold != None:
                                theta = theta[data.variance_explained > threshold]
                        elif retinotopic_map == 'pRFsize':
                            theta = np.abs(data.predicted_map - data.empirical_map)
                            if threshold != None:
                                theta = theta[data.variance_explained > threshold]
                        theta_withinsubj.append(np.mean(theta))
                    elif retinotopic_mapping == 'discrete':
                        if retinotopic_map == 'polarAngle' or retinotopic_map == 'eccentricity':
                            overlap = calculate_overlap(data.predicted_map, data.empirical_map, retinotopic_map, angle = 'rad', hemisphere = hemisphere)
                        else:
                            overlap = calculate_overlap(data.predicted_map, data.empirical_map, retinotopic_map, angle = 'original', hemisphere = hemisphere)
                        theta_withinsubj.append(overlap)

                # Inter-individual variability in predicted maps
                if i != j:
                    # Loading maps
                    data_1 = RetinotopyData(path, dev_set[i], hemisphere, retinotopic_map, model_index=str(model + 1))
                    data_1.apply_mask_to_maps(ROI_masked)
                    data_1.apply_mask_to_maps(mask)

                    data_2 = RetinotopyData(path, dev_set[j], hemisphere, retinotopic_map, model_index=str(model + 1))
                    data_2.apply_mask_to_maps(ROI_masked)
                    data_2.apply_mask_to_maps(mask)
                    
                    # Transform polar angle values from lh
                    if retinotopic_map == 'polarAngle' and hemisphere == 'lh':
                        data_1.apply_transform_polarangle()
                        data_2.apply_transform_polarangle()
                    assert np.min(data_1.predicted_map) >= 0.
                    assert np.min(data_2.predicted_map) >= 0.
                    
                    # Transforming to radians
                    if retinotopic_map == 'polarAngle' or retinotopic_map == 'eccentricity':
                        data_1.convert_to_radian()
                        data_2.convert_to_radian()
                    
                    # Calculating inter-individual variability
                    if retinotopic_mapping == 'continuous':
                        if retinotopic_map == 'polarAngle' or retinotopic_map == 'eccentricity':
                            theta_pred = smallest_angle(data_1.predicted_map,
                                                    data_2.predicted_map)
                        elif retinotopic_map == 'pRFsize':
                            theta_pred = np.abs(data_1.predicted_map - data_2.predicted_map)
                        theta_pred_across_temp.append(
                            np.mean(theta_pred))
                    elif retinotopic_mapping == 'discrete':
                        if retinotopic_map == 'polarAngle' or retinotopic_map == 'eccentricity':
                            overlap = calculate_overlap(data_1.predicted_map, data_2.predicted_map, retinotopic_map, angle = 'rad', hemisphere = hemisphere)
                        else:
                            overlap = calculate_overlap(data_1.predicted_map, data_2.predicted_map, retinotopic_map, angle = 'original', hemisphere = hemisphere)
                        theta_pred_across_temp.append(overlap)
            theta_acrosssubj_pred.append(theta_pred_across_temp)

        mean_theta_withinsubj = np.array(theta_withinsubj)
        mean_theta_acrosssubj_pred = np.mean(np.array(theta_acrosssubj_pred), axis=1)

        error.append(mean_theta_withinsubj)
        mean_ind_variability.append(mean_theta_acrosssubj_pred)

    error = np.reshape(np.array(error), (5, -1))
    mean_ind_variability = np.reshape(np.array(mean_ind_variability), (5, -1))

    data_ind_variability = [[mean_ind_variability[i], 
                             len(mean_ind_variability[i]) * ['Seed ' + str(i+1)],
                            len(mean_ind_variability[i]) * ['Individual variability']] for i in range(5)]
    data_error = [[error[i], 
                    len(mean_ind_variability[i]) * ['Seed ' + str(i+1)],
                    len(mean_ind_variability[i]) * ['Error']] for i in range(5)]
    data = np.concatenate(data_ind_variability + data_error, axis = 1)

    df = pd.DataFrame(columns=['$\Delta$$\t\Theta$', 'Seed', 'Metric'],
                      data=data.T)
    df['$\Delta$$\t\Theta$'] = df['$\Delta$$\t\Theta$'].astype(float)

    # Generate plot
    if os.path.isdir('../output/model_selection/') == False:
        os.makedirs('../output/model_selection/')
        
    sns.set_style("whitegrid")
    fig = plt.figure()
    ax = sns.catplot(y='$\Delta$$\t\Theta$', x='Seed',
                       col='Metric', data=df, palette="flare", kind="swarm", hue='Seed')
 
    ax.set_titles("{col_name}")
    if retinotopic_mapping == 'discrete':
        ax.set_axis_labels("Seed", "Mean Jaccard index")
    else:
        ax.set_axis_labels("Seed", "Mean $\Delta$$\t\Theta$")
    if str(retinotopic_map) == 'eccentricity':
        if retinotopic_mapping == 'continuous':
            plt.ylim([0, 2])
        else:
            plt.ylim([0, 1])
    elif str(retinotopic_map) == 'pRFsize':
        if retinotopic_mapping == 'continuous':
            plt.ylim([0, 1])
        else:
            plt.ylim([0, 1])
    else:
        if retinotopic_mapping == 'continuous':
            plt.ylim([0, 40])
        else:
            plt.ylim([0, 1])

    fig.suptitle('Early visual areas')
    if threshold != None and retinotopic_mapping == 'continuous':
         plt.savefig('../output/model_selection/ModelSelection_EarlyVisualAreas_' + retinotopic_map + '_' + retinotopic_mapping + '_' + hemisphere + '_' + str(threshold) + '.pdf')
    else:
         plt.savefig('../output/model_selection/ModelSelection_EarlyVisualAreas_' + retinotopic_map + '_' + retinotopic_mapping + '_' + hemisphere + '.pdf')
    plt.show()
    return df

def predicted_vs_empirical(path, dataset_name, retinotopic_maps, hemispheres, 
                                    threshold = None, experiment = None, 
                                    region_of_interest = 'all', pool_of_participants = 'all', 
                                    plot_type = 'hexbin'):
    """Generate hexbin plots of predicted and empirical maps for different models.
    
    Args:
        path (str): Path to the Freesurfer directory (both empirical and predicted maps should be 
            located at the deepRetinotopy folder within the individual's directory)
        dataset_name (str): Name of the dataset ('hcp', 'logbar', 'nyu', 'nsd', 'HCP' or 'stanford')
        retinotopic_maps (list): List of retinotopic maps ('polarAngle', 'eccentricity' or 'pRFsize')
        hemispheres (str): Hemisphere ('lh', 'rh' or 'both')
        threshold (float): Threshold for the explained variance
        experiment (str): Experiment name (only for the logbar dataset)
        region_of_interest (str): Region of interest ('all', 'earlyvisualcortex', 'V1', 'V2' or 'V3')
        pool_of_participants (str): Pool of participants ('all', 'children' or 'adults'; only for the stanford dataset)
        plot_type (str): Type of plot ('hexbin' or 'scatter')

    Returns:
        plt (matplotlib.pyplot): Scatter plot of the predicted vs empirical maps. 
            PNG file saved at the output folder.
    """

    # Use data from both hemispheres or not
    if hemispheres == 'both':
        hemispheres = ['lh', 'rh']
    elif hemispheres == 'lh' or hemispheres == 'rh':
        hemispheres = [hemispheres]

    # Region of interest
    print('Region of interest: ' + region_of_interest)
    ## Region of interest used for training
    final_mask_L_ROI, final_mask_R_ROI, index_L_mask, index_R_mask = roi(['ROI'])
    ## Early visual cortex
    final_mask_L, final_mask_R, index_L_mask, index_R_mask = roi_earlyvisualcortex(['ROI'])
    if region_of_interest != 'all' and region_of_interest != 'earlyvisualcortex':
        if region_of_interest == 'V1':
            areas = ['V1v', 'V1d']
        elif region_of_interest == 'V2':
            areas = ['V2v', 'V2d']
        elif region_of_interest == 'V3':
            areas = ['V3v', 'V3d']
        final_mask_L, final_mask_R, index_L_mask, index_R_mask = roi_parcel(areas)

    # Make output directory
    if os.path.isdir('../output/model_evaluation/' + dataset_name) == False:
        os.makedirs('../output/model_evaluation/' + dataset_name)

    # List of subjects
    list_of_sub_ids = return_list_of_subs(dataset_name, pool_of_participants)

    for retinotopic_map in retinotopic_maps:
        predicted_map = []
        empirical_map = []
        for hemisphere in hemispheres:
            predicted_map_hemi = []
            empirical_map_hemi = []

            # Create mask based on the region of interest and hemisphere
            ROI_masked, mask = create_mask(final_mask_L_ROI, final_mask_R_ROI, final_mask_L, final_mask_R, hemisphere)
            if region_of_interest == 'all':
                mask_final = ROI_masked
            else:
                mask_final = mask
            predicted_map_hemi, empirical_map_hemi = process_subjects(list_of_sub_ids, path, dataset_name, hemisphere, mask_final, retinotopic_map, experiment, threshold)
            predicted_map = predicted_map + list(np.array(predicted_map_hemi).flatten())
            empirical_map = empirical_map + list(np.array(empirical_map_hemi).flatten())
        
        predicted_map, empirical_map = remove_outliers(predicted_map, empirical_map, retinotopic_map, dataset_name)
        
        if threshold != None:
            print('Threshold: ' + str(threshold))
        data = {'DeepRetinotopy': predicted_map, 'Empirical': empirical_map}
        data = pd.DataFrame(data)
        fig = plt.figure(figsize=(10, 10))
        if plot_type == 'kde':
            sns.kdeplot(x = data['DeepRetinotopy'], y = data['Empirical'],cmap="Blues", fill=True,cbar=True)
        else:
            sns.jointplot(x =data['DeepRetinotopy'], y = data['Empirical'], color='blue', kind='hex')
        
        if retinotopic_map == 'polarAngle':
            tmp_empirical_map = np.array(data['Empirical'])
            tmp_predicted_map = np.array(data['DeepRetinotopy'])
            tmp_predicted_map[np.isnan(tmp_empirical_map)] = 0
            tmp_empirical_map[np.isnan(tmp_empirical_map)] = 0
            corr = circcorrcoef(tmp_predicted_map*u.deg, tmp_empirical_map*u.deg)
            print(corr)
        else:
            corr = scipy.stats.pearsonr(data['DeepRetinotopy'], data['Empirical'])
            print(corr)
        plt.xlabel('DeepRetinotopy ($^\circ$)', fontsize=17)
        plt.ylabel('Empirical ($^\circ$)', fontsize=17)

        if retinotopic_map ==  'polarAngle':
            plt.plot([0, 360], [0, 360], 'k--')
            plt.ylim(0,360)
            plt.xlim(0,360)
            if len(hemispheres) > 1:
                fig.suptitle('Polar Angle')
            else:
                fig.suptitle('Polar Angle - ' + hemispheres[0])
        elif retinotopic_map == 'eccentricity':
            plt.plot([0, 10], [0, 10], 'k--')
            plt.ylim(0,10)
            plt.xlim(0,10)
            if len(hemispheres) > 1:
                fig.suptitle('Eccentricity')
            else:
                fig.suptitle('Eccentricity - ' + hemispheres[0])
        elif retinotopic_map == 'pRFsize':
            plt.plot([0, 3], [0, 3], 'k--')
            plt.ylim(0,3)
            plt.xlim(0,3)
            if len(hemispheres) > 1:
                fig.suptitle('pRF size')
            else:
                fig.suptitle('pRF size - ' + hemispheres[0])

        # Create file name based on parameters
        ## Experiment name to file
        if dataset_name == 'logbar':
            experiment_name = '_' + experiment 
        else:
            experiment_name = ''
        ## Pool of participants to file
        if pool_of_participants == 'all':
            pool_of_participants_name = ''
        elif pool_of_participants == 'children':
            pool_of_participants_name = '_children'
        elif pool_of_participants == 'adults':
            pool_of_participants_name = '_adults'
        ## Base path
        base_path = f'../output/model_evaluation/{dataset_name}/PredictedVsEmpirical_{retinotopic_map}{experiment_name}'
        ## Create file suffix
        if len(hemispheres) > 1:
            hemisphere_name = 'both'
        else:
            hemisphere_name = hemispheres[0]
        threshold_str = f'_{threshold}' if threshold is not None else ''
        file_suffix = f'{hemisphere_name}{threshold_str}_{region_of_interest}{pool_of_participants_name}_{plot_type}'
        plt.savefig(f'{base_path}_{file_suffix}.png')
        plt.savefig(f'{base_path}_{file_suffix}.pdf')
        print(f'Saving plot to {base_path}_{file_suffix}.png and {base_path}_{file_suffix}.pdf')
        plt.show()
        # plt.show()
    return

def explainedvariance_vs_error(path, retinotopic_map, hemisphere, threshold = 10):
    """Calculate the inter-individual variability in predicted maps and the error
    between predicted and empirical maps for different models.
    
    Args:
        path (str): Path to the Freesurfer directory (both empirical and predicted maps should be 
            located at the deepRetinotopy folder within the individual's directory)
        retinotopic_map (str): Type of retinotopic map ('polarAngle', 'eccentricity' or 'pRFsize')
        hemisphere (str): Hemisphere ('lh' or 'rh')
        threshold (float): Threshold for the explained variance

    Returns:
        plt (matplotlib.pyplot): Plot of the error vs explained variance.
    """

    # Development dataset
    dev_set = ['186949', '169747', '826353', '825048', '671855',
                    '751550', '318637', '131722', '137128', '706040'] 
    
    ## Region of interest 
    # Region of interest used for training
    final_mask_L_ROI, final_mask_R_ROI, index_L_mask, index_R_mask = roi(['ROI'])
    ROI_masked = np.zeros((32492, 1))
    # Early visual cortex
    final_mask_L, final_mask_R, index_L_mask, index_R_mask = roi_earlyvisualcortex(['ROI'])
    earlyVisualCortex = np.zeros((32492, 1))
    # Hemisphere
    if hemisphere == 'lh':
        ROI_masked[final_mask_L_ROI == 1] = 1
        earlyVisualCortex[final_mask_L == 1] = 1
    else: 
        ROI_masked[final_mask_R_ROI == 1] = 1
        earlyVisualCortex[final_mask_R == 1] = 1
    # Final mask (selecting V1-V3 vertices)
    mask = ROI_masked + earlyVisualCortex
    mask = mask[ROI_masked == 1]
    mask = mask > 1

    # Number of nodes
    number_cortical_nodes = int(64984)
    number_hemi_nodes = int(number_cortical_nodes / 2)

    for model in range(5):
        errors_per_sub = []
        variance_explained_sub = []
        for j in range(len(dev_set)):
            for i in range(len(dev_set)):
                # Error
                if i == j:
                    # Loading maps
                    data = RetinotopyData(path, dev_set[i], hemisphere, retinotopic_map, model_index=model + 1)
                    data.apply_mask_to_maps(ROI_masked)
                    data.apply_mask_to_maps(mask)
                    
                    # Transform polar angle values from lh
                    if retinotopic_map == 'polarAngle' and hemisphere == 'lh':
                        data.apply_transform_polarangle()

                    assert np.min(data.predicted_map) >= 0.
                    assert np.min(data.empirical_map) >= 0.
                    
                    # Transforming to radians
                    if retinotopic_map == 'polarAngle' or retinotopic_map == 'eccentricity':
                        data.convert_to_radian()

                    # Calculating error
                    if retinotopic_map == 'polarAngle' or retinotopic_map == 'eccentricity':
                        theta = smallest_angle(data.predicted_map, data.empirical_map)
                        if threshold != None:
                            theta = theta[data.variance_explained > threshold]
                            data.variance_explained = data.variance_explained[data.variance_explained > threshold]
                    elif retinotopic_map == 'pRFsize':
                        theta = np.abs(data.predicted_map - data.empirical_map)
                        if threshold != None:
                            theta = theta[data.variance_explained > threshold]
                            data.variance_explained = data.variance_explained[data.variance_explained > threshold]
                    # data for histogram plot
                    errors_per_sub = errors_per_sub + list(theta)

                    variance_explained_sub = variance_explained_sub + list(data.variance_explained)
        
        errors_per_sub = np.array(errors_per_sub).flatten()
        variance_explained_sub = np.array(variance_explained_sub).flatten()

        data = [errors_per_sub[(variance_explained_sub < threshold + 10) & (variance_explained_sub > threshold)] for threshold in range(0, 100, 10)]
        data = [[np.mean(data[i]), i * 10] for i in range(len(data))]
        df = pd.DataFrame(columns=['Mean error', 'Variance explained'],
                      data=data)
        
        if os.path.isdir('../output/error_vs_explained_variance/') == False:
            os.makedirs('../output/error_vs_explained_variance/')

        sns.set_style("white")
        fig = plt.figure()
        ax = sns.barplot(y='Mean error', x='Variance explained', data=df, palette="Reds_r")
        ax.set_xticklabels(range(10, 110, 10))
        fig.suptitle('Early visual areas in the ' + hemisphere + ' - Seed ' + str(model + 1))
        if retinotopic_map == 'eccentricity':
            plt.ylim([0, 2])
        elif retinotopic_map == 'pRFsize':
            plt.ylim([0, 1])
        else:
            plt.ylim([0, 80])
        plt.savefig('../output/error_vs_explained_variance/' + retinotopic_map + '_' + hemisphere + '_Seed' + str(model + 1) + '.pdf')
        plt.show()
    return


def variance_explained(stim_time, x, y, stim_img, fmri_data, tr, y_coord, x_coord, pRFsize):
    '''
    Calculate the variance explained using pRF model to predict the fMRI data.

    Args:
        stim_time (np.array): Time of the stimulus
        x (np.array): x coordinate of the stimulus
        y (np.array): y coordinate of the stimulus
        stim_img (np.array): Stimulus image
        fmri_data (np.array): fMRI data
        tr (float): Repetition time
        y_coord (np.array): y coordinate of the pRF
        x_coord (np.array): x coordinate of the pRF
        pRFsize (np.array): pRF size

    Returns:
        variance_explained (np.array): Variance explained by the pRF model
    '''

    # Generate HRF
    t_max = int(np.max(stim_time)) + 1
    hrf = _gamma_difference_hrf(tr=tr, time_length=t_max, onset=0.1, 
                                  delay = 6.68, undershoot=14.66, dispersion=1.82, u_dispersion=3.15,
                                  ratio=1/3.08)
    
    # Interpolate HRF to match stimulus time
    hrf_function = scipy.interpolate.interp1d(np.linspace(0, t_max, len(hrf)), hrf)
    new_hrf = hrf_function(stim_time)
    new_hrf = np.reshape(new_hrf, (1, -1))

    # Generate x and y grid
    pos = np.dstack((x, y))
    
    # Compute variance explained of the pRF model with provided pRF parameters
    variance_explained_ = []
    for i in range(len(x_coord)):
        X0 = x_coord[i]
        Y0 = y_coord[i]
        sigma = pRFsize[i]
        mean = [X0, Y0]
        if sigma !=0:
            # Generate covariance matrix
            covariance = [[sigma, 0], [0, sigma]] 

            # Generate 2D Gaussian
            rv = multivariate_normal(mean, covariance)
            gaussian_2d = rv.pdf(pos)

            # Determine the overlap between the stimulus and the pRF and calculate the spatial summation
            overlap = stim_img * gaussian_2d[:,:,np.newaxis]
            overlap = np.sum(overlap, axis=(0,1))

            # Convolve the timeseries from the previous step with the HRF
            predicted_signal = np.convolve(overlap, new_hrf[0])[0:len(overlap)]

            # Interpolate the predicted signal to match the fmri data
            bold_func_pred = scipy.interpolate.interp1d(np.linspace(0, t_max, len(predicted_signal)), predicted_signal)
            predicted_signal = bold_func_pred(np.linspace(0, t_max, int(t_max/tr)))

            # Calculate variance explained
            if np.std(predicted_signal) == 0:
                corr = 0
            else:
                norm_predicted_signal = (predicted_signal - np.mean(predicted_signal)) / np.std(predicted_signal)
                norm_fmri_data = (fmri_data[i] - np.mean(fmri_data[i])) / np.std(fmri_data[i])
                corr = np.corrcoef(norm_predicted_signal, norm_fmri_data)[0,1]
        else:
            corr = 0    
        variance_explained_.append(corr**2)
    return np.array(variance_explained_)