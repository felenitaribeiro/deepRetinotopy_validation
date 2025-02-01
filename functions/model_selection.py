import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
import os.path as osp
import seaborn as sns
import pandas as pd
import sys
import scipy

sys.path.append('..')
from deepRetinotopy_TheToolbox.utils.metrics import smallest_angle
from deepRetinotopy_TheToolbox.utils.rois import ROI_WangParcelsPlusFovea as roi
from sklearn.metrics import jaccard_score


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

def calculate_overlap(map_1, map_2, retinotopic_map, mask, angle = 'rad'):
    if angle == 'rad':
        map_1 = map_1 * (180/np.pi)
        map_2 = map_2 * (180/np.pi)
    
    if retinotopic_map == 'polarAngle':
        map_1[(map_1 >= 0) & (map_1 <= 45)] = 0 
        map_1[(map_1 > 45) & (map_1 <= 180)] = 90 
        map_1[(map_1 >= 315) & (map_1 <= 360)] = 360 
        map_1[(map_1 > 180) & (map_1 < 315)] = 270

        map_2[(map_2 >= 0) & (map_2 <= 45)] = 0 
        map_2[(map_2 > 45) & (map_2 <= 180)] = 90 
        map_2[(map_2 >= 315) & (map_2 <= 360)] = 360 
        map_2[(map_2 > 180) & (map_2 < 315)] = 270
    elif retinotopic_map == 'eccentricity':
        map_1[(map_1 >= 0) & (map_1 <= 2)] = 0
        map_1[(map_1 > 2) & (map_1 <= 4)] = 2
        map_1[(map_1 > 4) & (map_1 <= 6)] = 4
        map_1[(map_1 > 6)] = 6

        map_2[(map_2 >= 0) & (map_2 <= 2)] = 0
        map_2[(map_2 > 2) & (map_2 <= 4)] = 2
        map_2[(map_2 > 4) & (map_2 <= 6)] = 4
        map_2[(map_2 > 6)] = 6

    # elif retinotopic_map == 'pRFsize':
        # TODO

    jaccard_1 = jaccard_score(map_1[mask],
                            map_2[mask],
                            average='weighted')
    jaccard_2 = jaccard_score(map_2[mask],
                            map_1[mask],
                            average='weighted')
    score = (jaccard_1 + jaccard_2) / 2
    return score

def transform_polarangle(data):
    subtract = data > 180
    add = data < 180
    data[subtract] = data[subtract] - 180
    data[add] = data[add] + 180
    return data

def metric_model_selection(path, retinotopic_map, hemisphere, retinotopic_mapping = 'continuous'):
    
    mean_ind_variability = []
    error = []

    # Development dataset
    dev_set = ['186949', '169747', '826353', '825048', '671855',
                    '751550', '318637', '131722', '137128', '706040'] 
    
    # Number of nodes
    number_cortical_nodes = int(64984)
    number_hemi_nodes = int(number_cortical_nodes / 2)

    for model in range(5):
        theta_withinsubj = []
        theta_acrosssubj_pred = []

        # Region of interest used for training
        ROI = ['ROI']
        final_mask_L, final_mask_R, index_L_mask, index_R_mask = roi(
            ROI)
        ROI_masked = np.zeros((32492, 1))
        # early visual cortex
        final_mask_L, final_mask_R, index_L_mask, index_R_mask = roi_earlyvisualcortex(['ROI'])
        dorsal_earlyVisualCortex = np.zeros((32492, 1))

        if hemisphere == 'lh':
            ROI_masked[final_mask_L == 1] = 1
            dorsal_earlyVisualCortex[final_mask_L == 1] = 1
        else: 
            ROI_masked[final_mask_R == 1] = 1
            dorsal_earlyVisualCortex[final_mask_R == 1] = 1

        # Final mask (selecting dorsal V1-V3 vertices)
        mask = ROI_masked + dorsal_earlyVisualCortex
        mask = mask[ROI_masked == 1]
        mask = mask > 1

        for j in range(len(dev_set)):
            theta_pred_across_temp = []
            for i in range(len(dev_set)):
                # Error
                if i == j:
                    # Loading predicted values
                    predicted_map = np.array(nib.load(osp.join(path,
                                     dev_set[i] + '/deepRetinotopy/' + dev_set[i] + '.fs_predicted_' + 
                                     retinotopic_map + '_' + hemisphere + '_curvatureFeat_model' + str(model + 1) + '.func.gii')).agg_data()).reshape(
                                     number_hemi_nodes, -1)[ROI_masked == 1]
                    empirical_map = np.array(nib.load(osp.join(path,
                                     dev_set[i] + '/deepRetinotopy/' + dev_set[i] + '.fs_empirical_' + 
                                     retinotopic_map + '_' + hemisphere + '_masked.func.gii')).agg_data()).reshape(
                                     number_hemi_nodes, -1)[ROI_masked == 1]
                    if retinotopic_map == 'polarAngle' and hemisphere == 'lh':
                        predicted_map = transform_polarangle(predicted_map)
                        predicted_map = np.array(predicted_map) * (np.pi / 180)

                        empirical_map = transform_polarangle(empirical_map)
                        empirical_map = np.array(empirical_map) * (np.pi / 180)

                    if retinotopic_mapping == 'continuous':
                        theta = smallest_angle(predicted_map[mask],
                                            empirical_map[mask])
                        theta_withinsubj.append(np.mean(theta))
                    elif retinotopic_mapping == 'discrete':
                        if retinotopic_map == 'polarAngle':
                            overlap = calculate_overlap(predicted_map, empirical_map, retinotopic_map, mask, angle = 'rad')
                            theta_withinsubj.append(overlap)
                # Inter-individual variability in predicted maps
                if i != j:
                    # Loading predicted values
                    predicted_map_1 = np.array(nib.load(osp.join(path,
                                        dev_set[i] + '/deepRetinotopy/' + dev_set[i] + '.fs_predicted_' + 
                                        retinotopic_map + '_' + hemisphere + '_curvatureFeat_model' + str(model + 1) + '.func.gii')).agg_data()).reshape(
                                        number_hemi_nodes, -1)[ROI_masked == 1]
                    predicted_map_2 = np.array(nib.load(osp.join(path,
                                        dev_set[j] + '/deepRetinotopy/' + dev_set[j] + '.fs_predicted_' + 
                                        retinotopic_map + '_' + hemisphere + '_curvatureFeat_model' + str(model + 1) + '.func.gii')).agg_data()).reshape(
                                        number_hemi_nodes, -1)[ROI_masked == 1]

                    if retinotopic_map == 'polarAngle' and hemisphere == 'lh':
                        predicted_map_1 = transform_polarangle(predicted_map_1)
                        predicted_map_1 = np.array(predicted_map_1) * (np.pi / 180)

                        predicted_map_2 = transform_polarangle(predicted_map_2)
                        predicted_map_2 = np.array(predicted_map_2) * (np.pi / 180)

                    if retinotopic_mapping == 'continuous':
                        theta_pred = smallest_angle(predicted_map_1[mask],
                                                    predicted_map_2[mask])
                        theta_pred_across_temp.append(
                            np.mean(theta_pred))
                    elif retinotopic_mapping == 'discrete':
                        if retinotopic_map == 'polarAngle':
                            overlap = calculate_overlap(predicted_map_1, predicted_map_2, retinotopic_map, mask, angle = 'rad')
                            theta_pred_across_temp.append(overlap)
            theta_acrosssubj_pred.append(theta_pred_across_temp)

        mean_theta_withinsubj = np.array(theta_withinsubj)
        # if retinotopic_mapping == 'continuous':
        #     mean_theta_acrosssubj_pred = np.mean(np.array(theta_acrosssubj_pred), axis=1)
        # else: 
        mean_theta_acrosssubj_pred = np.mean(np.array(theta_acrosssubj_pred), axis=1)

        error.append(mean_theta_withinsubj)
        mean_ind_variability.append(mean_theta_acrosssubj_pred)

    error = np.reshape(np.array(error), (5, -1))
    mean_ind_variability = np.reshape(np.array(mean_ind_variability), (5, -1))

    data_ind_variability = [[mean_ind_variability[i], len(mean_ind_variability[i]) * ['Model ' + str(i+1)],
                            len(mean_ind_variability[i]) * ['Individual variability']] for i in range(5)]
    data_error = [[error[i], len(mean_ind_variability[i]) * ['Model ' + str(i+1)],
                            len(mean_ind_variability[i]) * [
                                'Error']] for i in range(5)]
    data = np.concatenate(data_ind_variability + data_error, axis = 1)

    df = pd.DataFrame(columns=['$\Delta$$\t\Theta$', 'Model', 'Metric'],
                      data=data.T)
    df['$\Delta$$\t\Theta$'] = df['$\Delta$$\t\Theta$'].astype(float)

    palette = ['dimgray', 'lightgray']
    sns.set_style("whitegrid")
    fig = plt.figure()
    ax = sns.pointplot(y='$\Delta$$\t\Theta$', x='Model',
                       hue='Metric', data=df, palette=palette,
                       join=False, dodge=True, ci=95)
    ax.set_title('Early visual areas')
    if str(retinotopic_map) == 'eccentricity' or str(retinotopic_map) == 'pRFsize':
        plt.ylim([0, 1])
    else:
        if retinotopic_mapping == 'continuous':
            plt.ylim([0, 40])
        else:
            plt.ylim([0, 1])
    # plt.savefig('ModelSelection_DorsalV123.svg')
    plt.show()
    return df


