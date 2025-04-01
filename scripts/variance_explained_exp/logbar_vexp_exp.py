import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.interpolate
import nibabel as nib
import sys
import os
sys.path.append('../..')
from functions.evaluation import variance_explained

def data_load(data_directory, sub_id, num_of_runs, num_of_sessions, experiment = 'logbar', roi='all'):
    '''
    Load data for the logbar experiment
    
    Args:
        data_directory (str): path to the data directory
        sub_id (str): subject id
        num_of_runs (int): number of runs
        num_of_sessions (int): number of sessions
        experiment (str): experiment type (logbar or fixedbar)
        roi (str): region of interest (all or fovea)
        
    Returns:
        x_empirical (np.array): x coordinates of empirical data
        y_empirical (np.array): y coordinates of empirical data
        x_pred (np.array): x coordinates of predicted data
        y_pred (np.array): y coordinates of predicted data
        sigma_empirical (np.array): empirical pRF size
        sigma_pred (np.array): predicted pRF size
        fmri_data (np.array): fMRI data
        vexpl_empirical (np.array): variance explained of empirical data
        stim_img (np.array): stimulus images
        stim_time (np.array): stimulus time
        x (np.array): x coordinates of the stimulus
        y (np.array): y coordinates of the stimulus
        mask_final (np.array): mask for the data
    '''

    # Load stimulus
    stim_time = []
    stim_img = []
    t_final = 0
    for run in ['run-1', 'run-2', 'run-3']:
        if experiment == 'fixedbar':
            stimulus = scipy.io.loadmat(data_directory + 'derivatives/prf-estimation/stimuli/'
                                'task-logbar_' + run + '_desc-down_aperture.mat')
            print('Loading logbar stimulus')
        elif experiment == 'logbar':
            stimulus = scipy.io.loadmat(data_directory + 'derivatives/prf-estimation/stimuli/'
                                'task-fixedbar_' + run + '_desc-down_aperture.mat')
            print('Loading fixedbar stimulus')
        stim_time_tmp = stimulus['stim']['t'][0][0]

        stim_time.append(stim_time_tmp + t_final)
        stim_img.append(stimulus['stim']['img'][0][0])
        t_final += stim_time_tmp[0][-1] + 1 # stimulus starts at 0
        if run == 'run-3':
            x = stimulus['stim']['x'][0][0]
            y = stimulus['stim']['y'][0][0]
    
    if num_of_sessions > 1:
        stim_img = np.concatenate(stim_img, axis=2)
        stim_img = [stim_img] * num_of_sessions
        stim_img = np.concatenate(stim_img, axis=2)
        
        stim_time = list(np.concatenate(stim_time, axis=1))
        for i in range(1, num_of_sessions):
            stim_time.append(stim_time[0] + i + stim_time[0][-1] * i)
        stim_time = np.concatenate(stim_time, axis=0)
    else:
        stim_img = np.concatenate(stim_img, axis=2)
        stim_time = np.concatenate(stim_time, axis=1)
    assert np.unique(stim_time, return_counts=True)[1].max() == 1

    # Load fMRI data
    fmri_data = []
    runs_list = ['run-' + str(i) for i in range(1, num_of_runs + 1)]
    ses_list = ['ses-0' + str(i) for i in range(1, num_of_sessions + 1)]
    for ses in ses_list:
        for run in runs_list:
            if experiment == 'fixedbar':
                fmri_data_tmp = nib.load(data_directory + 'derivatives/prf-estimation/' + sub_id +
                                    '/func/' + sub_id + '_' + ses + '_task-logbar_' + run + 
                                    '_hemi-L_space-fsnative_desc-preproc_bold.func.gii').agg_data()
                print('Loading logbar fMRI data')
            elif experiment == 'logbar':
                fmri_data_tmp = nib.load(data_directory + 'derivatives/prf-estimation/' + sub_id +
                                        '/func/' + sub_id + '_' + ses + '_task-fixedbar_' + run + 
                                        '_hemi-L_space-fsnative_desc-preproc_bold.func.gii').agg_data()
                print('Loading fixedbar fMRI data')
            fmri_data_tmp = np.array(fmri_data_tmp)
            fmri_data.append(fmri_data_tmp)
    fmri_data = np.concatenate(fmri_data, axis=1)

    # Load empirical data    
    polar_angle_empirical = nib.load(data_directory + 'derivatives/prf-estimation/' + sub_id + '/' + '/braincoder/' +
                        sub_id + '_ses-concat_task-' + experiment + '_run-concat_hemi-l_iter_pol.gii').agg_data()
    eccentricity_empirical = nib.load(data_directory + 'derivatives/prf-estimation/' + sub_id + '/braincoder/' + 
                        sub_id + '_ses-concat_task-' + experiment + '_run-concat_hemi-l_iter_ecc.gii').agg_data()
    vexpl_empirical = nib.load(data_directory + 'derivatives/prf-estimation/' + sub_id + '/braincoder/' + 
                        sub_id + '_ses-concat_task-' + experiment + '_run-concat_hemi-l_iter_r2.gii').agg_data()
    pRFsize_empirical = nib.load(data_directory + 'derivatives/prf-estimation/' + sub_id + '/braincoder/' + 
                        sub_id + '_ses-concat_task-' + experiment + '_run-concat_hemi-l_iter_sd.gii').agg_data()
    
    
    # Load predicted data
    polar_angle = nib.load(data_directory + 'derivatives/freesurfer/' + sub_id + 
                           '/deepRetinotopy/' + sub_id + '.predicted_polarAngle_model.lh.native.func.gii').agg_data()
    eccentricity = nib.load(data_directory + 'derivatives/freesurfer/' + sub_id + 
                            '/deepRetinotopy/' + sub_id + '.predicted_eccentricity_model.lh.native.func.gii').agg_data()
    pRFsize = nib.load(data_directory + 'derivatives/freesurfer/' + sub_id + 
                       '/deepRetinotopy/' + sub_id + '.predicted_pRFsize_model.lh.native.func.gii').agg_data()


    # Create mask 
    mask = pRFsize > 0 # mask based on ROI from predicted maps
    mask = mask * 1
    mask_final = mask

    # Apply mask
    ## polar angle
    polar_angle_empirical = polar_angle_empirical[mask_final == 1]
    polar_angle = polar_angle[mask_final == 1]
    ## eccentricity
    eccentricity_empirical = eccentricity_empirical[mask_final == 1]
    eccentricity = eccentricity[mask_final == 1]
    ## pRF size
    sigma_pred = pRFsize[mask_final == 1]
    sigma_empirical = pRFsize_empirical[mask_final == 1]
    ## fMRI data
    fmri_data = fmri_data[mask_final == 1]
    vexpl_empirical = vexpl_empirical[mask_final == 1] 

    # Transforms
    ## Transform polar angle to radians
    polar_angle_empirical = polar_angle_empirical * np.pi/180
    polar_angle = polar_angle * np.pi/180
    ## Transform data to cartesian coordinates
    x_empirical = eccentricity_empirical * np.cos(polar_angle_empirical)
    y_empirical = eccentricity_empirical * np.sin(polar_angle_empirical)

    x_pred = eccentricity * np.cos(polar_angle)
    y_pred = eccentricity * np.sin(polar_angle)
    
    assert np.shape(x_pred) == np.shape(x_empirical)
    assert np.shape(y_pred) == np.shape(y_empirical)
    assert np.shape(sigma_pred) == np.shape(sigma_empirical)
    print('Data loaded')
    return x_empirical, y_empirical, x_pred, y_pred, sigma_empirical, sigma_pred, fmri_data, vexpl_empirical, stim_img, stim_time, x, y, mask_final

def main(args):
    x_empirical, y_empirical, x_pred, y_pred, sigma_empirical, sigma_pred, fmri_data, vexpl_empirical, stim_img, stim_time, x, y, mask = data_load(
        args.data_directory, args.sub_id, args.num_of_runs, args.num_of_sessions, args.experiment)
    if os.path.isdir('../../output/model_evaluation/logbar/explained_var/') == False:
        os.mkdir('../../output/model_evaluation/logbar/explained_var/')
    
    variance_explained_pred = variance_explained(stim_time, x, y, stim_img, fmri_data, args.tr, y_pred, x_pred, sigma_pred)
    np.save('../../output/model_evaluation/logbar/explained_var/' + args.sub_id + '_expvar_deepRetinotopy_' + args.experiment + '.npy', variance_explained_pred)
    
    variance_explained_empirical = variance_explained(stim_time, x, y, stim_img, fmri_data, args.tr, y_empirical, x_empirical, sigma_empirical)
    np.save('../../output/model_evaluation/logbar/explained_var/' + args.sub_id + '_expvar_empirical_' + args.experiment + '.npy', variance_explained_empirical)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Variance explained experiment')
    parser.add_argument('--data_directory', type=str, help='Path to the data directory', 
                        default='/BULK/LABDATA/openneuro/ds004698/')
    parser.add_argument('--sub_id', type=str, help='Subject ID', default='sub-01')
    parser.add_argument('--num_of_runs', type=int, help='Number of runs', default=3)
    parser.add_argument('--num_of_sessions', type=int, help='Number of sessions', default=2)
    parser.add_argument('--experiment', type=str, help='Experiment type (logbar or fixedbar)')
    parser.add_argument('--tr', type=float, help='Repetition time', default=1.2)
    args = parser.parse_args()
    main(args)