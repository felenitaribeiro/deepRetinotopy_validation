import numpy as np
import nibabel as nib
import sys
sys.path.append('..')


def transform_polarangle_logbar(path, hemisphere): # TODO edit
    """
    Transform the polar angle maps from -180 to 180 degrees where the origin in the positive y-axis, to 0 to 360 degrees where
    the origin is the positive x-axis.
    
    Parameters
    ----------
    path : str
        The path to polar angle map file.
        
    Returns
    -------
    numpy.ndarray
        The transformed polar angle in degrees.
    """

    data = nib.load(path)
    angle = data.agg_data()
    # Step 0: Create a mask for angles greater than 180 degrees
    mask = angle > 180 
    # Step 1: Switch signal
    angle = -angle
    # Step 2: Add 180 degrees
    angle = (angle + 180) % 360
    # Step 3: Rotate by 90 degrees using a rotation matrix
    x = np.cos(np.radians(angle))
    y = np.sin(np.radians(angle))
    rotation_matrix = np.array([[0, 1], [-1, 0]])
    rotated_coords = rotation_matrix @ np.array([x, y])
    rotated_angle = np.degrees(np.arctan2(rotated_coords[1], rotated_coords[0]))
    final_angle = rotated_angle % 360
    
    # left hemisphere data in the same range as the predicted data
    if hemisphere == 'lh':
        # Rescaling polar angle values
        sum_180 = final_angle < 180
        minus_180 = final_angle > 180
        final_angle[sum_180] = final_angle[sum_180] + 180
        final_angle[minus_180] = final_angle[minus_180] - 180
    # Step 4: Apply the mask
    final_angle[mask] = -1
    data.agg_data()[:] = final_angle
    file_name = path[:-4] + '_transformed.gii'
    nib.save(data, file_name)
    return print('Polar angle map has been transformed and saved as ' + file_name)

def transform_eccensigma_logbar(path): # TODO edit
    """
    Mask areas without estimates with -1.
    
    Parameters
    ----------
    path : str
        The path to eccentricity map file.
        
    Returns
    -------
    numpy.ndarray
        The transformed eccentricity in degrees.
    """

    data = nib.load(path)
    angle = data.agg_data()
    # Step 0: Create a mask for angles greater than 180 degrees
    mask = angle > 10 
    # Step 1: Apply the mask
    angle[mask] = -1
    data.agg_data()[:] = angle
    file_name = path[:-4] + '_transformed.gii'
    nib.save(data, file_name)
    return print('Eccentricity map has been transformed and saved as ' + file_name)

def transform_angle_nyu(path_to_empirical_data, hemisphere):
    """
    Transform the polar angle maps from -180 to 180 degrees where the origin in the positive x-axis, to 0 to 360 degrees where
    the origin is the positive x-axis. The angles of the left hemisphere will be shifted by 180 degrees.
    """
    path_to_empirical_data = str(path_to_empirical_data)
    path_to_save = path_to_empirical_data[:-4] + '_transformed.gii'
    
    # Load the empirical data
    template = nib.load(path_to_empirical_data)
    data = template.agg_data()
    data = data * 180 / np.pi

    # shift values to be between 0 and 360
    sum = data < 0
    data[sum] = data[sum] + 360

    if hemisphere == 'lh':
        # Rescaling polar angle values
        sum_180 = data < 180
        minus_180 = data > 180
        data[sum_180] = data[sum_180] + 180
        data[minus_180] = data[minus_180] - 180

    template.agg_data()[:] = data

    nib.save(template, path_to_save)
    return 'Transformed data saved as ' + path_to_save