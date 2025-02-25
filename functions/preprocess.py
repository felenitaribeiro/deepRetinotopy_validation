import numpy as np
import nibabel as nib
import sys
sys.path.append('..')


def transform_polarangle_logbar(path):
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
    # Step 4: Apply the mask
    final_angle[mask] = -1
    data.agg_data()[:] = final_angle
    if path[-9:] == '.func.gii':
        file_name = path[:-9] + '_transformed.func.gii'
    else:
        file_name = path[:-4] + '_transformed.gii'
    nib.save(data, file_name)
    return print('Polar angle map has been transformed and saved as ' + file_name)

def transform_eccentricity_logbar(path):
    """
    Mask areas without estimates with -1.
    
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
    mask = angle > 10
    # Step 1: Apply the mask
    angle[mask] = -1
    data.agg_data()[:] = angle
    if path[-9:] == '.func.gii':
        file_name = path[:-9] + '_transformed.func.gii'
    else:
        file_name = path[:-4] + '_transformed.gii'
    nib.save(data, file_name)
    return print('Polar angle map has been transformed and saved as ' + file_name)