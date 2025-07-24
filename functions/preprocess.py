import numpy as np
import nibabel as nib
import sys
sys.path.append('..')

def polarcoord(x_data, y_data):
    path_to_save = x_data[:-6]
    # Load data
    template = nib.load(x_data)
    xs = template.agg_data()
    ys = nib.load(y_data).agg_data()

    # Transform to polar coordinates
    theta = np.arctan2(ys,xs) * 180 / np.pi 
    sum = theta < 0 # shift values to be between 0 and 360
    theta[sum] = theta[sum] + 360

    r = np.sqrt(xs**2 + ys**2)
    print(theta.max(), theta.min())
    # Save data
    template.agg_data()[:] = theta
    nib.save(template, path_to_save + 'angle_new.gii')
    template.agg_data()[:] = r
    nib.save(template, path_to_save + 'eccen_new.gii')
    return print('Data in cartesian coordinates transformed to polar coordinates')


def transform_angle(path_to_empirical_data, hemisphere, radians = False, left_hemi_shift = True):
    """
    Transform the polar angle maps from -180 to 180 degrees where the origin in the positive x-axis, to 0 to 360 degrees where
    the origin is the positive x-axis. The angles of the left hemisphere will be shifted by 180 degrees.
    """
    path_to_empirical_data = str(path_to_empirical_data)
    if left_hemi_shift:
        path_to_save = path_to_empirical_data[:-4] + '_transformed.gii'
    else:
        path_to_save = path_to_empirical_data[:-4] + '_0-360_transformed.gii'
    
    # Load the empirical data
    template = nib.load(path_to_empirical_data)
    data = template.agg_data()
    if radians:
        data = data * 180 / np.pi

    # shift values to be between 0 and 360
    sum = data < 0
    data[sum] = data[sum] + 360

    if hemisphere == 'lh' and left_hemi_shift:
        # Rescaling polar angle values
        sum_180 = data < 180
        minus_180 = data > 180
        data[sum_180] = data[sum_180] + 180
        data[minus_180] = data[minus_180] - 180
    template.agg_data()[:] = data

    nib.save(template, path_to_save)
    return 'Transformed data saved as ' + path_to_save

def transform_polarangle_benson14(path, hemisphere = 'lh'): 
    """
    Transform the polar angle maps from Neuropythy convention (LH: 0-180 referring to UVM -> RHM -> LVM; 
      RH: 0-180 referring to UVM -> LHM -> LVM) to standard angle representation from 0 to 360 degrees where
      the origin is the positive x-axis.
    
    Parameters
    ----------
    path : str
        The path to polar angle map file.
    hemisphere : str, optional
        The hemisphere of the polar angle map, either 'lh' for left hemisphere or 'rh' for right hemisphere.
        Default is 'lh'.
        
    Returns
    -------
    numpy.ndarray
        The transformed polar angle in degrees.
    """

    data = nib.load(path)
    angle = data.agg_data()
    if hemisphere == 'lh':
        angle = - angle
    mask = angle == 0
    # Step 1: Rotate by 90 degrees using a rotation matrix
    angle = angle / 180 * np.pi
    x = np.cos(angle)
    y = np.sin(angle)
    rotation_matrix = np.array([[0, -1], [1, 0]])
    rotated_coords = rotation_matrix @ np.array([x, y])
    rotated_angle = np.degrees(np.arctan2(rotated_coords[1], rotated_coords[0]))

    # Step 2: Shift values to be between 0 and 360
    rotated_angle[rotated_angle <= 0] = np.abs(rotated_angle[rotated_angle <= 0] + 360)
    rotated_angle[mask] = 0
    data.agg_data()[:] = rotated_angle
    file_name = path[:-4] + '_neuropythy.gii'

    nib.save(data, file_name)

    return print('Polar angle map has been transformed and saved as ' + file_name)

def transform_polarangle_to_benson14(path, hemisphere = 'lh'):
    """
    Transform the polar angle maps from standard angle representation from 0 to 360 degrees where
      the origin is the positive x-axis to Neuropythy convention (LH: 0-180 referring to UVM -> RHM -> LVM; 
      RH: 0-180 referring to UVM -> LHM -> LVM).
    
    Parameters
    ----------
    path : str
        The path to polar angle map file.
    hemisphere : str, optional
        The hemisphere of the polar angle map, either 'lh' for left hemisphere or 'rh' for right hemisphere.
        Default is 'lh'.
        
    Returns
    -------
    numpy.ndarray
        The transformed polar angle in degrees.
    """
    data = nib.load(path)
    angle = data.agg_data()

    #switch 180-360 degrees to -180-0 degrees
    over_180 = angle > 180
    angle[over_180] = angle[over_180] - 360

    if hemisphere == 'lh':
        angle = - angle
    mask = angle == 0
    # Step 1: Rotate by 90 degrees using a rotation matrix
    angle = angle * np.pi /180
    x = np.cos(angle)
    y = np.sin(angle)
    rotation_matrix = np.array([[0, -1], [1, 0]])
    rotated_coords = rotation_matrix @ np.array([x, y])
    rotated_angle = np.degrees(np.arctan2(rotated_coords[1], rotated_coords[0]))

    # Step 2: Shift values to be between 0 and 360
    rotated_angle[mask] = 0
    data.agg_data()[:] = rotated_angle
    file_name = path[:-22] + '_180-180_neuropythy.gii'

    nib.save(data, file_name)

    return print('Polar angle map has been transformed and saved as ' + file_name)