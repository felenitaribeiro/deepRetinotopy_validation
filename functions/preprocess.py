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
    theta = np.arctan2(ys,xs)
    r = np.sqrt(xs**2 + ys**2)
    print(theta.max(), theta.min())
    # Save data
    template.agg_data()[:] = theta
    nib.save(template, path_to_save + 'angle_new.gii')
    template.agg_data()[:] = r
    nib.save(template, path_to_save + 'eccen_new.gii')
    return print('Data in cartesian coordinates transformed to polar coordinates')


def transform_angle(path_to_empirical_data, hemisphere):
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

def transform_angle_lh_nsd(path_to_empirical_data):
    """
    Transform angles of the left hemisphere to avoid discontinuity.
    """
    path_to_empirical_data = str(path_to_empirical_data)
    path_to_save = path_to_empirical_data[:-4] + '_transformed.gii'
    
    # Load the empirical data
    template = nib.load(path_to_empirical_data)
    data = template.agg_data()

    # Rescaling polar angle values
    sum_180 = data < 180
    minus_180 = data > 180
    data[sum_180] = data[sum_180] + 180
    data[minus_180] = data[minus_180] - 180
    template.agg_data()[:] = data

    nib.save(template, path_to_save)
    return 'Transformed data saved as ' + path_to_save