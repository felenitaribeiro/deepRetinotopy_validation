import numpy as np
import nibabel as nib
import sys
import argparse
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


def transform_angle(path_to_empirical_data, hemisphere, radians = False, left_hemi_shift = False):
    """
    Transform the polar angle maps from -180 to 180 degrees, to 0 to 360 degrees, where the
    origin is the positive x-axis in both cases. This leaves the maps in the natural
    convention, in which the left hemisphere represents the right visual field, matching what
    the toolbox predicts.

    left_hemi_shift additionally shifts the left hemisphere by 180 degrees, which is the
    convention the legacy models were trained on. It is only there to reproduce the empirical
    maps as they used to be generated, and should be left off.
    """
    path_to_empirical_data = str(path_to_empirical_data)
    path_to_save = path_to_empirical_data[:-4] + '_transformed.gii'

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

def transform_polarangle_benson14(path, hemisphere = 'lh', output_path = None): 
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
    output_path : str, optional
        Where to save the transformed map. Default is the input path with a '_neuropythy.gii' suffix.

    Notes
    -----
    Apply this AFTER resampling to the 32k surface, not before. In the natural convention the
    left hemisphere wraps at 0/360 on the horizontal meridian, inside the represented hemifield,
    and barycentric resampling across that wrap produces values from the wrong hemifield. The
    neuropythy convention (0-180) has no wrap inside the hemifield, so resample that instead
    and convert the resampled map (see scripts/regenerate_benson14_polarangle.sh).

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
    file_name = output_path if output_path is not None else path[:-4] + '_neuropythy.gii'

    nib.save(data, file_name)

    return print('Polar angle map has been transformed and saved as ' + file_name)

def transform_polarangle_to_benson14(path, path_to_save = None, deepretinotopy_data = False, hemisphere = 'lh'):
    """
    Transform the polar angle maps from standard angle representation from 0 to 360 degrees where
      the origin is the positive x-axis to Neuropythy convention (LH: 0-180 referring to UVM -> RHM -> LVM; 
      RH: 0-(-)180 referring to UVM -> LHM -> LVM).
    
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
    
    mask = angle == 0

    #switch 180-360 degrees to -180-0 degrees
    over_180 = angle > 180
    angle[over_180] = angle[over_180] - 360

    # Step 1: Rotate by 90 degrees using a rotation matrix
    angle = angle * np.pi /180
    x = np.cos(angle)
    y = np.sin(angle)
    rotation_matrix = np.array([[0, -1], [1, 0]])
    rotated_coords = rotation_matrix @ np.array([x, y])
    rotated_angle = np.degrees(np.arctan2(rotated_coords[1], rotated_coords[0]))
    
    # Step 2: Reverse (mirror) coordinates
    if hemisphere == 'lh':
        rotated_angle = np.abs(rotated_angle - 180)
    if hemisphere == 'rh':
        rotated_angle = rotated_angle + 180
        
    # Step 3: Apply mask
    rotated_angle[mask] = 0
    data.agg_data()[:] = rotated_angle
    if path_to_save == None:
        file_name = path[:-22] + '_180-180_neuropythy.gii'
        if deepretinotopy_data:
            file_name = path[:-16] + '_180-180_neuropythy.gii'
    else:
        file_name = path_to_save
    nib.save(data, file_name)

    return print('Polar angle map has been transformed and saved as ' + file_name)


def polarangle_to_components(path, path_prefix):
    """Write the cosine and sine of a polar angle map as two metric files,
    {path_prefix}_cos.func.gii and {path_prefix}_sin.func.gii.

    Resample these instead of the angle itself: interpolating an angle across its 0/360
    wrap (the horizontal meridian of the left hemisphere in the natural convention)
    produces values from the wrong hemifield, whereas the components interpolate safely
    in any convention. Recombine with components_to_polarangle() after resampling.
    NaN vertices stay NaN, as they would when resampling the angle directly.
    """
    data = nib.load(path)
    angle = np.radians(data.agg_data().astype(float))
    for name, values in (('cos', np.cos(angle)), ('sin', np.sin(angle))):
        data.agg_data()[:] = values
        nib.save(data, f"{path_prefix}_{name}.func.gii")


def components_to_polarangle(cos_path, sin_path, path_to_save):
    """Recombine resampled cosine and sine maps into a polar angle map in degrees, 0-360,
    in the same convention as the map given to polarangle_to_components()."""
    data = nib.load(cos_path)
    cos_values = data.agg_data().astype(float)
    sin_values = nib.load(sin_path).agg_data().astype(float)
    angle = np.degrees(np.arctan2(sin_values, cos_values)) % 360
    data.agg_data()[:] = angle
    nib.save(data, path_to_save)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("method", choices=["transform_polarangle_to_benson14",
                                           "polarangle_to_components", "components_to_polarangle"])
    parser.add_argument("--path_to_use", type=str)
    parser.add_argument("--path_to_save", type=str)
    parser.add_argument("--hemisphere", type=str)
    parser.add_argument("--path_prefix", type=str, help="polarangle_to_components: output prefix")
    parser.add_argument("--cos_path", type=str, help="components_to_polarangle: resampled cosine map")
    parser.add_argument("--sin_path", type=str, help="components_to_polarangle: resampled sine map")
    args = parser.parse_args()

    
    if args.method == "transform_polarangle_to_benson14":
        transform_polarangle_to_benson14(args.path_to_use, args.path_to_save, hemisphere=args.hemisphere)
    elif args.method == "polarangle_to_components":
        polarangle_to_components(args.path_to_use, args.path_prefix)
    elif args.method == "components_to_polarangle":
        components_to_polarangle(args.cos_path, args.sin_path, args.path_to_save)