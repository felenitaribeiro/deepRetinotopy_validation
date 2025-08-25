import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import argparse



def transform_angle_native_neuropythy(path_to_use, path_to_save, hemisphere):
    img = nib.load(path_to_use)
    template =img.agg_data()

    #switch 180-360 degrees to -180-0 degrees
    over_180 = template > 180
    template[over_180] = template[over_180] - 360

    # template * -1
    template = template * -1

    # change from degrees to radians
    template = template * np.pi / 180
    # rotate about -90 degrees
    x = np.cos(template)
    y = np.sin(template)
    rotation_matrix = np.array([[0, -1], [1, 0]])
    rotated_coords = rotation_matrix @ np.array([x, y])
    rotated_angle = np.degrees(np.arctan2(rotated_coords[1], rotated_coords[0]))

    img.agg_data()[:] = rotated_angle
    nib.save(img, path_to_save)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("method", choices=["transform_polarangle_neuropythy", "retransform_angle_360_180_180", "transform_angle_native_neuropythy"])
    parser.add_argument("--path_to_use", type=str)
    parser.add_argument("--path_to_save", type=str)
    parser.add_argument("--hemisphere", type=str)
    args = parser.parse_args()

    
    if args.method == "transform_angle_native_neuropythy":
        transform_angle_native_neuropythy(args.path_to_use, args.path_to_save, args.hemisphere)


