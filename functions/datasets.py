import os.path as osp
import os
import numpy as np
import nibabel as nib
import sys

sys.path.append(os.path.join(osp.dirname(osp.abspath(__file__)), '..'))

from functions.visualization import roi, roi_earlyvisualcortex

class RetinotopyData:
    def __init__(self, path, subject_id, hemisphere, 
                 retinotopic_map, number_hemi_nodes=int(32492), model = 'deepRetinotopy25', model_index=None, split_half=None):
        self.path = path
        self.subject_id = subject_id
        self.hemisphere = hemisphere
        self.retinotopic_map = retinotopic_map
        self.number_hemi_nodes = number_hemi_nodes
        self.model = model  # model type, e.g., 'deepRetinotopy' or 'benson14'
        self.model_index = model_index  # optional model index for specific seeds from deepRetinotopy
        self.split_half = split_half  # optional parameter for split-half analysis

        # Load maps during initialization
        self.predicted_map = self._load_map("predicted")
        self.empirical_map = self._load_map("empirical")
        self.variance_explained = self._load_map("variance_explained")
        self.curvature = self._load_map("curvature")

        if self.split_half is not None:
            self.empirical_map_split2 = self._load_map("empirical_split2")
            self.empirical_map_split3 = self._load_map("empirical_split3")

    def _load_map(self, map_type):
        """Load a retinotopic map (predicted or empirical)."""
        if map_type == 'predicted':
            if self.model == 'deepRetinotopy25':
                if self.model_index is not None:
                    file_name = f"deepRetinotopy/{self.subject_id}.fs_predicted_{self.retinotopic_map}_{self.hemisphere}_curvatureFeat_model{self.model_index}.func.gii"
                else:
                    file_name = f"deepRetinotopy/{self.subject_id}.fs_predicted_{self.retinotopic_map}_{self.hemisphere}_curvatureFeat_model.func.gii"
            elif self.model == 'deepRetinotopy21':
                file_name = f"predicted_deepRetinotopy_21/{self.subject_id}.fs_predicted_{self.retinotopic_map}_{self.hemisphere}_curvatureMyelinFeat_model.func.gii"
            elif self.model == 'benson14':
                file_name = f"surf/{self.hemisphere}.benson14_{self.retinotopic_map}.gii"
            elif self.model == 'noise_ceiling':
                file_name = f"deepRetinotopy/{self.subject_id}.fs_predicted_{self.retinotopic_map}_{self.hemisphere}_curvatureFeat_model.func.gii" # This won't be used
        elif map_type == 'empirical':
            file_name = f"surf/{self.subject_id}.fs_empirical_{self.retinotopic_map}_{self.hemisphere}.func.gii"
        elif map_type == 'variance_explained':
            file_name = f"surf/{self.subject_id}.fs_empirical_variance_explained_{self.hemisphere}.func.gii"
        elif map_type == 'curvature':
            file_name = f"surf/{self.subject_id}.curvature-midthickness.{self.hemisphere}.32k_fs_LR.func.gii"
        elif map_type == 'empirical_split2':
            file_name = f"surf/{self.subject_id}.fs_empirical_{self.retinotopic_map}_{self.hemisphere}_fit2.func.gii"
        elif map_type == 'empirical_split3':
            file_name = f"surf/{self.subject_id}.fs_empirical_{self.retinotopic_map}_{self.hemisphere}_fit3.func.gii"
        else:
            raise ValueError("Invalid map type specified.")

        if self.model == 'deepRetinotopy21' and map_type == 'predicted':
            file_path = osp.join(self.path, '..', file_name)
        else:
            file_path = osp.join(self.path, self.subject_id, file_name)
        if not osp.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
        
        data = np.array(nib.load(file_path).agg_data()).reshape(self.number_hemi_nodes, -1)
        return data
            
    def _apply_mask(self, data):
        """Apply the ROI mask to the data."""
        return data[self.mask == 1]
    
    def _transform_polarangle(self, data):
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
    
    def  _convert_to_radian(self, data):
        """Convert polar angle values from degrees to radians."""
        return data * np.pi / 180

    def apply_mask_to_maps(self, mask):
        """Apply the mask to the predicted, empirical, and variance explained maps."""
        self.mask = mask
        self.predicted_map = self._apply_mask(self.predicted_map)
        self.empirical_map = self._apply_mask(self.empirical_map)
        self.variance_explained = self._apply_mask(self.variance_explained)
        if self.split_half is not None:
            self.empirical_map_split2 = self._apply_mask(self.empirical_map_split2)
            self.empirical_map_split3 = self._apply_mask(self.empirical_map_split3)
    
    def apply_transform_polarangle(self):
        """Transform the polar angle values in the empirical and predicted maps."""
        if self.retinotopic_map == 'polarAngle':
            self.empirical_map = self._transform_polarangle(self.empirical_map)
            self.predicted_map = self._transform_polarangle(self.predicted_map)
            if self.split_half is not None:
                self.empirical_map_split2 = self._transform_polarangle(self.empirical_map_split2)
                self.empirical_map_split3 = self._transform_polarangle(self.empirical_map_split3)
        else:
            raise ValueError("Polar angle transformation is only applicable to this map.")
    
    def convert_to_radian(self):
        """Convert polar angle values in the empirical and predicted maps from degrees to radians."""
        if self.retinotopic_map == 'polarAngle' or self.retinotopic_map == 'eccentricity':
            self.empirical_map = self._convert_to_radian(self.empirical_map)
            self.predicted_map = self._convert_to_radian(self.predicted_map)
            if self.split_half is not None:
                self.empirical_map_split2 = self._convert_to_radian(self.empirical_map_split2)
                self.empirical_map_split3 = self._convert_to_radian(self.empirical_map_split3)
        else:
            raise ValueError("Conversion to radians is only applicable to polar angle and eccentricity maps.")
    
    def binarize_curvature_map(self):
        """Binarize the curvature map for background."""
        background = self.curvature
        nocurv = np.isnan(background)
        background[nocurv == 1] = 0
        background[background < 0] = 0
        background[background > 0] = 1
        return background

    def plot_maps(self, plot_type='predicted', surface_template_path=None,  region_of_interest=None, save_html = False):
        """Plot the specified map type."""
        import matplotlib.pyplot as plt
        from nilearn import plotting
        
        # Background of the surface plot
        background = self.binarize_curvature_map()
        threshold = 1  # threshold for the curvature map

        if self.retinotopic_map == 'polarAngle': 
            self.apply_transform_polarangle()
            self.predicted_map = self.predicted_map + threshold
            self.empirical_map = self.empirical_map + threshold
            max_value = 360 + threshold
        elif self.retinotopic_map == 'eccentricity':
            self.predicted_map = self.predicted_map + threshold
            self.empirical_map = self.empirical_map + threshold
            max_value = 8 + threshold
        elif self.retinotopic_map == 'pRFsize':
            self.predicted_map = self.predicted_map + threshold
            self.empirical_map = self.empirical_map + threshold
            max_value = 2 + threshold

        # Apply mask
        if region_of_interest == 'visualcortex':
            label_primary_visual_areas = ['ROI']
            final_mask_L, final_mask_R, index_L_mask, index_R_mask = roi(
                label_primary_visual_areas)
        elif region_of_interest == 'earlyvisualcortex':
            label_primary_visual_areas = ['ROI']
            final_mask_L, final_mask_R, index_L_mask, index_R_mask = roi_earlyvisualcortex(
                label_primary_visual_areas)
            
        if self.hemisphere == 'lh':
            self.mask = final_mask_L
        elif self.hemisphere == 'rh':
            self.mask = final_mask_R
        
        self.predicted_map[self.mask == 0] = 0
        self.empirical_map[self.mask == 0] = 0
        self.variance_explained[self.mask == 0] = 0

        if plot_type == 'predicted':
            data = self.predicted_map
        elif plot_type == 'empirical':
            data = self.empirical_map
        elif plot_type == 'variance_explained':
            data = self.variance_explained
        else:
            raise ValueError("Invalid plot type specified. Choose from 'predicted', 'empirical', or 'variance_explained'.")

        colour = 'gist_rainbow_r'  # Default colormap
        if self.hemisphere == 'lh':
            surface = osp.join(surface_template_path,'fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii')
            colour = 'gist_rainbow_r'
        elif self.hemisphere == 'rh':
            surface = osp.join(surface_template_path,'fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii')
            if self.retinotopic_map == 'polarAngle':
                colour = 'gist_rainbow'

        view = plotting.view_surf(
            surf_mesh=surface,
            surf_map=np.reshape(data[0:32492], (-1)), bg_map=background,
            cmap=colour, black_bg=False, symmetric_cmap=False,
            threshold=threshold, vmax=max_value)
        if save_html:
            save_path = os.path.join(osp.dirname(osp.abspath(__file__)), '..', 'output/benchmarking/maps')
            if not osp.exists(save_path):
                os.makedirs(save_path) 
            if plot_type == 'predicted':
                view.save_as_html(osp.join(save_path, f"{self.model}_{self.retinotopic_map}_predicted_{self.subject_id}_{self.hemisphere}.html"))
            elif plot_type == 'empirical':
                view.save_as_html(osp.join(save_path, f"empirical_{self.retinotopic_map}_empirical_{self.subject_id}_{self.hemisphere}.html"))
        return view

class RetinotopyData_logbar(RetinotopyData):
    def __init__(self, path, subject_id, hemisphere, retinotopic_map, 
                 number_hemi_nodes=int(32492), experiment=None):
        # Use the experiment parameter to differentiate logbar data
        self.experiment = experiment
        # Initialize the parent class with the provided parameters
        super().__init__(path, subject_id, hemisphere, 
                         retinotopic_map, number_hemi_nodes)
        # Load maps during initialization
        self.predicted_map = self._load_map("predicted")
        self.empirical_map = self._load_map("empirical")
        self.variance_explained = self._load_map("variance_explained")

    def _load_map(self, map_type):
        """Override to load logbar data."""
        if map_type == 'predicted':
            file_name = f"deepRetinotopy/{self.subject_id}.fs_predicted_{self.retinotopic_map}_{self.hemisphere}_curvatureFeat_model.func.gii"
        elif map_type == 'empirical':
            file_name = f"surf/{self.subject_id}.fs_empirical_{self.retinotopic_map}_{self.experiment}_{self.hemisphere}.func.gii"
        elif map_type == 'variance_explained':
            file_name = f"surf/{self.subject_id}.fs_empirical_variance_explained_{self.experiment}_{self.hemisphere}.func.gii"
        else:
            raise ValueError("Invalid map type specified.")

        file_path = osp.join(self.path, self.subject_id, file_name)
        data = np.array(nib.load(file_path).agg_data()).reshape(self.number_hemi_nodes, -1)
        return data

class RetinotopyData_training(RetinotopyData):
    def __init__(self, path, subject_id, hemisphere, ROI_masked, retinotopic_map, 
                 number_hemi_nodes=int(32492), encoding_model='AnalyzePRF'):
        self.path = path
        self.subject_id = subject_id
        self.hemisphere = hemisphere
        self.retinotopic_map = retinotopic_map
        self.number_hemi_nodes = number_hemi_nodes
        self.encoding_model = encoding_model

        # Load maps during initialization
        self.empirical_map = self._load_map("empirical")
        self.variance_explained = self._load_map("variance_explained")

    def _load_map(self, map_type):
        if map_type == 'empirical':
            if self.encoding_model == 'AnalyzePRF':
                file_name = f"surf/{self.subject_id}.fs_empirical_{self.retinotopic_map}_{self.hemisphere}.func.gii"
            elif self.encoding_model == 'SamSrf':
                file_name = f"samsrf/{self.subject_id}.fs_empirical_samsrf_{self.retinotopic_map}_{self.hemisphere}.func.gii"
            else:
                raise ValueError("Invalid encoding model specified.")
        elif map_type == 'variance_explained':
            if self.encoding_model == 'AnalyzePRF':
                file_name = f"surf/{self.subject_id}.fs_empirical_variance_explained_{self.hemisphere}.func.gii"
            elif self.encoding_model == 'SamSrf':
                file_name = f"samsrf/{self.subject_id}.fs_empirical_samsrf_variance_explained_{self.hemisphere}.func.gii"
            else:
                raise ValueError("Invalid encoding model specified.")

        file_path = osp.join(self.path, self.subject_id, file_name)
        data = np.array(nib.load(file_path).agg_data()).reshape(self.number_hemi_nodes, -1)
        return data
    
    def apply_mask_to_maps(self, mask):
        """Apply the mask to the predicted, empirical, and variance explained maps."""
        self.mask = mask
        self.empirical_map = self._apply_mask(self.empirical_map)
        self.variance_explained = self._apply_mask(self.variance_explained)