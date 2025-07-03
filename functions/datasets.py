import os.path as osp
import numpy as np
import nibabel as nib

class RetinotopyData:
    def __init__(self, path, subject_id, hemisphere, ROI_masked, 
                 retinotopic_map, number_hemi_nodes=int(32492), model_index=None):
        self.path = path
        self.subject_id = subject_id
        self.hemisphere = hemisphere
        self.ROI_masked = ROI_masked
        self.retinotopic_map = retinotopic_map
        self.number_hemi_nodes = number_hemi_nodes
        self.model_index = model_index  # optional model index for specific models

        # Load maps during initialization
        self.predicted_map = self._load_map("predicted")
        self.empirical_map = self._load_map("empirical")
        self.variance_explained = self._load_map("variance_explained")

    def _load_map(self, map_type):
        """Load a retinotopic map (predicted or empirical)."""
        if map_type == 'predicted':
            if self.model_index is not None:
                file_name = f"deepRetinotopy/{self.subject_id}.fs_predicted_{self.retinotopic_map}_{self.hemisphere}_curvatureFeat_model{self.model_index}.func.gii"
            else:
                file_name = f"deepRetinotopy/{self.subject_id}.fs_predicted_{self.retinotopic_map}_{self.hemisphere}_curvatureFeat_model.func.gii"
        elif map_type == 'empirical':
            file_name = f"surf/{self.subject_id}.fs_empirical_{self.retinotopic_map}_{self.hemisphere}.func.gii"
        elif map_type == 'variance_explained':
            file_name = f"surf/{self.subject_id}.fs_empirical_variance_explained_{self.hemisphere}.func.gii"
        else:
            raise ValueError("Invalid map type specified.")

        file_path = osp.join(self.path, self.subject_id, file_name)
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
    
    def apply_transform_polarangle(self):
        """Transform the polar angle values in the empirical and predicted maps."""
        if self.retinotopic_map == 'polarAngle':
            self.empirical_map = self._transform_polarangle(self.empirical_map)
            self.predicted_map = self._transform_polarangle(self.predicted_map)
        else:
            raise ValueError("Polar angle transformation is only applicable to this map.")
    
    def convert_to_radian(self):
        """Convert polar angle values in the empirical and predicted maps from degrees to radians."""
        if self.retinotopic_map == 'polarAngle' or self.retinotopic_map == 'eccentricity':
            self.empirical_map = self._convert_to_radian(self.empirical_map)
            self.predicted_map = self._convert_to_radian(self.predicted_map)
        else:
            raise ValueError("Conversion to radians is only applicable to polar angle and eccentricity maps.")

class RetinotopyData_logbar(RetinotopyData):
    def __init__(self, path, subject_id, hemisphere, ROI_masked, retinotopic_map, 
                 number_hemi_nodes=int(32492), experiment=None):
        # Use the experiment parameter to differentiate logbar data
        self.experiment = experiment
        # Initialize the parent class with the provided parameters
        super().__init__(path, subject_id, hemisphere, ROI_masked, 
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
            file_name = f"surf/{self.subject_id}.fs_empirical_{self.retinotopic_map}_{self.hemisphere}.func.gii"
        elif map_type == 'variance_explained':
            file_name = f"surf/{self.subject_id}.fs_empirical_variance_explained_{self.hemisphere}.func.gii"
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