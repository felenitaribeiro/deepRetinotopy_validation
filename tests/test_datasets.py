import pytest
import numpy as np
import sys
sys.path.append('..')
from functions.datasets import *

@pytest.fixture
def mock_data():
    path = "/home/ribeiro/Projects/deepRetinotopy_validation/HCP/freesurfer/"
    subject_id = "100610"
    hemisphere = "lh"
    ROI_masked = np.ones((32492, 1))
    retinotopic_map = "polarAngle"
    return RetinotopyData(path, subject_id, hemisphere, ROI_masked, retinotopic_map)

@pytest.fixture
def mock_logbar_data():
    path = "/BULK/LABDATA/openneuro/logbar-updated/ds004698/derivatives/freesurfer/"
    subject_id = "sub-01"
    hemisphere = "lh"
    ROI_masked = np.ones((32492, 1))
    retinotopic_map = "polarAngle"
    experiment = "logbar"
    return RetinotopyData_logbar(path, subject_id, hemisphere, ROI_masked, retinotopic_map, experiment=experiment)

@pytest.mark.parametrize("data_class_fixture", ["mock_data", "mock_logbar_data"])
def test_load_map(request, data_class_fixture):
    # Test loading maps for both RetinotopyData and RetinotopyData_logbar
    data_instance = request.getfixturevalue(data_class_fixture)
    try:
        assert isinstance(data_instance.predicted_map, np.ndarray)
        assert isinstance(data_instance.empirical_map, np.ndarray)
        assert isinstance(data_instance.variance_explained, np.ndarray)
    except FileNotFoundError:
        pytest.skip("Mock file not found. Ensure test files are available.")

@pytest.mark.parametrize("data_class_fixture", ["mock_data", "mock_logbar_data"])
def test_apply_mask_to_maps(request, data_class_fixture):
    # Test applying mask to maps for both RetinotopyData and RetinotopyData_logbar
    data_instance = request.getfixturevalue(data_class_fixture)
    mask = np.ones((32492,), dtype=bool)  # Mock mask
    mask[:1000] = False  # Set first 1000 to False
    data_instance.apply_mask_to_maps(mask)
    assert data_instance.predicted_map.shape[0] == np.sum(mask)
    assert data_instance.empirical_map.shape[0] == np.sum(mask)
    assert data_instance.variance_explained.shape[0] == np.sum(mask)

@pytest.mark.parametrize("data_class_fixture", ["mock_data", "mock_logbar_data"])
def test_transform_polarangle(request, data_class_fixture):
    # Test transforming polar angle values
    data_instance = request.getfixturevalue(data_class_fixture)
    data_instance.apply_transform_polarangle()
    mask = data_instance.empirical_map >= 0
    assert np.all(data_instance.empirical_map[mask == 1] >= -1)  # Check for valid polar angle values
    assert np.all(data_instance.predicted_map[mask == 1] >= -1)  # Check for valid polar angle values

@pytest.mark.parametrize("data_class_fixture", ["mock_data", "mock_logbar_data"])
def test_convert_to_radian(request, data_class_fixture):
    data_instance = request.getfixturevalue(data_class_fixture)
    # Test converting polar angle values to radians
    mask = data_instance.predicted_map > -1
    data_instance.apply_mask_to_maps(mask)
    data_instance.convert_to_radian()
    assert np.all(data_instance.empirical_map >= 0)  # Check for valid radian values
    assert np.all(data_instance.predicted_map >= 0)  # Check for valid radian values
    assert np.all(data_instance.empirical_map <= 2 * np.pi)  # Check upper limit for radians
    assert np.all(data_instance.predicted_map <= 2 * np.pi)  # Check upper limit for radians