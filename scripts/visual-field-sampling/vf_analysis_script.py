#!/usr/bin/env python3
"""
Visual Field Surface Area Analysis Script
=========================================

This script analyzes surface area in the primary visual cortex dedicated to 
representing different parts of the visual field using retinotopic mapping data.

The script performs three main steps:
1. Resample ROI templates to individual native space
2. Calculate surface area metrics for different visual field regions
3. Save results to a spreadsheet for further analysis

Usage:
    python vf_analysis_script.py --config config.yaml
    or
    python vf_analysis_script.py --freesurfer_dir /path/to/freesurfer --template_dir /path/to/templates --hcp_dir /path/to/hcp --output_file results.csv

Requirements:
    - nibabel
    - numpy
    - subprocess (for wb_command calls)
    - argparse
    - yaml (optional, for config file)
"""

import os
import subprocess
import argparse
import numpy as np
import nibabel as nib
from pathlib import Path
import sys
from typing import Dict, List, Tuple, Optional
from multiprocessing import Pool, cpu_count
from functools import partial
import time
import csv


def check_dependencies():
    """Check if required dependencies are available"""
    try:
        import nibabel
        import numpy
    except ImportError as e:
        print(f"Error: Missing required Python package: {e}")
        sys.exit(1)
    
    # Check if wb_command is available
    try:
        subprocess.run(['wb_command', '-version'], 
                      capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Error: wb_command not found. Please install Connectome Workbench.")
        sys.exit(1)


def resample_roi_to_native(freesurfer_dir: str, subject_id: str, 
                          template_dir: str, hcp_dir: str, verbose: bool = True) -> Dict:
    """
    Resample ROI templates to individual native space using wb_command.
    
    Parameters:
    -----------
    freesurfer_dir : str
        Path to FreeSurfer derivatives directory
    subject_id : str
        Subject identifier
    template_dir : str
        Directory containing ROI template files
    hcp_dir : str
        Directory containing HCP template files
    verbose : bool
        Whether to print progress messages
        
    Returns:
    --------
    dict
        Results dictionary with subject_id, success status, and any error messages
    """
    
    result = {
        'subject_id': subject_id,
        'success': True,
        'errors': []
    }
    
    try:
        subject_path = Path(freesurfer_dir) / subject_id
        deepret_path = subject_path / "deepRetinotopy"
        
        # Create deepRetinotopy directory if it doesn't exist
        deepret_path.mkdir(exist_ok=True)
        
        hemispheres = ['lh', 'rh']
        portions = ['ventral', 'dorsal']
        
        for hemisphere in hemispheres:
            hemi = 'L' if hemisphere == 'lh' else 'R'
            
            source_sphere = Path(hcp_dir) / f"fs_LR-deformed_to-fsaverage.{hemi}.sphere.32k_fs_LR.surf.gii"
            target_sphere = subject_path / "surf" / f"{hemisphere}.sphere.reg.surf.gii"
            source_surface = subject_path / "surf" / f"{subject_id}.{hemisphere}.midthickness.32k_fs_LR.surf.gii"
            target_surface = subject_path / "surf" / f"{hemisphere}.midthickness.surf.gii"

            # Resample dorsal and ventral ROIs
            for portion in portions:
                # Define paths
                template_file = Path(template_dir) / f"V12{portion}.{hemisphere}.32k_fs_LR.label.gii"
                output_file = deepret_path / f"{subject_id}.V12{portion}-roi.{hemisphere}.native.label.gii"
                
                # Check if required input files exist
                required_files = [template_file, source_sphere, target_sphere, source_surface, target_surface]
                missing_files = [f for f in required_files if not f.exists()]
                
                if missing_files:
                    error_msg = f"Missing required files for {subject_id} {hemisphere} {portion}: {missing_files}"
                    result['errors'].append(error_msg)
                    if verbose:
                        print(f"Warning: {error_msg}")
                    continue
                
                # Skip if output already exists
                if output_file.exists():
                    if verbose:
                        print(f"Output already exists for {subject_id} {hemisphere} {portion}, skipping")
                    continue
                
                # Build wb_command
                cmd = [
                    'wb_command', '-label-resample',
                    str(template_file),
                    str(source_sphere),
                    str(target_sphere),
                    'ADAP_BARY_AREA',
                    str(output_file),
                    '-area-surfs',
                    str(source_surface),
                    str(target_surface)
                ]
                
                try:
                    subprocess.run(cmd, check=True, capture_output=True)
                    if verbose:
                        print(f"Successfully resampled {portion} ROI for {subject_id} {hemisphere}")
                except subprocess.CalledProcessError as e:
                    error_msg = f"wb_command failed for {subject_id} {hemisphere} {portion}: {e}"
                    result['errors'].append(error_msg)
                    if verbose:
                        print(f"Error: {error_msg}")
            
            # Resample V1-3 ROI
            ################################################################################
            # Define paths
            template_file = Path(template_dir) / f"V123plusfovea.{hemisphere}.32k_fs_LR.label.gii"
            output_file = deepret_path / f"{subject_id}.V123plusfovea-roi.{hemisphere}.native.label.gii"
            
            # Check if required input files exist
            required_files = [template_file, source_sphere, target_sphere, source_surface, target_surface]
            missing_files = [f for f in required_files if not f.exists()]
            
            if missing_files:
                error_msg = f"Missing required files for {subject_id} {hemisphere} V123plusfovea: {missing_files}"
                result['errors'].append(error_msg)
                if verbose:
                    print(f"Warning: {error_msg}")
                continue
            
            # Skip if output already exists
            if output_file.exists():
                if verbose:
                    print(f"Output already exists for {subject_id} {hemisphere} V123plusfovea, skipping")
                continue
            
            # Build wb_command
            cmd = [
                'wb_command', '-label-resample',
                str(template_file),
                str(source_sphere),
                str(target_sphere),
                'ADAP_BARY_AREA',
                str(output_file),
                '-area-surfs',
                str(source_surface),
                str(target_surface)
            ]
            
            try:
                subprocess.run(cmd, check=True, capture_output=True)
                if verbose:
                    print(f"Successfully resampled V1-3 ROI for {subject_id} {hemisphere}")
            except subprocess.CalledProcessError as e:
                error_msg = f"wb_command failed for {subject_id} {hemisphere} V123plusfovea: {e}"
                result['errors'].append(error_msg)
                if verbose:
                    print(f"Error: {error_msg}")
                    
    except Exception as e:
        result['success'] = False
        result['errors'].append(f"General error processing {subject_id}: {str(e)}")
        if verbose:
            print(f"Error processing {subject_id}: {e}")
    
    if result['errors']:
        result['success'] = False
    
    return result


def resample_roi_wrapper(args):
    """Wrapper function for parallel processing of ROI resampling."""
    subject_id, freesurfer_dir, template_dir, hcp_dir = args
    return resample_roi_to_native(freesurfer_dir, subject_id, template_dir, hcp_dir, verbose=False)


def resample_all_subjects_parallel(freesurfer_dir: str, subjects: List[str], 
                                 template_dir: str, hcp_dir: str, n_jobs: int = None) -> List[Dict]:
    """
    Resample ROI templates for all subjects in parallel.
    
    Parameters:
    -----------
    freesurfer_dir : str
        Path to FreeSurfer derivatives directory
    subjects : list
        List of subject identifiers
    template_dir : str
        Directory containing ROI template files
    hcp_dir : str
        Directory containing HCP template files
    n_jobs : int
        Number of parallel processes to use. If None, uses all available CPUs.
        
    Returns:
    --------
    list
        List of result dictionaries from each resampling operation
    """
    
    if n_jobs is None:
        n_jobs = cpu_count()
    
    print(f"Starting parallel resampling with {n_jobs} processes for {len(subjects)} subjects...")
    
    # Prepare arguments for parallel processing
    args_list = [(subject, freesurfer_dir, template_dir, hcp_dir) for subject in subjects]
    
    start_time = time.time()
    
    # Run parallel processing
    with Pool(processes=n_jobs) as pool:
        results = pool.map(resample_roi_wrapper, args_list)
    
    end_time = time.time()
    
    # Summary
    successful = sum(1 for r in results if r['success'])
    failed = len(results) - successful
    
    print(f"\nResampling completed in {end_time - start_time:.1f} seconds")
    print(f"Successful: {successful}/{len(results)}")
    if failed > 0:
        print(f"Failed: {failed}/{len(results)}")
        print("Failed subjects:")
        for result in results:
            if not result['success']:
                print(f"  - {result['subject_id']}: {'; '.join(result['errors'])}")
    
    return results


def compute_surface_area(vertices: np.ndarray, faces: np.ndarray, mask: np.ndarray) -> float:
    """
    Compute surface area for masked region.
    
    Parameters:
    -----------
    vertices : np.ndarray
        Array of vertex coordinates
    faces : np.ndarray
        Array of triangle face indices
    mask : np.ndarray
        Boolean mask indicating which vertices to include
        
    Returns:
    --------
    float
        Total surface area of masked region
    """
    # Get triangles that have at least one vertex in the mask
    triangles = vertices[faces]
    
    # Calculate triangle areas using cross product
    v1 = triangles[:, 1] - triangles[:, 0]
    v2 = triangles[:, 2] - triangles[:, 0]
    areas = 0.5 * np.linalg.norm(np.cross(v1, v2), axis=1)
    
    # Only include triangles where all vertices are in the mask
    valid_faces = np.all(mask[faces], axis=1)
    return np.sum(areas[valid_faces])


def generate_masks(freesurfer_directory: str, subject_id: str, hemisphere: str, 
                  retinotopic_mapping: str = 'deepRetinotopy', group: str = 'all',
                  wedge_size: int = 35) -> Tuple[np.ndarray, ...]:
    """
    Generate masks for wedges of specified size.
    
    Parameters:
    -----------
    freesurfer_directory : str
        Path to FreeSurfer directory
    subject_id : str
        Subject identifier
    hemisphere : str
        Hemisphere ('lh' or 'rh')
    retinotopic_mapping : str
        Type of retinotopic mapping ('deepRetinotopy')
    wedge_size : int
        Size of visual field wedge in degrees
        
    Returns:
    --------
    tuple
        Masks for upper, lower, horizontal lower, horizontal upper visual field regions,
        plus polar angle data and V1 ROI
    """
    
    # Load retinotopic maps based on mapping type
    if retinotopic_mapping == 'deepRetinotopy':
        polar_angle = nib.load(f'{freesurfer_directory}/{subject_id}/deepRetinotopy/{subject_id}.predicted_polarAngle_model.{hemisphere}.native.func.gii').darrays[0].data
        eccentricity = nib.load(f'{freesurfer_directory}/{subject_id}/deepRetinotopy/{subject_id}.predicted_eccentricity_model.{hemisphere}.native.func.gii').darrays[0].data
    elif retinotopic_mapping == 'empirical':
        # polar_angle = nib.load(f'{freesurfer_directory}/../prfanalyze-vista/{group}/{subject_id}/{hemisphere}.angle_0-360_transformed.gii').darrays[0].data
        # eccentricity = nib.load(f'{freesurfer_directory}/../prfanalyze-vista/{group}/{subject_id}/bayesian_retinotopy/{hemisphere}.inferred_eccen.gii').darrays[0].data
        polar_angle = nib.load(f'{freesurfer_directory}/{subject_id}/deepRetinotopy/{subject_id}.empirical_polarAngle.{hemisphere}.native.func.gii').darrays[0].data
        eccentricity = nib.load(f'{freesurfer_directory}/{subject_id}/deepRetinotopy/{subject_id}.empirical_eccentricity.{hemisphere}.native.func.gii').darrays[0].data
    else:
        raise ValueError(f"Unknown retinotopic mapping type: {retinotopic_mapping}")
    
    # Load ROI masks
    dorsal_roi = nib.load(f'{freesurfer_directory}/{subject_id}/deepRetinotopy/{subject_id}.V12dorsal-roi.{hemisphere}.native.label.gii').darrays[0].data == 1
    ventral_roi = nib.load(f'{freesurfer_directory}/{subject_id}/deepRetinotopy/{subject_id}.V12ventral-roi.{hemisphere}.native.label.gii').darrays[0].data == 1

    v123_roi = nib.load(f'{freesurfer_directory}/{subject_id}/deepRetinotopy/{subject_id}.V123plusfovea-roi.{hemisphere}.native.label.gii').darrays[0].data
    v1_roi = v123_roi == 1
    v2_roi = v123_roi == 2
    v1_plus = v1_roi | v2_roi

    # Define visual field masks
    upper_vf_mask = v1_plus & (polar_angle <= 90 + wedge_size) & (polar_angle > 90 - wedge_size) & ventral_roi
    lower_vf_mask = v1_plus & (polar_angle <= 270 + wedge_size) & (polar_angle > 270 - wedge_size) & dorsal_roi
    
    if hemisphere == 'lh':
        horizontal_upper_vf_mask = v1_plus & (polar_angle > 0) & (polar_angle <= 0 + wedge_size) & v1_roi
        horizontal_lower_vf_mask = v1_plus & (polar_angle < 360) & (polar_angle >= 360 - wedge_size) & v1_roi
    elif hemisphere == 'rh':
        horizontal_upper_vf_mask = v1_plus & (polar_angle <= 180 + wedge_size) & (polar_angle >= 180) & v1_roi
        horizontal_lower_vf_mask = v1_plus & (polar_angle >= 180 - wedge_size) & (polar_angle <= 180) & v1_roi

    # Apply eccentricity constraints
    if eccentricity is not None:
        lower_ecc_thr = 0
        upper_ecc_thr = 6
        ecc_mask = (eccentricity >= lower_ecc_thr) & (eccentricity <= upper_ecc_thr)
        upper_vf_mask = upper_vf_mask & ecc_mask
        lower_vf_mask = lower_vf_mask & ecc_mask
        horizontal_lower_vf_mask = horizontal_lower_vf_mask & ecc_mask
        horizontal_upper_vf_mask = horizontal_upper_vf_mask & ecc_mask
        
    # Mask out invalid data points
    valid_data = ~np.isnan(polar_angle) & ~np.isnan(eccentricity)
    upper_vm_mask = upper_vf_mask & valid_data
    lower_vm_mask = lower_vf_mask & valid_data
    horizontal_lower_mask = horizontal_lower_vf_mask & valid_data
    horizontal_upper_mask = horizontal_upper_vf_mask & valid_data

    return upper_vm_mask, lower_vm_mask, horizontal_lower_mask, horizontal_upper_mask, polar_angle, v1_roi


def analyze_subject_vf_areas(freesurfer_directory: str, subject: str, 
                           retinotopic_mapping: str = 'deepRetinotopy', 
                           wedge_size: int = 35, 
                           group: str = 'all',
                           hemispheres: str = 'both') -> Dict:
    """
    Analyze visual field areas for one subject across hemispheres.
    
    Parameters:
    -----------
    freesurfer_directory : str
        Path to FreeSurfer directory
    subject : str
        Subject identifier
    retinotopic_mapping : str
        Type of retinotopic mapping to use
    wedge_size : int
        Size of visual field wedge in degrees
    hemispheres : str or list
        Hemispheres to analyze ('both', 'lh', 'rh', or list)
        
    Returns:
    --------
    dict
        Dictionary containing surface area measurements and derived metrics
    """
    
    # Initialize totals
    total_upper_area = 0
    total_lower_area = 0
    total_horizontal_area = 0
    total_v1_area = 0
    right_horizontal_area = 0
    left_horizontal_area = 0
    
    # Process hemispheres
    if hemispheres == 'both':
        hemispheres = ['lh', 'rh']
    elif isinstance(hemispheres, str):
        hemispheres = [hemispheres]

    for hemisphere in hemispheres:
        # Load surface geometry
        surface = nib.load(f'{freesurfer_directory}/{subject}/surf/{hemisphere}.midthickness.surf.gii')
        vertices = surface.darrays[0].data  # vertex coordinates
        faces = surface.darrays[1].data     # triangle faces

        upper_vm_mask, lower_vm_mask, horizontal_l_mask, horizontal_u_mask, _, v1_roi = generate_masks(
            freesurfer_directory, subject, hemisphere, retinotopic_mapping, group,
            wedge_size=wedge_size)

        # Calculate surface areas
        upper_area_hemi = compute_surface_area(vertices, faces, upper_vm_mask)
        lower_area_hemi = compute_surface_area(vertices, faces, lower_vm_mask)
        horizontal_upper_hemi = compute_surface_area(vertices, faces, horizontal_u_mask)
        horizontal_lower_hemi = compute_surface_area(vertices, faces, horizontal_l_mask)
        v1_area_hemi = compute_surface_area(vertices, faces, v1_roi)
        
        # Save hemisphere-specific data
        if hemisphere == 'lh':
            left_horizontal_area = horizontal_lower_hemi + horizontal_upper_hemi
        elif hemisphere == 'rh':
            right_horizontal_area = horizontal_lower_hemi + horizontal_upper_hemi

        # Add to totals
        total_upper_area += upper_area_hemi / 2  # divide by 2 to account for V1 and V2
        total_lower_area += lower_area_hemi / 2  # divide by 2 to account for V1 and V2
        total_horizontal_area += horizontal_lower_hemi + horizontal_upper_hemi
        total_v1_area += v1_area_hemi

    # Calculate derived metrics
    mean_area_lower_upper = (total_lower_area + total_upper_area) / 2
    cortical_vma = ((total_lower_area - total_upper_area) / mean_area_lower_upper) * 100 if mean_area_lower_upper > 0 else np.nan
    
    total_vertical_area = total_upper_area + total_lower_area
    mean_area_vertical_horizontal = (total_vertical_area + total_horizontal_area) / 2
    cortical_hva = ((total_horizontal_area - total_vertical_area) / mean_area_vertical_horizontal) * 100 if mean_area_vertical_horizontal > 0 else np.nan
    
    return {
        'subject': subject,
        'upper_area': total_upper_area,
        'lower_area': total_lower_area,
        'vertical_area': total_upper_area + total_lower_area,
        'horizontal_area': total_horizontal_area,
        'total_v1_area': total_v1_area,
        'ratio': total_upper_area/total_lower_area if total_lower_area > 0 else np.nan,
        'cortical_vma': cortical_vma,
        'cortical_hva': cortical_hva,
        'upper_pct_v1': (total_upper_area/total_v1_area)*100 if total_v1_area > 0 else np.nan,
        'lower_pct_v1': (total_lower_area/total_v1_area)*100 if total_v1_area > 0 else np.nan,
        'right_horizontal_area': right_horizontal_area,
        'left_horizontal_area': left_horizontal_area
    }


def save_results_to_csv(results: List[Dict], output_file: str) -> None:
    """
    Save results to a CSV file using the csv module.
    
    Parameters:
    -----------
    results : list
        List of result dictionaries
    output_file : str
        Path to output CSV file
    """
    if not results:
        print("No results to save.")
        return
    
    # Get all unique keys from all results
    fieldnames = set()
    for result in results:
        fieldnames.update(result.keys())
    fieldnames = sorted(fieldnames)
    
    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)


def get_subjects_from_directory(freesurfer_dir: str) -> List[str]:
    """
    Get subject list from FreeSurfer directory.
    
    Parameters:
    -----------
    freesurfer_dir : str
        Path to FreeSurfer directory containing subject folders
        
    Returns:
    --------
    list
        List of subject identifiers
    """
    
    freesurfer_path = Path(freesurfer_dir)
    if not freesurfer_path.exists():
        raise ValueError(f"FreeSurfer directory does not exist: {freesurfer_dir}")
    
    # Get all directories that look like subject folders
    subjects = []
    for item in freesurfer_path.iterdir():
        if item.is_dir() and not item.name.startswith('.'):
            # Check if it has required FreeSurfer structure (surf directory)
            surf_dir = item / 'surf'
            if surf_dir.exists():
                subjects.append(item.name)
    
    return sorted(subjects)


def main():
    """Main function to run the analysis pipeline."""
    
    parser = argparse.ArgumentParser(
        description='Analyze visual field surface areas from retinotopic mapping data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    # Basic usage:
    python vf_analysis_script.py --freesurfer_dir /path/to/freesurfer --template_dir /path/to/templates --hcp_dir /path/to/hcp --output_file results.csv

    # Analyzing multiple methods and groups (space-separated):
    python vf_analysis_script.py --freesurfer_dir /path/to/freesurfer --methods deepRetinotopy empirical --groups adults children --output_file results.csv

    # Or using repeated flags:
    python vf_analysis_script.py --freesurfer_dir /path/to/freesurfer --methods deepRetinotopy --methods empirical --groups adults --groups children
        """)
    
    parser.add_argument('--freesurfer_dir', type=str, required=True, help='Path to FreeSurfer derivatives directory')
    parser.add_argument('--template_dir', type=str, required=True, help='Path to ROI template directory')
    parser.add_argument('--hcp_dir', type=str, required=True, help='Path to HCP template directory')
    parser.add_argument('--output_file', type=str, default='vf_analysis_results.csv', help='Output CSV file')
    parser.add_argument('--wedge_size', type=int, default=45, help='Visual field wedge size in degrees')
    parser.add_argument('--methods', nargs='+', default=['deepRetinotopy'], 
                       help='Retinotopic mapping methods to analyze')
    parser.add_argument('--groups', nargs='+', default=['all'],
                       help='Groups to analyze (e.g., adults children). Use "all" for no group separation')
    parser.add_argument('--n_jobs', type=int, help='Number of parallel processes for resampling (default: all CPUs)')
    parser.add_argument('--skip_resampling', action='store_true', help='Skip ROI resampling step')
    parser.add_argument('--subjects', nargs='+', help='Specific subjects to analyze (optional)')
    
    args = parser.parse_args()
    
    # Check dependencies
    check_dependencies()
    
    # Validate paths exist
    if not Path(args.freesurfer_dir).exists():
        print(f"Error: FreeSurfer directory does not exist: {args.freesurfer_dir}")
        sys.exit(1)
    
    if not args.skip_resampling:
        if not Path(args.template_dir).exists():
            print(f"Error: Template directory does not exist: {args.template_dir}")
            sys.exit(1)
        
        if not Path(args.hcp_dir).exists():
            print(f"Error: HCP directory does not exist: {args.hcp_dir}")
            sys.exit(1)
    
    print("Starting visual field surface area analysis...")
    print(f"FreeSurfer directory: {args.freesurfer_dir}")
    print(f"Output file: {args.output_file}")
    print(f"Methods: {args.methods}")
    print(f"Wedge size: {args.wedge_size} degrees")
    
    # Get subjects
    if args.subjects:
        # Use specific subjects provided
        subjects = args.subjects
        print(f"Using {len(subjects)} specified subjects")
    else:
        # Auto-discover subjects from FreeSurfer directory
        subjects = get_subjects_from_directory(args.freesurfer_dir)
        print(f"Found {len(subjects)} subjects in FreeSurfer directory")
    
    if not subjects:
        print("No subjects found! Please check your FreeSurfer directory or specify subjects manually.")
        sys.exit(1)
    
    # Step 1: Resample ROI templates (if not skipped)
    if not args.skip_resampling:
        print("\nStep 1: Resampling ROI templates to native space...")
        
        n_jobs = args.n_jobs
        if n_jobs is None:
            n_jobs = cpu_count()
            print(f"Using all available CPUs: {n_jobs}")
        else:
            print(f"Using {n_jobs} parallel processes")
        
        resampling_results = resample_all_subjects_parallel(
            args.freesurfer_dir, subjects, 
            args.template_dir, args.hcp_dir, 
            n_jobs=n_jobs)
        
        # Check if any critical failures occurred
        failed_subjects = [r['subject_id'] for r in resampling_results if not r['success']]
        if failed_subjects:
            print(f"\nWarning: {len(failed_subjects)} subjects failed resampling:")
            for subject in failed_subjects[:10]:  # Show first 10
                print(f"  - {subject}")
            if len(failed_subjects) > 10:
                print(f"  ... and {len(failed_subjects) - 10} more")
            
            response = input("\nContinue with surface area analysis? (y/n): ")
            if response.lower() != 'y':
                print("Analysis stopped by user.")
                return
    else:
        print("Skipping ROI resampling step")
    
    # Step 2: Analyze surface areas
    print("\nStep 2: Analyzing surface areas...")
    results = []
    
    total_analyses = len(args.methods) * len(args.groups) * len(subjects)
    current_analysis = 0
    
    for retinotopic_mapping in args.methods:
        for group in args.groups:
            for subject in subjects:
                current_analysis += 1
                print(f"Analyzing {subject} with {retinotopic_mapping} (group: {group}) ({current_analysis}/{total_analyses})")
                
                try:
                    result = analyze_subject_vf_areas(
                        args.freesurfer_dir, subject,
                        retinotopic_mapping=retinotopic_mapping,
                        group=group,
                        wedge_size=args.wedge_size)
                    result['method'] = retinotopic_mapping
                    result['group'] = group
                    
                    # Extract age from subject ID if possible
                    try:
                        result['age'] = int(subject[-2:])
                    except (ValueError, IndexError):
                        result['age'] = 'N/A'
                    
                    results.append(result)
                
                except Exception as e:
                    print(f"Error processing {subject} with {retinotopic_mapping} (group: {group}): {e}")
    
    # Step 3: Save results
    print(f"\nStep 3: Saving {len(results)} results to {args.output_file}...")
    save_results_to_csv(results, args.output_file)
    
    print("Analysis complete!")
    print(f"Results saved to: {args.output_file}")
    
    # Calculate summary statistics
    unique_subjects = set(r['subject'] for r in results)
    print(f"Total subjects analyzed: {len(unique_subjects)}")
    print(f"Total analyses: {len(results)}")


if __name__ == "__main__":
    main()