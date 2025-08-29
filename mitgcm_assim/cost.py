"""
Functions for loading MITgcm cost function files into xarray Dataset.

This module provides functions to read costfunction* files and create xarray Datasets
with proper structure for analysis and visualization.
"""

import numpy as np
import xarray as xr
import glob
import os
import re
from pathlib import Path


def parse_costfunction_file(filepath):
    """
    Parse a single costfunction file and return dictionary of data.
    
    Parameters
    ----------
    filepath : str
        Path to the costfunction file
        
    Returns
    -------
    dict
        Dictionary with dataset names as keys and (cost, n_obs) tuples as values
    """
    data = {}
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            # Skip f_ob* rows
            if line.startswith('f_ob'):
                continue
                
            # Parse fc row (total cost)
            if line.startswith('fc'):
                parts = line.split('=')[1].strip().split()
                cost = float(parts[0].replace('D', 'E'))
                n_obs = float(parts[1].replace('D', 'E'))
                data['fc'] = (cost, n_obs)
                continue
            
            # Parse other rows
            if '=' in line:
                left_part, right_part = line.split('=', 1)
                left_part = left_part.strip()
                right_part = right_part.strip()
                
                # Parse the right side (cost and n_obs values)
                values = right_part.split()
                if len(values) >= 2:
                    cost = float(values[0].replace('D', 'E'))
                    n_obs = float(values[1].replace('D', 'E'))
                    
                    # Parse the left side to get dataset name
                    # Handle cases like:
                    # "NODC_TP_2012_MRB_w             prof_T"
                    # "sst-MW          (gencost  1)"
                    # "xx_atemp        (gentim2d  1)"
                    
                    # Remove parenthetical information
                    left_clean = re.sub(r'\s*\([^)]+\)\s*', '', left_part)
                    left_clean = left_clean.strip()
                    
                    # Split by whitespace
                    parts = left_clean.split()
                    
                    if len(parts) == 1:
                        # Simple case like "sst-MW"
                        dataset_name = parts[0]
                    elif len(parts) == 2:
                        # Case like "NODC_TP_2012_MRB_w prof_T"
                        dataset_name = f"{parts[0]}_{parts[1]}"
                    else:
                        # Fallback: join all parts
                        dataset_name = "_".join(parts)
                    
                    data[dataset_name] = (cost, n_obs)
    
    return data


def read_costfunction(directory='.'):
    """
    Read all costfunction files in a directory into xarray Dataset.
    
    This function automatically looks for all files matching the pattern 'costfunction*'
    and loads them into a single xarray Dataset. The iteration number is extracted
    from the number at the end of each filename.
    
    Parameters
    ----------
    directory : str
        Directory containing costfunction files (default: current directory)
        
    Returns
    -------
    xarray.Dataset
        Dataset with cost function data containing:
        - 'cost' data variable: cost values for each dataset and iteration
        - 'n_obs' data variable: number of observations for each dataset and iteration
        - 'iter' coordinate: iteration numbers
        - 'dataset' coordinate: dataset names
        
    Raises
    ------
    FileNotFoundError
        If no costfunction* files are found in the directory
        
    Examples
    --------
    >>> import mitgcm_assim.load_cost as load_cost
    >>> ds = load_cost.read_costfunction('.')
    >>> print(ds)
    >>> # Access total cost for all iterations
    >>> fc_cost = ds.cost.sel(dataset='fc')
    >>> print(fc_cost)
    """
    # Find all costfunction files
    filepath_pattern = os.path.join(directory, 'costfunction*')
    files = sorted(glob.glob(filepath_pattern))
    
    if not files:
        raise FileNotFoundError(f"No files found matching pattern: {filepath_pattern}")
    
    # Extract iteration numbers from filenames
    iterations = []
    file_data = {}
    
    for filepath in files:
        filename = os.path.basename(filepath)
        # Extract number at end of filename
        match = re.search(r'(\d+)$', filename)
        if match:
            iter_num = int(match.group(1))
        else:
            # If no number found, use 0 as default
            iter_num = 0
        
        iterations.append(iter_num)
        file_data[iter_num] = parse_costfunction_file(filepath)
    
    # Get all unique dataset names across all files
    all_datasets = set()
    for data in file_data.values():
        all_datasets.update(data.keys())
    all_datasets = sorted(list(all_datasets))
    
    # Create arrays for cost and n_obs
    iterations = sorted(iterations)
    n_iters = len(iterations)
    n_datasets = len(all_datasets)
    
    cost_array = np.full((n_iters, n_datasets), np.nan)
    n_obs_array = np.full((n_iters, n_datasets), np.nan)
    
    # Fill arrays
    for i, iter_num in enumerate(iterations):
        data = file_data[iter_num]
        for j, dataset_name in enumerate(all_datasets):
            if dataset_name in data:
                cost_array[i, j] = data[dataset_name][0]
                n_obs_array[i, j] = data[dataset_name][1]
    
    # Create xarray Dataset
    ds = xr.Dataset(
        {
            'cost': (['iter', 'dataset'], cost_array),
            'n_obs': (['iter', 'dataset'], n_obs_array)
        },
        coords={
            'iter': iterations,
            'dataset': all_datasets
        },
        attrs={
            'description': 'MITgcm cost function data',
            'source_files': [os.path.basename(f) for f in files],
            'source_directory': os.path.abspath(directory)
        }
    )
    
    return ds


# Legacy function name for backward compatibility
read_costfunction_files = read_costfunction