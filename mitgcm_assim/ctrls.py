"""
Functions for loading MITgcm control variables and sensitivities into xarray DataArrays.

This module provides efficient, reusable functions to combine grid information from
open_mdsdataset with control/sensitivity data from MITgcmutils.rdmds.
"""

import xarray as xr
import numpy as np
try:
    from MITgcmutils import mds
except ImportError:
    # Try alternative import pattern
    import MITgcmutils.mds as mds
from xmitgcm import open_mdsdataset


# Known units for MITgcm control variables
CONTROL_VARIABLE_UNITS = {
    # Ocean velocities
    'xx_uvel': 'm/s',
    'xx_vvel': 'm/s',
    'adxx_uvel': 'm/s',
    'adxx_vvel': 'm/s',
    
    # Wind velocities
    'xx_uwind': 'm/s',
    'xx_vwind': 'm/s', 
    'adxx_uwind': 'm/s',
    'adxx_vwind': 'm/s',
    
    # Temperature
    'xx_theta': 'deg C',
    'xx_atemp': 'K',
    'adxx_theta': 'deg C',
    'adxx_atemp': 'K',
    
    # Salinity
    'xx_salt': 'psu',
    'adxx_salt': 'psu',
    
    # Specific humidity
    'xx_aqh': 'kg/kg',
    'adxx_aqh': 'kg/kg',
    
    # Radiation
    'xx_lwdown': 'W/m^2',
    'xx_swdown': 'W/m^2',
    'adxx_lwdown': 'W/m^2',
    'adxx_swdown': 'W/m^2',
    
    # Precipitation
    'xx_precip': 'm/s',
    'adxx_precip': 'm/s'
}


def _get_variable_units(var_name):
    """
    Get the units for a control variable, including effective versions.
    
    Parameters
    ----------
    var_name : str
        Variable name (e.g., 'xx_theta', 'adxx_salt.effective')
        
    Returns
    -------
    str
        Units string, or 'unknown' if not found
    """
    # For .effective variables, look up the base variable units
    # but keep the original name intact
    if '.effective' in var_name:
        base_name = var_name.replace('.effective', '')
        return CONTROL_VARIABLE_UNITS.get(base_name, 'unknown')
    else:
        return CONTROL_VARIABLE_UNITS.get(var_name, 'unknown')


def read_controls_and_sensitivities(data_dir, grid_dir, iteration, 
                                    control_vars=None, sensitivity_vars=None, 
                                    effective_vars=None, load_grid_ds=True):
    """
    Load MITgcm control variables and sensitivities into an xarray Dataset with proper grid coordinates.
    
    Parameters
    ----------
    data_dir : str
        Directory containing the control/sensitivity files
    grid_dir : str  
        Directory containing grid files
    iteration : int
        Iteration number to load
    control_vars : list of str, optional
        List of control variable names to load (e.g., ['xx_theta', 'xx_salt'])
        If None, loads common control variables
    sensitivity_vars : list of str, optional
        List of sensitivity variable names to load (e.g., ['adxx_theta', 'adxx_salt'])
        If None, loads common sensitivity variables  
    effective_vars : list of str, optional
        List of effective variable names to load (e.g., ['xx_theta.effective'])
        If None, loads effective versions of control_vars
    load_grid_ds : bool, optional
        Whether to load a grid dataset using open_mdsdataset. Default True.
        
    Returns
    -------
    xarray.Dataset
        Dataset containing control variables, sensitivities, and model grid coordinates
        
    Examples
    --------
    >>> import mitgcm_assim.load_ctrls as load_ctrls
    >>> # Load all default variables
    >>> ds = load_ctrls.read_controls_and_sensitivities('./data', './grid', 0)
    >>> 
    >>> # Load specific variables only
    >>> ds = load_ctrls.read_controls_and_sensitivities(
    ...     './data', './grid', 0,
    ...     control_vars=['xx_theta', 'xx_salt'],
    ...     sensitivity_vars=['adxx_theta', 'adxx_salt']
    ... )
    """
    
    # Default variable lists
    if control_vars is None:
        control_vars = ['xx_theta', 'xx_salt', 'xx_uwind', 'xx_vwind', 
                       'xx_atemp', 'xx_aqh', 'xx_lwdown', 'xx_precip', 'xx_swdown','xx_vvel','xx_uvel']
    
    if sensitivity_vars is None:
        sensitivity_vars = ['adxx_theta', 'adxx_salt', 'adxx_uwind', 'adxx_vwind',
                           'adxx_atemp', 'adxx_aqh', 'adxx_lwdown', 'adxx_precip', 'adxx_swdown','adxx_vvel','adxx_uvel']
    
    if effective_vars is None:
        effective_vars = [var + '.effective' for var in control_vars]
    
    # Load grid information
    grid_ds = None
    if load_grid_ds:
        try:
            grid_ds = open_mdsdataset(data_dir, grid_dir=grid_dir, iters=iteration, 
                                    ignore_unknown_vars=True)
        except Exception as e:
            print(f"Warning: Could not load grid dataset: {e}")
            print("Proceeding without grid coordinates...")
    
    data_vars = {}
    coords = {}
    
    # Load all requested variables
    all_vars = control_vars + sensitivity_vars + effective_vars
    
    for var_name in all_vars:
        try:
            # Load data using rdmds
            data, _, meta = mds.rdmds(f'{data_dir}/{var_name}', iteration, returnmeta=True)
            
            # Create DataArray with appropriate coordinates
            da = _create_dataarray_from_rdmds(data, meta, var_name, grid_ds)
            data_vars[var_name] = da
            
            # Don't collect generic coordinates - we'll replace them with grid coordinates
            
        except Exception as e:
            print(f"Warning: Could not load {var_name}: {e}")
            continue
    
    # Add model grid coordinates if available
    if grid_ds is not None:
        try:
            # Add grid coordinates with proper names
            if 'XC' in grid_ds:
                coords['XC'] = grid_ds.XC
            if 'YC' in grid_ds:
                coords['YC'] = grid_ds.YC  
            if 'XG' in grid_ds:
                coords['XG'] = grid_ds.XG
            if 'YG' in grid_ds:
                coords['YG'] = grid_ds.YG
            if 'Z' in grid_ds:
                coords['Z'] = grid_ds.Z
            if 'Zl' in grid_ds:
                coords['Zl'] = grid_ds.Zl
            if 'drF' in grid_ds:
                coords['drF'] = grid_ds.drF
            if 'dxC' in grid_ds:
                coords['dxC'] = grid_ds.dxC
            if 'dyC' in grid_ds:
                coords['dyC'] = grid_ds.dyC
            if 'dxG' in grid_ds:
                coords['dxG'] = grid_ds.dxG
            if 'dyG' in grid_ds:
                coords['dyG'] = grid_ds.dyG
        except Exception as e:
            print(f"Warning: Could not add all grid coordinates: {e}")
    
    # Replace coordinate names in data variables to use proper grid names
    for var_name, da in data_vars.items():
        da = _rename_coordinates_to_grid(da, coords)
        data_vars[var_name] = da
    
    # Create the Dataset
    ds = xr.Dataset(data_vars=data_vars, coords=coords)
    
    return ds


def _rename_coordinates_to_grid(da, grid_coords):
    """
    Replace generic coordinate names with proper MITgcm grid coordinate names.
    
    Parameters
    ----------
    da : xarray.DataArray
        DataArray with generic coordinate names
    grid_coords : dict
        Dictionary of grid coordinates
        
    Returns
    -------
    xarray.DataArray
        DataArray with proper grid coordinates
    """
    
    # Mapping from generic names to grid names
    coord_mapping = {
        'k': 'Z',
        'j': 'YC', 
        'i': 'XC',
        'j_g': 'YG',
        'i_g': 'XG'
    }
    
    # Create new dimension names
    new_dims = []
    for dim in da.dims:
        new_dims.append(coord_mapping.get(dim, dim))
    
    # Create new coordinates dictionary - only use grid coordinates, not generic ones
    new_coords = {}
    for i, dim in enumerate(da.dims):
        if dim in coord_mapping:
            grid_name = coord_mapping[dim]
            if grid_name in grid_coords:
                # Use the actual grid coordinate
                new_coords[grid_name] = grid_coords[grid_name]
        elif dim == 'time':
            # Keep time coordinate
            if 'time' in da.coords:
                new_coords['time'] = da.coords['time']
            else:
                new_coords['time'] = np.arange(da.shape[i])
    
    # Create new DataArray with proper grid dimensions and coordinates
    da_renamed = xr.DataArray(
        da.values,
        dims=new_dims,
        coords=new_coords,
        name=da.name,
        attrs=da.attrs
    )
    
    return da_renamed


def _create_dataarray_from_rdmds(data, meta, var_name, grid_ds=None):
    """
    Create an xarray DataArray from rdmds output with proper coordinates.
    
    Parameters
    ----------
    data : numpy.ndarray
        Data array from rdmds
    meta : dict
        Metadata dictionary from rdmds  
    var_name : str
        Variable name
    grid_ds : xarray.Dataset, optional
        Grid dataset for coordinate information
        
    Returns
    -------
    xarray.DataArray
        DataArray with proper coordinates and attributes
    """
    
    # Get dimensions from metadata
    dims = meta['dimlist'] if 'dimlist' in meta else data.shape
    meta_ndims = meta['ndims'][0] if 'ndims' in meta else len(data.shape)
    nrecords = meta['nrecords'][0] if 'nrecords' in meta else 1
    
    # Use actual data shape rather than metadata dims to avoid conflicts
    actual_shape = data.shape
    actual_ndims = len(actual_shape)
    
    # Determine coordinate names and values based on nrecords and dimensions
    coords = {}
    dim_names = []
    
    if nrecords > 1:
        # Multiple records = time dimension + spatial dimensions
        if actual_ndims == 3:
            # time + 2D space
            nt, ny, nx = actual_shape
            dim_names = ['time', 'j', 'i']
            coords['time'] = np.arange(nt)
        elif actual_ndims == 4:
            # time + 3D space
            nt, nz, ny, nx = actual_shape
            dim_names = ['time', 'k', 'j', 'i']
            coords['time'] = np.arange(nt)
        else:
            # Fallback for other dimensions with time
            dim_names = ['time'] + [f'dim_{i}' for i in range(1, actual_ndims)]
            coords['time'] = np.arange(actual_shape[0])
            for i in range(1, actual_ndims):
                coords[dim_names[i]] = np.arange(actual_shape[i])
    elif actual_ndims == 2:
        # 2D field (likely surface or forcing field)
        ny, nx = actual_shape
        dim_names = ['j', 'i']
        if grid_ds is not None:
            try:
                # Try different grid coordinate possibilities
                if hasattr(grid_ds, 'YC'):
                    if grid_ds.YC.ndim == 2 and grid_ds.YC.shape[0] == ny:
                        coords['j'] = grid_ds.YC[:, 0].values
                    elif grid_ds.YC.ndim == 1 and len(grid_ds.YC) == ny:
                        coords['j'] = grid_ds.YC.values
                    else:
                        coords['j'] = np.arange(ny)
                else:
                    coords['j'] = np.arange(ny)
                    
                if hasattr(grid_ds, 'XC'):
                    if grid_ds.XC.ndim == 2 and grid_ds.XC.shape[1] == nx:
                        coords['i'] = grid_ds.XC[0, :].values
                    elif grid_ds.XC.ndim == 1 and len(grid_ds.XC) == nx:
                        coords['i'] = grid_ds.XC.values
                    else:
                        coords['i'] = np.arange(nx)
                else:
                    coords['i'] = np.arange(nx)
            except:
                coords['j'] = np.arange(ny)
                coords['i'] = np.arange(nx)
        else:
            coords['j'] = np.arange(ny)
            coords['i'] = np.arange(nx)
                
    elif actual_ndims == 3:
        # 3D field without time (depth + space)
        nz, ny, nx = actual_shape
        dim_names = ['k', 'j', 'i']
        
    elif actual_ndims == 4:
        # 4D field without time - unusual case
        nw, nz, ny, nx = actual_shape
        dim_names = ['w', 'k', 'j', 'i']
        coords['w'] = np.arange(nw)
    else:
        # Fallback for other dimensions
        dim_names = [f'dim_{i}' for i in range(actual_ndims)]
        for i, dim_size in enumerate(actual_shape):
            coords[dim_names[i]] = np.arange(dim_size)
    
    # Add grid coordinates for spatial dimensions if not already set
    if grid_ds is not None and nrecords == 1:
        try:
            # Handle depth coordinate for 3D fields
            if 'k' in dim_names:
                k_idx = dim_names.index('k')
                nz = actual_shape[k_idx]
                if hasattr(grid_ds, 'Z') and len(grid_ds.Z) == nz:
                    coords['k'] = grid_ds.Z.values
                else:
                    coords['k'] = np.arange(nz)
                    
            # Handle spatial coordinates
            if 'j' in dim_names and 'i' in dim_names:
                j_idx = dim_names.index('j')
                i_idx = dim_names.index('i')
                ny = actual_shape[j_idx]
                nx = actual_shape[i_idx]
                
                if hasattr(grid_ds, 'YC'):
                    if grid_ds.YC.ndim == 2 and grid_ds.YC.shape[0] == ny:
                        coords['j'] = grid_ds.YC[:, 0].values
                    elif grid_ds.YC.ndim == 1 and len(grid_ds.YC) == ny:
                        coords['j'] = grid_ds.YC.values
                    else:
                        coords['j'] = np.arange(ny)
                else:
                    coords['j'] = np.arange(ny)
                    
                if hasattr(grid_ds, 'XC'):
                    if grid_ds.XC.ndim == 2 and grid_ds.XC.shape[1] == nx:
                        coords['i'] = grid_ds.XC[0, :].values
                    elif grid_ds.XC.ndim == 1 and len(grid_ds.XC) == nx:
                        coords['i'] = grid_ds.XC.values
                    else:
                        coords['i'] = np.arange(nx)
                else:
                    coords['i'] = np.arange(nx)
        except:
            # Fallback to indices if grid coordinates don't work
            pass
    
    # Create the DataArray
    da = xr.DataArray(
        data,
        dims=dim_names,
        coords=coords,
        name=var_name,
        attrs={
            'units': _get_variable_units(var_name),
            'long_name': var_name,
            'source': 'MITgcm rdmds',
            'dataprec': meta.get('dataprec', ['unknown'])[0],
            'timestepnumber': meta.get('timestepnumber', [0])[0]
        }
    )
    
    return da


def read_single_control(data_dir, var_name, iteration, grid_ds=None):
    """
    Load a single control variable into an xarray DataArray.
    
    Parameters
    ---------- 
    data_dir : str
        Directory containing the control file
    var_name : str
        Variable name (e.g., 'xx_theta', 'adxx_salt.effective')
    iteration : int
        Iteration number to load
    grid_ds : xarray.Dataset, optional
        Grid dataset for coordinate information
        
    Returns
    -------
    xarray.DataArray
        DataArray with proper coordinates and attributes
        
    Examples
    --------
    >>> import mitgcm_assim.load_ctrls as load_ctrls
    >>> # Load a single control variable
    >>> da = load_ctrls.read_single_control('./data', 'xx_theta', 0)
    >>> print(da)
    """
    
    # Load data using rdmds
    data, _, meta = mds.rdmds(f'{data_dir}/{var_name}', iteration, returnmeta=True)
    
    # Create DataArray with appropriate coordinates
    da = _create_dataarray_from_rdmds(data, meta, var_name, grid_ds)
    
    return da


def get_available_controls(data_dir, iteration):
    """
    Get list of available control/sensitivity files in a directory.
    
    Parameters
    ----------
    data_dir : str
        Directory to search
    iteration : int
        Iteration number
        
    Returns
    -------
    dict
        Dictionary with keys 'controls', 'sensitivities', 'effective' containing
        lists of available variable names
        
    Examples
    --------
    >>> import mitgcm_assim.load_ctrls as load_ctrls
    >>> available = load_ctrls.get_available_controls('./data', 0)
    >>> print("Available controls:", available['controls'])
    >>> print("Available sensitivities:", available['sensitivities'])
    """
    import os
    import glob
    
    # Look for files matching iteration pattern
    pattern = f"{data_dir}/*{iteration:010d}.data"
    files = glob.glob(pattern)
    
    controls = []
    sensitivities = []
    effective = []
    
    for file in files:
        basename = os.path.basename(file)
        var_name = basename.replace(f'.{iteration:010d}.data', '')
        
        if var_name.startswith('adxx_'):
            if '.effective' in var_name:
                effective.append(var_name)
            else:
                sensitivities.append(var_name)
        elif var_name.startswith('xx_'):
            if '.effective' in var_name:
                effective.append(var_name)  
            else:
                controls.append(var_name)
    
    return {
        'controls': sorted(controls),
        'sensitivities': sorted(sensitivities), 
        'effective': sorted(effective)
    }


# Legacy function names for backward compatibility
load_controls_and_sensitivities = read_controls_and_sensitivities
load_single_control = read_single_control