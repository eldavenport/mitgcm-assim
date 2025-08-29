# MITgcm Data Assimilation Utilities

This package provides utilities for loading and analyzing MITgcm data assimilation files from TPOSE6 velocity data assimilation project.

## Installation

Install the package in development mode:

```bash
cd mitgcm-assim
pip install -e .
```

## Usage

### Loading Cost Function Data

```python
import mitgcm_assim.cost as cost

# Read all costfunction* files in current directory
ds = cost.read_costfunction('.')

# Access total cost for all iterations
fc_cost = ds.cost.sel(dataset='fc')
print(fc_cost.values)

# Access specific dataset costs
tao_u_cost = ds.cost.sel(dataset='TAO_WO_2012_ADCP_prof_U')
print(f"U-velocity cost: {tao_u_cost.values}")
```

### Loading Control Variables and Sensitivities

```python
import mitgcm_assim.ctrls as ctrls

# Load all default variables
ds = ctrls.read_controls_and_sensitivities('./data', './grid', 0)

# Load specific variables only
ds = ctrls.read_controls_and_sensitivities(
    './data', './grid', 0,
    control_vars=['xx_theta', 'xx_salt'],
    sensitivity_vars=['adxx_theta', 'adxx_salt']
)

# Load a single control variable
da = ctrls.read_single_control('./data', 'xx_theta', 0)

# Check what variables are available
available = ctrls.get_available_controls('./data', 0)
print("Available controls:", available['controls'])
```

## Modules

### `mitgcm_assim.cost`

Functions for reading MITgcm cost function files into xarray Datasets:

- `read_costfunction(directory)`: Read all costfunction* files in a directory
- `parse_costfunction_file(filepath)`: Parse a single costfunction file

### `mitgcm_assim.ctrls` 

Functions for reading MITgcm control variables and sensitivities:

- `read_controls_and_sensitivities()`: Load control variables and sensitivities with grid coordinates
- `read_single_control()`: Load a single control variable
- `get_available_controls()`: List available control files in a directory

## Requirements

- numpy
- xarray  
- MITgcmutils
- xmitgcm