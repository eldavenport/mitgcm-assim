"""
MITgcm Data Assimilation Utilities

This package provides utilities for loading and analyzing MITgcm data assimilation files:
- cost: Functions for reading cost function files
- ctrls: Functions for reading control variables and sensitivities
"""

from . import cost
from . import ctrls

__version__ = "0.1.0"
__author__ = "Ellen Davenport"
__description__ = "Utilities for MITgcm data assimilation analysis"
