"""
Setup script for mitgcm_assim package.

This package provides utilities for loading and analyzing MITgcm data assimilation files.
"""

from setuptools import setup, find_packages
import os

# Read the contents of README file if it exists
this_directory = os.path.abspath(os.path.dirname(__file__))
readme_path = os.path.join(this_directory, 'README.md')
if os.path.exists(readme_path):
    with open(readme_path, encoding='utf-8') as f:
        long_description = f.read()
else:
    long_description = """
    MITgcm Data Assimilation Utilities
    
    This package provides utilities for loading and analyzing MITgcm data assimilation files:
    - load_cost: Functions for reading cost function files
    - load_ctrls: Functions for reading control variables and sensitivities
    """

setup(
    name='mitgcm_assim',
    version='0.1.0',
    author='Ellen Davenport',
    author_email='ellen.davenport7@gmail.com',
    description='Utilities for MITgcm data assimilation analysis',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='',
    packages=find_packages(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    python_requires='>=3.7',
    install_requires=[
        'numpy',
        'xarray',
        'MITgcmutils',
        'xmitgcm',
    ],
    extras_require={
        'dev': [
            'pytest',
            'pytest-cov',
            'flake8',
            'black',
        ],
    },
    entry_points={
        'console_scripts': [
            # Add command line scripts here if needed
        ],
    },
    package_data={
        'mitgcm_assim': [
            # Add any package data files here if needed
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
