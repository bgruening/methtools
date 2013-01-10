#!/usr/bin/env python
"""
Python setup file for the mehttools project.
"""
from setuptools import setup, find_packages

__author__ = "Bjoern A. Gruening"
__copyright__ = "Copyright 2013, Bjoern A. Gruening"
__credits__ = ["Bjoern Gruening", "Dr. Ralf Gilsbach"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Bjoern Gruening"
__email__ = "bjoern.gruening@gmail.com"
__status__ = "Development"

setup(name = "methtools",
    version = __version__,
    author = __author__,
    author_email = __email__,
    description = "",
    packages = find_packages(),
    #install_requires = ['foo>=0.4']
    entry_points = {
    'console_scripts': [ 
        #'calling = calling:main',
        #'dmr = dmr:main',
        #'tiling = tiling:main',
        #'filter = filter:main',
        #'destrand = destrand:main',
        #'plot = plot:main',
        'methtools = methtools:main'
        ]
        }
    )
