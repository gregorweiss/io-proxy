"""
/*---------------------------------------------------------------------------*\
 ####                 ####     #   #             |
 #                    #  #     ##  #             | FoamPython
 ##  #### ####  ##### #### # # #   #### #### ### | v1.0
 #   #  # #  #  # # # #     #  # # #  # #  # # # |
 #   #### ##### # # # #    #    #  #  # #### # # |
-------------------------------------------------------------------------------

Author
	Robert Anderluh, 2021
	
Description
    Setup file for Cythonizing the Gauss-Seidel solver
    
\*---------------------------------------------------------------------------*/
"""

import os

from distutils.core import setup
from Cython.Build import cythonize

directives = {'language_level': 3}

currDir = os.path.dirname(os.path.abspath(__file__))
setup(ext_modules = cythonize(currDir + '/GaussSeidel.pyx', compiler_directives=directives))

"""
// ************************************************************************* //
"""
