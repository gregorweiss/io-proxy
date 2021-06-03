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
    Setup file for Cythonizing volVectorField
    
\*---------------------------------------------------------------------------*/
"""

import os

from distutils.core import setup
from Cython.Build import cythonize

directives = {'language_level': 3}

currDir = os.path.dirname(os.path.abspath(__file__))
setup(ext_modules = cythonize(currDir + '/volVectorField.pyx', compiler_directives=directives))

"""
// ************************************************************************* //
"""
