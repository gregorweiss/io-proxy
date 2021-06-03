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
    General useful functions
    
\*---------------------------------------------------------------------------*/
"""

import os
import psutil
import numpy as np
from src.OpenFOAM.fields.include import *

def printMemoryUsageInMB():
	print("\nMemory usage is: " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2) + ' MB\n')


# def weightedAverage(field, weightField):
	
	# if (np.size(field) != np.size(weightField)):
		# raise RuntimeError("Field sizes in weightedAverage function are not the same!")
	
	# weightFieldAverage = np.average(weightField)
	
	# result = np.average(field * weightField / weightFieldAverage)
	
	# return result



	
