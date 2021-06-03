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
    Code for continuity error-checking
    
\*---------------------------------------------------------------------------*/
"""

contErr = volScalarField.div(phi)

sumLocalContErr = deltaT * np.sum(np.absolute(contErr.data.cellValues_[0]) * np.average(mesh.geometryData.V_))

globalContErr = deltaT * np.sum(contErr.data.cellValues_[0] * np.average(mesh.geometryData.V_))

cumulativeContErr += globalContErr

print("Time step continuity errors: sum local = " + str(sumLocalContErr) + ", global = " + str(globalContErr) + ", cumulative = " + str(cumulativeContErr))
