# Case setup file for laplacianFoam.py

"""--------------------------------------------------------------------------"""

# Cavity case
boundaryDefinition_T = {
'movingWall': ( 'fixedValue'     	, np.array([373.0])), 
'fixedWalls': ( 'fixedValue' 		, np.array([273.0])),
'frontAndBack': ( 'empty'		    , None   )
}

"""--------------------------------------------------------------------------"""

internalField_T = np.array([273.0])

# Conductivity
DT = 2e-4

# Time control
startTime = 0.0

writeTime = 0.1

deltaT = 0.005

endTime = 3.0

# Finite volume schemes

ddtScheme = "Euler"

laplacianScheme = "linearOrthogonal"

# Update case setup if it is modified
runTimeModifiable = True

# fvSolution parameters
fvSolution_T = {
'solver'		:	'GaussSeidel',
# 'solver'		:	'PointJacobi',
# 'minIter'		:	0,
# 'maxIter'		:	5,
'tolerance'		:	1e-6,
# 'relTol'		:	0
}
