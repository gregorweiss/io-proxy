# Case setup file for simpleFoam.py

"""--------------------------------------------------------------------------"""

# 2D cavity case
boundaryDefinition_U = {
'movingWall': ( 'fixedValue'     	, np.array([1.0, 0.0, 0.0])), 
'fixedWalls': ( 'fixedValue' 		, np.array([0.0, 0.0, 0.0])),
'frontAndBack': ( 'empty'		    , None   )
}

boundaryDefinition_p = {
'movingWall': ( 'fixedGradient'    	, np.array([0.0])), 
'fixedWalls': ( 'fixedGradient' 	, np.array([0.0])),
'frontAndBack': ( 'empty'		    , None)
}


"""--------------------------------------------------------------------------"""

internalField_U = np.array([0.0, 0.0, 0.0])
internalField_p = np.array([0.0])

# Kinematic viscosity
nu = 0.01

# Time control
startTime = 0.0

writeTime = 1

deltaT = 1

endTime = 50

# Finite volume schemes

ddtScheme = "steadyState"

laplacianScheme = "linearOrthogonal"

divUScheme = "upwind"

gradPScheme = "linear"

divHbyAScheme = "linear"

# Update case setup if it is modified
runTimeModifiable = True

# fvSolution parameters
fvSolution_U = {
'solver'		:	'GaussSeidel',
# 'solver'		:	'PointJacobi',
# 'minIter'		:	0,
# 'maxIter'		:	10,
'tolerance'		:	1e-6,
'relTol'		:	0.01,
'impUR'			:	0.7
}

fvSolution_p = {
'solver'		:	'GaussSeidel',
# 'solver'		:	'PointJacobi',
# 'minIter'		:	0,
# 'maxIter'		:	100000,
'tolerance'		:	1e-5,
'relTol'		:	0.1,
'expUR'			:	0.3,
'refCell'		:	0,
'refValue'		:	np.array([0.0])
}
