# Case setup file for icoFoam.py

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

# # 3D cavity case
# boundaryDefinition_U = {
# 'movingWall': ( 'fixedValue'     	, np.array([1.0, 0.0, 0.0])), 
# 'fixedWalls': ( 'fixedValue' 		, np.array([0.0, 0.0, 0.0])),
# 'frontAndBack': ( 'fixedValue'     	, np.array([1.0, 0.0, 0.0])), 
# }

# boundaryDefinition_p = {
# 'movingWall': ( 'fixedGradient'    	, np.array([0.0])), 
# 'fixedWalls': ( 'fixedGradient' 	, np.array([0.0])),
# 'frontAndBack': ( 'fixedGradient'    	, np.array([0.0])), 
# }

"""--------------------------------------------------------------------------"""

# # 2D inlet-outlet cavity case
# boundaryDefinition_U = {
# 'inlet':		( 'fixedValue'     		, np.array([0.2, 0.0, 0.0])), 
# # 'inlet':		( 'fixedGradient' 		, np.array([0.0, 0.0, 0.0])),
# 'outlet':		( 'fixedGradient' 		, np.array([0.0, 0.0, 0.0])),
# 'walls':		( 'fixedValue'   		, np.array([0.0, 0.0, 0.0])), 
# 'frontAndBack':	( 'empty'				, None					   )
# }

# boundaryDefinition_p = {
# 'inlet':		( 'fixedGradient' 		, np.array([0.0])),
# # 'inlet':		( 'fixedValue'   		, np.array([1.84])), 
# 'outlet':		( 'fixedValue'   		, np.array([0.0])), 
# 'walls':		( 'fixedGradient' 		, np.array([0.0])), 
# 'frontAndBack':	( 'empty'			    , None			 )
# }

"""--------------------------------------------------------------------------"""

# # 3D inlet-outlet cavity case
# boundaryDefinition_U = {
# 'inlet':		( 'fixedValue'     		, np.array([0.2, 0.0, 0.0])),
# 'outlet':		( 'fixedGradient' 		, np.array([0.0, 0.0, 0.0])),
# 'walls':		( 'fixedValue'   		, np.array([0.0, 0.0, 0.0])), 
# }

# boundaryDefinition_p = {
# 'inlet':		( 'fixedGradient' 		, np.array([0.0])),
# 'outlet':		( 'fixedValue'   		, np.array([0.0])), 
# 'walls':		( 'fixedGradient' 		, np.array([0.0])), 
# }

"""--------------------------------------------------------------------------"""

# # Pitz Daily
# boundaryDefinition_U = {
# 'inlet':		( 'fixedValue'     		, np.array([1.0, 0.0, 0.0])), 
# 'outlet':		( 'fixedGradient' 		, np.array([0.0, 0.0, 0.0])),
# 'upperWall':	( 'fixedValue'   		, np.array([0.0, 0.0, 0.0])), 
# 'lowerWall':	( 'fixedValue'   		, np.array([0.0, 0.0, 0.0])), 
# 'frontAndBack':	( 'empty'				, None					   )
# }

# boundaryDefinition_p = {
# 'inlet':		( 'fixedGradient' 		, np.array([0.0])),
# 'outlet':		( 'fixedValue'   		, np.array([0.0])), 
# 'upperWall':	( 'fixedGradient' 		, np.array([0.0])), 
# 'lowerWall':	( 'fixedGradient' 		, np.array([0.0])), 
# 'frontAndBack':	( 'empty'			    , None			 )
# }

"""--------------------------------------------------------------------------"""

internalField_U = np.array([0.0, 0.0, 0.0])
internalField_p = np.array([0.0])

# Kinematic viscosity
nu = 0.01

# Time control
startTime = 0.0

writeTime = 0.0001

deltaT = 0.0001

endTime = 0.5

# Finite volume schemes

ddtScheme = "Euler"

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
'relTol'		:	0
}

fvSolution_p = {
'solver'		:	'GaussSeidel',
# 'solver'		:	'PointJacobi',
# 'minIter'		:	0,
# 'maxIter'		:	10000,
'tolerance'		:	1e-5,
'relTol'		:	0,
'refCell'		:	0,
'refValue'		:	np.array([0.0])
}

PISO = {
'nCorrectors'	:	3
}
