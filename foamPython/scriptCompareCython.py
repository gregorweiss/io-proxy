import numpy as np
import timeit
import time as timeModule
import os

solver = "simpleFoam.py"

noRuns = 100

resultsFile = "gignore/cythonComparisonResults.dat"


# Calculate average times without Cython

os.system("sh scriptUnsetCython.sh > /dev/null 2>&1")

noCyTimes = np.empty(0, dtype=float)

for i in range(noRuns):
	
	os.system("rm -r [1-9]* > /dev/null 2>&1 ; rm -r 0.[0-9]* > /dev/null 2>&1")
	
	start = timeModule.perf_counter()
	# os.system("python3 " + solver)
	os.system("python3 " + solver + "> /dev/null 2>&1")
	end = timeModule.perf_counter()
	
	noCyTimes = np.append(noCyTimes, (end - start))
	
	print("No Cython run", i+1, "done.")



# Calculate average times with Cython

os.system("sh scriptSetCython.sh > /dev/null 2>&1")

CyTimes = np.empty(0, dtype=float)

for i in range(noRuns):
	
	os.system("rm -r [1-9]* > /dev/null 2>&1 ; rm -r 0.[0-9]* > /dev/null 2>&1")
	
	start = timeModule.perf_counter()
	os.system("python3 " + solver + "> /dev/null 2>&1")
	end = timeModule.perf_counter()
	
	CyTimes = np.append(CyTimes, (end - start))
	
	print("Cython run", i+1, "done.")

print("Cython speedup =", np.average(noCyTimes) / np.average(CyTimes))



convergenceArray = np.empty(noRuns, dtype = float)

for i in range(noRuns):
	convergenceArray[i] = np.average(noCyTimes[:(i+1)]) / np.average(CyTimes[:(i+1)])

print("Convergence history is being written to " + resultsFile)

file = open(resultsFile, 'w')

# file.write("#noRuns    speed-up\n")

for i in range(noRuns):
	file.write(str(i+1) + " " + str(convergenceArray[i]) + "\n")
  
