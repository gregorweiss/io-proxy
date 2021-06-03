#!/bin/bash

# Note: This script used to work when the case setup was in this folder.
# To use, it should be added to case folders and modified accordingly

python3 -m cProfile -s cumulative simpleFoam.py > profileLog
