#!/bin/bash

# Quick incremental compilation - only compiles modified files
# Weria Pezeshkian
# Niels Bohr International Academy
# Niels Bohr Institute
# University of Copenhagen

cd dts_src

# Only compile the files we've modified
echo "Compiling modified files only..."
g++ -c -O3 -std=c++11 EvolveVerticesByMetropolisAlgorithm.cpp
g++ -c -O3 -std=c++11 CurvatureByShapeOperatorType1.cpp
g++ -c -O3 -std=c++11 MC_Simulation.cpp

# Link all object files together
echo "Linking..."
g++ -o DTS *.o
mv DTS ../
cd ..

echo "Quick compilation complete!"

