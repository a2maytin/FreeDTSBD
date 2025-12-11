#!/bin/bash

# Quick incremental compilation with OpenMP - only compiles modified files
# Weria Pezeshkian
# Niels Bohr International Academy
# Niels Bohr Institute
# University of Copenhagen

cd dts_src

# Only compile the files we've modified
echo "Compiling modified files only (OpenMP)..."
g++ -c -O3 -fopenmp -std=c++11 EvolveVerticesByMetropolisAlgorithm.cpp
g++ -c -O3 -fopenmp -std=c++11 EvolveVerticesByMetropolisAlgorithmWithOpenMPType1.cpp
g++ -c -O3 -fopenmp -std=c++11 CurvatureByShapeOperatorType1.cpp
g++ -c -O3 -fopenmp -std=c++11 MC_Simulation.cpp
g++ -c -O3 -fopenmp -std=c++11 DNARepulsion.cpp
g++ -c -O3 -fopenmp -std=c++11 HarmonicBondsList.cpp
g++ -c -O3 -fopenmp -std=c++11 angle.cpp
g++ -c -O3 -fopenmp -std=c++11 VertexHarmonicBounds.cpp

# Link all object files together
echo "Linking..."
g++ -fopenmp -o DTS *.o
mv DTS ../
cd ..

echo "Quick compilation complete!"

