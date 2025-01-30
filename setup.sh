#!/bin/bash

rm -rf build

mkdir build
cd build
echo "Running from $(pwd)"

cmake ..

# cmake --build . --target optimised_moments
cmake --build . --target SGEN_Py

# sudo make install

export PYTHONPATH=$(pwd):$PYTHONPATH
