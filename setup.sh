#!/bin/bash

# Install required dependencies
apt install -y cmake
pip install pybind11 numpy pyvista -qq

# Find where Pybind11 is installed
python3 -m pybind11 --cmakedir

pip install .
