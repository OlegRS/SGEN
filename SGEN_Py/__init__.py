import sys

# Ensure we can find the compiled C++ module
sys.path.append('../build')

# Import the compiled C++ module
from _SGEN_Py import *  # Import everything from the compiled module

# Import Python wrapper classes
from .Neuron import Neuron  # Import Neuron class
