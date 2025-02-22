import sys
sys.path.append('../')

import SGEN_Py as sg
import pyvista as pv
import numpy as np

# Load the neuron morphology from the SWC file
# neuron = sg.Neuron("../data/morphologies/DD13-10-c6-2.CNG.swc")

# neuron = sg.Neuron("../data/morphologies/DD13-10-c5-1.CNG_.swc")

neuron = sg.Neuron("../data/morphologies/10_2REDO-850-GM18-Ctl-Ctl-Chow-BNL16A-CA1_Finished2e.swc")

expectations = neuron.expected_counts()

neuron.draw_3d(expectations['prot']/neuron.volumes())
