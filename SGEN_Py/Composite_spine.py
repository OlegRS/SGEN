import sys
sys.path.append('../build')

import _SGEN_Py as _sg
import numpy as np
import pyvista as pv

class Composite_spine:
    def __init__(self, parent, name='no_name', neck_length=10, neck_radius=1, head_length=3, head_radius=5, PSD_length=.1, PSD_radius=4, PSD_binding_rate=.6, PSD_unbinding_rate=6, d_theta=np.pi/2, d_phi=0):
        """Dendritic spine consisting of a neck, head and post-synaptic density (PSD)"""
        self.neck = _sg.Dendritic_segment(parent, name+"_neck", neck_length, neck_radius, d_theta, d_phi)
        self.head = _sg.Dendritic_segment(self.neck, name+"_head", head_length, head_radius)
        self.PSD = _sg.Spine(self.head, name+"_PSD", PSD_length, PSD_radius, PSD_binding_rate, PSD_unbinding_rate, 0, 0)
