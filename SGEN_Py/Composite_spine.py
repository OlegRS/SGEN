import sys
sys.path.append('../build')

import _SGEN_Py as _sg
import numpy as np
import pyvista as pv

class Composite_spine:
    def __init__(self, parent, name='no_name', neck_length=10, neck_radius=1, neck_mRNA_decay_rate=0.0432, neck_translation_rate=75.6, neck_protein_decay_rate=0.004356, neck_mRNA_diffusion_constant=3.4e-3, neck_protein_diffusion_constant=.24, neck_mRNA_forward_trafficking_velocity=.5e-2, neck_mRNA_backward_trafficking_velocity=.1e-2, neck_protein_forward_trafficking_velocity=0, neck_protein_backward_trafficking_velocity=0, head_length=3, head_radius=5, head_mRNA_decay_rate=0.0432, head_translation_rate=75.6, head_protein_decay_rate=0.004356, head_mRNA_diffusion_constant=3.4e-3, head_protein_diffusion_constant=.24, head_mRNA_forward_trafficking_velocity=.5e-2, head_mRNA_backward_trafficking_velocity=.1e-2, head_protein_forward_trafficking_velocity=0, head_protein_backward_trafficking_velocity=0, PSD_length=.1, PSD_radius=4, PSD_binding_rate=.6, PSD_unbinding_rate=6, d_theta=np.pi/2, d_phi=0):
        """Dendritic spine consisting of a neck, head and post-synaptic density (PSD)"""
        self.neck = _sg.Dendritic_segment(parent, name+"_neck", neck_length, neck_radius, d_theta, d_phi, neck_mRNA_decay_rate, neck_translation_rate, neck_protein_decay_rate, neck_mRNA_diffusion_constant, neck_protein_diffusion_constant, neck_mRNA_forward_trafficking_velocity, neck_mRNA_backward_trafficking_velocity, neck_protein_forward_trafficking_velocity, neck_protein_backward_trafficking_velocity, middle_placement=True)
        self.head = _sg.Dendritic_segment(self.neck, name+"_head", head_length, head_radius, 0, 0, head_mRNA_decay_rate, head_translation_rate, head_protein_decay_rate, head_mRNA_diffusion_constant, head_protein_diffusion_constant, head_mRNA_forward_trafficking_velocity, head_mRNA_backward_trafficking_velocity, head_protein_forward_trafficking_velocity, head_protein_backward_trafficking_velocity)
        self.PSD = _sg.Spine(self.head, name+"_PSD", PSD_length, PSD_radius, PSD_binding_rate, PSD_unbinding_rate, 0, 0)
