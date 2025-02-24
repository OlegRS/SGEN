import sys
sys.path.append('../build')

import SGEN_Py as sg

import numpy as np

Dendrite_length = 200 #um
N_dendritic_segments = 10

soma_length = 20
soma = sg.Soma("soma",
               length=soma_length,
               radius=10,
               translation_rate=75.6*200/soma_length,
               protein_diffusion_constant=0.24)

primary_branch = [sg.Dendritic_segment(protein_diffusion_constant=0.24,
                                       translation_rate=0,
                                       parent=soma,
                                       name = "d_1",
                                       length = Dendrite_length/N_dendritic_segments)]
for i in np.arange(1,N_dendritic_segments):
    primary_branch.append(sg.Dendritic_segment(protein_diffusion_constant=0.24,
                                               translation_rate=0,
                                               parent=primary_branch[-1],
                                               name="d_" + str(i+1),
                                               length=Dendrite_length/N_dendritic_segments))
  
secondary_branch_1 = [sg.Dendritic_segment(protein_diffusion_constant=0.24,
                                           translation_rate=0,
                                           parent=primary_branch[-1],
                                           name="d_1_1",
                                           length=Dendrite_length/N_dendritic_segments,
                                           d_theta=30*np.pi/360,
                                           d_phi=0)]
for i in np.arange(1,N_dendritic_segments):
    secondary_branch_1.append(sg.Dendritic_segment(protein_diffusion_constant=0.24,
                                                   translation_rate=0,
                                                   parent=secondary_branch_1[-1],
                                                   name="d_1_" + str(i+1),
                                                   length=Dendrite_length/N_dendritic_segments))

secondary_branch_2 = [sg.Dendritic_segment(protein_diffusion_constant=0.24,
                                           translation_rate=0,
                                           parent=primary_branch[-1],
                                           name="d_1_2",
                                           length=Dendrite_length/N_dendritic_segments,
                                           d_theta=-30*np.pi/360,
                                           d_phi=0)]
for i in np.arange(1,N_dendritic_segments):
    secondary_branch_2.append(sg.Dendritic_segment(protein_diffusion_constant=0.24,
                                                   translation_rate=0,
                                                   parent=secondary_branch_2[-1],
                                                   name="d_1_1-" + str(i+1),
                                                   length=Dendrite_length/N_dendritic_segments))
                                                   # radius=5*np.exp(-1/200*i)))

# Placing "ribosomes"
# primary_branch[int(N_dendritic_segments/3)].set_translation_rate(75.6*200/(Dendrite_length/N_dendritic_segments))
secondary_branch_1[int(N_dendritic_segments/3)].set_translation_rate(75.6*200/(Dendrite_length/N_dendritic_segments))

# Creating dendritic spines
s_1_1 = sg.Spine(parent=primary_branch[int(N_dendritic_segments/3)],
                 name="s_1_1",
                 length=10,
                 radius=1)
s_1_2 = sg.Spine(parent=primary_branch[int(2*N_dendritic_segments/3)],
                 name="s_1_2",
                 length=10,
                 radius=1,
                 d_theta=-np.pi/2)
s_11_1 = sg.Spine(parent=secondary_branch_1[int(N_dendritic_segments/3)],
                  name = "s_11_1",
                  length=10,
                  radius=1)
s_11_2 = sg.Spine(parent=secondary_branch_1[int(2*N_dendritic_segments/3)],
                  name = "s_11_2",
                  length=10,
                  radius=1,
                  d_theta=-np.pi/2)
s_12_1 = sg.Spine(parent=secondary_branch_2[int(N_dendritic_segments/3)],
                  name="s_12_1",
                  length=10,
                  radius=1,
                  d_theta=-np.pi/2)
s_12_2 = sg.Spine(parent=secondary_branch_2[int(2*N_dendritic_segments/3)],
                  name="s_12_2",
                  length=10,
                  radius=1)

# Initialising neuron
neuron = sg.Neuron(soma, "Test_neuron")

expectations = neuron.expected_counts()

neuron.draw_3d(expectations['prot']/neuron.volumes())
