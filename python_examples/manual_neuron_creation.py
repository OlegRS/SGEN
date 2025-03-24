import sys
sys.path.append('../')

import sys
import os

# Add the parent directory (../) to the start of sys.path
sys.path.insert(0, os.path.abspath('../'))

import SGEN_Py as sg
import numpy as np

Dendrite_length = 200 #um
N_dendritic_segments = 5

soma = sg.Soma("soma", length=20, x=0, y=0, z=0, radius=20)

primary_branch = [sg.Dendritic_segment(parent=soma,
                                       name = "d_1-1",
                                       length = Dendrite_length/N_dendritic_segments)]

for i in np.arange(1,N_dendritic_segments):
    primary_branch.append(sg.Dendritic_segment(parent=primary_branch[i-1],
                                               name="d_1-" + str(i+1),
                                               length=Dendrite_length/N_dendritic_segments))

secondary_branch_1 = [sg.Dendritic_segment(parent=primary_branch[N_dendritic_segments-1],
                                           name="d_1_1-0",
                                           length=Dendrite_length/N_dendritic_segments,
                                           d_theta=30*np.pi/360,
                                           d_phi=0)]
for i in np.arange(1,N_dendritic_segments):
    secondary_branch_1.append(sg.Dendritic_segment(parent=secondary_branch_1[i-1],
                                                   name="d_1_1-" + str(i+1),
                                                   length=Dendrite_length/N_dendritic_segments))

secondary_branch_2 = [sg.Dendritic_segment(parent=primary_branch[N_dendritic_segments-1],
                                           name="d_1_1-0",
                                           length=Dendrite_length/N_dendritic_segments,
                                           d_theta=-30*np.pi/360,
                                           d_phi=0)]
for i in np.arange(1,N_dendritic_segments):
    secondary_branch_2.append(sg.Dendritic_segment(parent=secondary_branch_2[i-1],
                                                   name="d_1_1-" + str(i+1),
                                                   length=Dendrite_length/N_dendritic_segments,
                                                   radius=5*np.exp(-1/50*i)))


# Creating dendritic spines
s_1_1 = sg.Spine(parent=primary_branch[int(N_dendritic_segments/3)],
                   name="s_1_1",
                   length=10,
                   radius=1)
s_1_2 = sg.Spine(parent=primary_branch[int(2*N_dendritic_segments/3)],
                   name="s_1_2",
                   length=10,
                   radius=1)
s_11_1 = sg.Spine(parent=secondary_branch_1[int(N_dendritic_segments/3)],
                    name = "s_11_1",
                    length=10,
                    radius=1)
s_11_2 = sg.Spine(parent=secondary_branch_1[int(2*N_dendritic_segments/3)],
                    name = "s_11_2",
                    length=10,
                    radius=1)
s_12_1 = sg.Spine(parent=secondary_branch_2[int(N_dendritic_segments/3)],
                  name="s_12_1",
                  length=10,
                  radius=1,
                  d_theta=-np.pi/2)
s_12_2 = sg.Spine(parent=secondary_branch_2[int(2*N_dendritic_segments/3)],
                  name="s_12_2",
                  length=10,
                  radius=1,
                  d_theta=-np.pi/2)

neuron = sg.Neuron(soma, "neuron")

expectations = neuron.expected_counts()

correlation = neuron.correlations()

sim_results = neuron.stationary_Gillespie_sim(range(0,100,5))

# Creating a new neuron specified in .swc file
# neuron = sg.Neuron("../data/morphologies/10_2REDO-850-GM18-Ctl-Ctl-Chow-BNL16A-CA1_Finished2e.swc")
swc_neuron = sg.Neuron("../data/morphologies/DD13-10-c6-2.CNG.swc")
# Computing expectations
expectations = swc_neuron.expected_counts()
# Visualising
swc_neuron.draw_3d()

print('Printing soma...')
print(swc_neuron.soma())
