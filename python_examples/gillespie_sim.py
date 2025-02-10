import sys
sys.path.append('../build')

import SGEN_Py as sg

import pyvista as pv
import numpy as np

Dendrite_length = 200 #um
N_dendritic_segments = 15

soma = sg.Soma("soma", 20, x=0, y=0, z=0, radius=20)

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

neuron = sg.Neuron(soma, "Test_neuron")

ge = sg.Gillespie_engine(neuron) # Initialise Gillespie_engine for a given neuron 
# Run Gillespie algorithm for 100 hours recording every hour to test_Gillespie_out.csv.
# time_offset=-10 tells that the algorithm starts at t=-10 (useful for plotting sometimes)
ge.run_Gillespie(record_times=np.arange(0,100), file_name="test_Gillespie_out.csv", time_offset=-10)
