import sys
sys.path.append('../')

import SGEN_Py as sg
import pyvista as pv
import numpy as np

import matplotlib.pyplot as plt

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

# Computing expectations and correlations
expectations = neuron.expected_counts(dict_return=True)
correlations = neuron.correlations()

# sim_results = neuron.stationary_Gillespie_sim(record_times=range(0,3000,1), output_file_name='test_Gillespie')

sim_results = neuron.load_Gillespie_sim(file_name='test_Gillespie.csv')

# Define compartment to plot for
compartment = s_12_2

plt.plot(sim_results['time'], sim_results[compartment.name() + '_prot'], color='blue', label='Gillespie trajectory')
plt.axhline(expectations['prot'][compartment.name()], color='cyan', alpha=.7, linewidth=2, label='Stationary expectation')

# Plotting stationary standard deviation
std = neuron.standard_deviation(compartment)
plt.axhline(expectations['prot'][compartment.name()] + std, color='cyan', alpha=.7, linewidth=2, label='Stationary expectation' + r'$\pm\sigma$', linestyle='--')
plt.axhline(expectations['prot'][compartment.name()] - std, color='cyan', alpha=.7, linewidth=2, linestyle='--')

plt.xlabel('Time in hours')
plt.ylabel('Protein count')

plt.legend()
plt.show()
