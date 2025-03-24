import SGEN_Py as sg

import numpy as np
import pyvista as pv
import imageio

def create_neuron_gif(neuron, record_times, output_gif="neuron_simulation.gif"):
    # Run Gillespie simulation
    gillespie_data = neuron.Gillespie_sim(record_times)
    
    # Extract protein concentration data over time
    protein_data = {key: value for key, value in gillespie_data.items() if "_prot" in key}
    protein_concentrations = np.array(list(protein_data.values()))
    
    # PyVista setup
    plotter = pv.Plotter(off_screen=True)
    frames = []
    
    for t_index, time in enumerate(record_times):
        # Create visualization values from the simulation data at time t_index
        visualisation_values = protein_concentrations[:,t_index]/neuron.volumes()
        # print(f"Time {time}: {visualisation_values}")
        
        # Draw neuron with the current protein levels
        neuron.draw_3d(visualisation_values=visualisation_values, plotter=plotter)
        
        # Capture screenshot
        img = plotter.screenshot(return_img=True)
        frames.append(img)
        plotter.clear()
    
    # Save as GIF
    imageio.mimsave(output_gif, frames, fps=10)
    print(f"GIF saved as {output_gif}")

### Creating a neuron
Dendrite_length = 200 #um
N_dendritic_segments = 5

soma = sg.Soma("soma", 20, x=0, y=0, z=0, radius=20)

primary_branch = [sg.Dendritic_segment(parent=soma,
                                       name = "d_1-1",
                                       length = Dendrite_length/N_dendritic_segments)]

for i in np.arange(1,N_dendritic_segments):
    primary_branch.append(sg.Dendritic_segment(parent=primary_branch[i-1],
                                               name="d_1-" + str(i+1),
                                               length=Dendrite_length/N_dendritic_segments))

secondary_branch_1 = [sg.Dendritic_segment(parent=primary_branch[N_dendritic_segments-1],
                                           name="d_1_1-1",
                                           length=Dendrite_length/N_dendritic_segments,
                                           d_theta=30*np.pi/360,
                                           d_phi=0)]
for i in np.arange(1,N_dendritic_segments):
    secondary_branch_1.append(sg.Dendritic_segment(parent=secondary_branch_1[i-1],
                                                   name="d_1_1-" + str(i+1),
                                                   length=Dendrite_length/N_dendritic_segments))

secondary_branch_2 = [sg.Dendritic_segment(parent=primary_branch[N_dendritic_segments-1],
                                           name="d_1_2-1",
                                           length=Dendrite_length/N_dendritic_segments,
                                           d_theta=-30*np.pi/360,
                                           d_phi=0)]
for i in np.arange(1,N_dendritic_segments):
    secondary_branch_2.append(sg.Dendritic_segment(parent=secondary_branch_2[i-1],
                                                   name="d_1_2-" + str(i+1),
                                                   length=Dendrite_length/N_dendritic_segments,
                                                   radius=5*np.exp(-1/50*i)))

# Creating dendritic spines
s_1_1 = sg.Spine(parent=primary_branch[N_dendritic_segments//3],
                 name="s_1_1",
                 length=10,
                 radius=1)
s_1_2 = sg.Spine(parent=primary_branch[2*N_dendritic_segments//3],
                 name="s_1_2",
                 length=10,
                 radius=1,
                 d_theta=-np.pi/2)
s_11_1 = sg.Spine(parent=secondary_branch_1[N_dendritic_segments//3],
                  name = "s_11_1",
                  length=10,
                  radius=1)
s_11_2 = sg.Spine(parent=secondary_branch_1[2*N_dendritic_segments//3],
                  name = "s_11_2",
                  length=10,
                  radius=1,
                  d_theta=-np.pi/2)
s_12_1 = sg.Spine(parent=secondary_branch_2[N_dendritic_segments//3],
                  name="s_12_1",
                  length=10,
                  radius=1,
                  d_theta=-np.pi/2)
s_12_2 = sg.Spine(parent=secondary_branch_2[2*N_dendritic_segments//3],
                  name="s_12_2",
                  length=10,
                  radius=1)

neuron = sg.Neuron(soma, "Test_neuron")

# Example Usage
record_times = np.linspace(0, 1000, num=100)  # Define time points for simulation
create_neuron_gif(neuron, record_times)
