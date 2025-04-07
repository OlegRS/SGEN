import sys
import os

import SGEN_Py as sg

import numpy as np
import pyvista as pv
import imageio

import SGEN_Py as sg
import numpy as np
import pyvista as pv
import imageio

def create_neuron_gif(neuron, record_times, output_gif="neuron_sim_br_p2_with_mRNA.gif", camera_position=None, frames_dir=None):
    expectations = neuron.expected_counts()
    gillespie_data = neuron.stationary_Gillespie_sim(record_times)
    
    # Extract protein concentration data over time
    protein_data = {key: value for key, value in gillespie_data.items() if "_prot" in key}
    protein_counts = np.array(list(protein_data.values()))
    protein_concentrations = np.array([protein_counts[:, t_ind] / neuron.volumes() for t_ind in range(record_times.shape[0])]).T

    # Extract protein concentration data over time
    mRNA_data = {key: value for key, value in gillespie_data.items() if "_mRNA" in key}
    mRNA_counts = np.array(list(mRNA_data.values()))
    
    # Compute global min/max for consistent color bar
    global_min = np.min(np.log(protein_concentrations.flatten()))    
    print("global_min = ", global_min)
    global_max = np.max(np.log(protein_concentrations.flatten()))
    print("global_max = ", global_max)

    # PyVista setup
    frames = []
    
    for t_index, time in enumerate(record_times):
        plotter = pv.Plotter(off_screen=True)

        if camera_position is not None:
            plotter.camera_position = camera_position
        
        print(f"Time {time}:\n{protein_concentrations[:, t_index]}")  # Debugging output
        # Draw neuron with fixed color bar range
        mRNA_counts = {key: value[t_index] for key, value in mRNA_data.items()}
        active_gene_count = next((value for key, value in mRNA_data.items() if "gene" in key), None)
        neuron.draw_3d(visualisation_values=np.log(protein_concentrations[:, t_index]), plotter=plotter, clim=[global_min, global_max], mRNA_visualisation=True, mRNA_data=mRNA_counts, gene_visualisation=True, active_gene_count=active_gene_count)
        
        # Add time annotation
        plotter.add_text(f"Time: {time:.0f} hours", position='upper_left', font_size=12, color="black")

        # Capture screenshot
        img = plotter.screenshot(return_img=True, window_size=[1920, 1080])
        plotter.render()
        if frames_dir is not None:
            imageio.imwrite(f"{frames_dir}/frame_{t_index:04d}.png", img)
        frames.append(img)
        plotter.close()
        del plotter
    
    # Save as GIF
    imageio.mimsave(output_gif, frames, fps=10)
    print(f"GIF saved as {output_gif}")

### Creating a neuron
Dendrite_length = 200 #um
N_dendritic_segments = 15

soma = sg.Soma("soma", 20, x=0, y=0, z=0, radius=20, transcription_rate=.216)

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
                 radius=1,
                 binding_rate=.2)
s_1_2 = sg.Spine(parent=primary_branch[2*N_dendritic_segments//3],
                 name="s_1_2",
                 length=10,
                 radius=1,
                 d_theta=-np.pi/2,
                 binding_rate=.2)
s_11_1 = sg.Spine(parent=secondary_branch_1[N_dendritic_segments//3],
                  name = "s_11_1",
                  length=10,
                  radius=1,
                  binding_rate=.2)
s_11_2 = sg.Spine(parent=secondary_branch_1[2*N_dendritic_segments//3],
                  name = "s_11_2",
                  length=10,
                  radius=1,
                  d_theta=-np.pi/2,
                  binding_rate=.2)
s_12_1 = sg.Spine(parent=secondary_branch_2[N_dendritic_segments//3],
                  name="s_12_1",
                  length=10,
                  radius=1,
                  d_theta=-np.pi/2,
                  binding_rate=.2)
s_12_2 = sg.Spine(parent=secondary_branch_2[2*N_dendritic_segments//3],
                  name="s_12_2",
                  length=10,
                  radius=1,
                  binding_rate=.2)

neuron = sg.Neuron(soma, "Test_neuron")

# Example Usage
record_times = np.linspace(0, 100, num=500)  # Define time points for simulation
create_neuron_gif(neuron, record_times, camera_position=[(0, Dendrite_length*3.2, Dendrite_length), (0, 0, Dendrite_length), (1, 0, 0)])
