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

def create_neuron_gif(neuron, record_times, output_gif="neuron_sim_br_p2_with_mRNA.gif", camera_position=None, frames_dir=None, sim_file_name=None, cbar_scale="linear"):
    expectations = neuron.expected_counts()['prot']
    if sim_file_name is None:
        gillespie_data = neuron.stationary_Gillespie_sim(record_times, output_file_name="Gillespie_ordinary_spines", seed=1)
    else:
        gillespie_data = neuron.load_Gillespie_sim(sim_file_name)
        record_times = gillespie_data["time"]
    
    # Extract protein concentration data over time
    protein_data = {key: value for key, value in gillespie_data.items() if "_prot" in key}
    protein_counts = np.array(list(protein_data.values()))
    protein_concentrations = np.array([protein_counts[:, t_ind] / neuron.volumes() for t_ind in range(record_times.shape[0])]).T

    # Compute global min/max for consistent color bar
    global_min = np.min(protein_concentrations.flatten())
    global_max = np.max(protein_concentrations.flatten())
    print("global_min = ", global_min)
    print("global_max = ", global_max)

    if cbar_scale == "uniform":
        expected_concentrations = expectations / neuron.volumes()
        expected_concentrations_appended = np.append(expected_concentrations, [global_min - 1e-7, global_max + 1e-7])
        sorted_ec = np.sort(expected_concentrations_appended)
        rescaled_pc = np.zeros(protein_concentrations.shape)
        for t_ind in range(protein_concentrations.shape[1]):
            for i in range(protein_concentrations.shape[0]):
                m = np.searchsorted(sorted_ec, protein_concentrations[i,t_ind])-1 # Find index: sorted_ec[m] < protein_concentrations[i,t_ind] < sorted_ec[m+1]
                rescaled_pc[i, t_ind] = m + (protein_concentrations[i,t_ind] - sorted_ec[m])/(sorted_ec[m+1] - sorted_ec[m])
        ticks = []
        tick_labels = []
        for i, ec in enumerate(sorted_ec):
            if i == len(sorted_ec) - 1:
                continue
            mid_rescaled = i + .5
            ticks.append(mid_rescaled)
            tick_labels.append(f"{(sorted_ec[i] + sorted_ec[i+1])/2:.2f}")
            
        scalar_bar_args = {
            "title": "Protein Concentration",
            "n_labels": len(ticks),
            "vertical": False,
            "fmt": "",  # disables auto-format
            "color": "black",
            "title_font_size": 18,
            "label_font_size": 12,
            "labels": True,
            "custom_labels": dict(zip(ticks, tick_labels))  # CRUCIAL
        }
    elif cbar_scale == "log":
        rescaled_pc = np.log(protein_concentrations)
    elif cbar_scale == "linear":
        rescaled_pc = protein_concentrations
    else:
        raise ValueError("cbar_scale can be either \"uniform\", \"log\" or \"linear\"")

    rescaled_global_min = np.min(rescaled_pc.flatten())
    rescaled_global_max = np.max(rescaled_pc.flatten())
    print("rescaled_global_min = ", rescaled_global_min)
    print("rescaled_global_max = ", rescaled_global_max)


    # Extract protein concentration data over time
    mRNA_data = {key: value for key, value in gillespie_data.items() if "_mRNA" in key}
    active_gene_counts = next((value for key, value in gillespie_data.items() if "gene" in key), None)
    
    # r_max = 5
    # spines = neuron.spines()

    # PyVista setup
    frames = []
    
    for t_index, time in enumerate(record_times):
        # for spine in spines:
        #     prot_count = next((value for key, value in gillespie_data.items() if spine.name() in key), None)
        #     print("r_max*spine.protein_count()/global_max/(np.pi*.1*4**2)", r_max*spine.protein_count()/global_max/(np.pi*.1*4**2))
        #     spine.set_radius(r_max*prot_count/global_max/(np.pi*.1*4**2))

        plotter = pv.Plotter(off_screen=True)

        if camera_position is not None:
            plotter.camera_position = camera_position
        
        # Draw neuron with fixed color bar range
        mRNA_counts = {key: value[t_index] for key, value in mRNA_data.items()}
        active_gene_count = active_gene_counts[t_index]
        print(f"Frame: {t_index}; Time: {time}")  # Debugging output
        
        neuron.draw_3d(visualisation_values=rescaled_pc[:, t_index], plotter=plotter, scalar_bar_args=scalar_bar_args, clim=[rescaled_global_min, rescaled_global_max], mRNA_visualisation=True, mRNA_data=mRNA_counts, gene_visualisation=True, active_gene_count=active_gene_count, mRNA_radius=1, opacity=.5)
        
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

soma = sg.Soma("soma", length=20, x=0, y=0, z=0, radius=20, transcription_rate=.216, n_gene_copies=2, translation_rate=75.6/10, protein_forward_trafficking_velocity=1.6e-05*30)

primary_branch = [sg.Dendritic_segment(parent=soma,
                                       name = "d_1-1",
                                       length = Dendrite_length/N_dendritic_segments, translation_rate=75.6/10)]

for i in np.arange(1,N_dendritic_segments):
    primary_branch.append(sg.Dendritic_segment(parent=primary_branch[i-1],
                                               name="d_1-" + str(i+1),
                                               length=Dendrite_length/N_dendritic_segments, translation_rate=75.6/10))

secondary_branch_1 = [sg.Dendritic_segment(parent=primary_branch[N_dendritic_segments-1],
                                           name="d_1_1-1",
                                           length=Dendrite_length/N_dendritic_segments,
                                           d_theta=30*np.pi/360,
                                           d_phi=0, translation_rate=75.6/10)]
for i in np.arange(1,N_dendritic_segments):
    secondary_branch_1.append(sg.Dendritic_segment(parent=secondary_branch_1[i-1],
                                                   name="d_1_1-" + str(i+1),
                                                   length=Dendrite_length/N_dendritic_segments, translation_rate=75.6/10))

secondary_branch_2 = [sg.Dendritic_segment(parent=primary_branch[N_dendritic_segments-1],
                                           name="d_1_2-1",
                                           length=Dendrite_length/N_dendritic_segments,
                                           d_theta=-30*np.pi/360,
                                           d_phi=0, translation_rate=75.6/10)]
for i in np.arange(1,N_dendritic_segments):
    secondary_branch_2.append(sg.Dendritic_segment(parent=secondary_branch_2[i-1],
                                                   name="d_1_2-" + str(i+1),
                                                   length=Dendrite_length/N_dendritic_segments, translation_rate=75.6/10))

# Creating dendritic spines
# cs_1_1 = sg.Composite_spine(parent=primary_branch[N_dendritic_segments//3],
#                             name="s_1_1",
#                             d_phi=np.pi/4,
#                             PSD_binding_rate=.2)
# cs_1_2 = sg.Composite_spine(parent=primary_branch[2*N_dendritic_segments//3],
#                             name="s_1_2",
#                             d_phi=np.pi/4,
#                             PSD_binding_rate=.2)

s_1_1 = sg.Spine(parent=primary_branch[N_dendritic_segments//3],
                 name="s_1_1",
                 length=10,
                 radius=1,
                 placement="middle")
s_1_2 = sg.Spine(parent=primary_branch[2*N_dendritic_segments//3],
                 name="s_1_2",
                 length=10,
                 radius=1,
                 d_theta=-np.pi/2,
                 placement="middle")
s_11_1 = sg.Spine(parent=secondary_branch_1[N_dendritic_segments//3],
                  name = "s_11_1",
                  length=10,
                  radius=1,
                  binding_rate=.2,
                  placement="middle")
s_11_2 = sg.Spine(parent=secondary_branch_1[2*N_dendritic_segments//3],
                  name = "s_11_2",
                  length=10,
                  radius=1,
                  d_theta=-np.pi/2,
                  binding_rate=.2,
                  placement="middle")
s_12_1 = sg.Spine(parent=secondary_branch_2[N_dendritic_segments//3],
                  name="s_12_1",
                  length=10,
                  radius=1,
                  d_theta=-np.pi/2,
                  binding_rate=.2,
                  placement="middle")
s_12_2 = sg.Spine(parent=secondary_branch_2[2*N_dendritic_segments//3],
                  name="s_12_2",
                  length=10,
                  radius=1,
                  binding_rate=.2,
                  placement="middle")

neuron = sg.Neuron(soma, "neuron")

# Example Usage
record_times = np.linspace(0, 1000, num=500)  # Define time points for simulation
# expectations = neuron.expected_counts()
# create_neuron_gif(neuron, record_times, camera_position=[(0, Dendrite_length*3.2, Dendrite_length), (0, 0, Dendrite_length), (1, 0, 0)], sim_file_name="Gillespie_stat.csv", cbar_scale="uniform")
create_neuron_gif(neuron, record_times, camera_position=[(0, Dendrite_length*3.2, Dendrite_length), (0, 0, Dendrite_length), (1, 0, 0)], cbar_scale="uniform")
