import sys
sys.path.append('../build')

import SGEN_Py as sg

import pyvista as pv
import numpy as np

def PCC(comp1, comp2, analytic_engine):
    return (analytic_engine.protein_protein_correlation(comp1,comp2) - analytic_engine.protein_expectation(comp1)*analytic_engine.protein_expectation(comp2))/(np.sqrt(analytic_engine.protein_protein_correlation(comp1,comp1)-analytic_engine.protein_expectation(comp1)**2)*np.sqrt(analytic_engine.protein_protein_correlation(comp2,comp2)-analytic_engine.protein_expectation(comp2)**2))
    
Dendrite_length = 200 #um
N_dendritic_segments = 100
N_soma_segments = N_dendritic_segments/4

soma_length = 20
soma = sg.Soma("soma",
               length=soma_length/N_dendritic_segments,
               radius=soma_length,
               translation_rate=75.6*200/soma_length*N_dendritic_segments,
               protein_diffusion_constant=0.24/(np.pi*100))

soma_transition_front = [sg.Dendritic_segment(protein_diffusion_constant=0.24/(np.pi*100),
                                        translation_rate=0,
                                        parent=soma,
                                        name = "soma_transition_front_0",
                                        length = soma_length/N_soma_segments,
                                        radius = soma_length)]
for i in np.arange(1,N_soma_segments):
    soma_transition_front.append(sg.Dendritic_segment(parent=soma_transition_front[-1],
                                                protein_diffusion_constant=0.24/(np.pi*100),
                                                translation_rate=0,
                                                name="soma_transition_front_" + str(i+1),
                                                length=soma_length/N_soma_segments,
                                                radius=np.sqrt(soma_length**2-(i*soma_length/N_soma_segments)**2)))

soma_transition_back = [sg.Dendritic_segment(protein_diffusion_constant=0.24/(np.pi*100),
                                             translation_rate=0,
                                             parent=soma,
                                             name = "soma_transition_back_0",
                                             length = soma_length/N_soma_segments,
                                             radius = soma_length,
                                             d_theta=np.pi)]
for i in np.arange(1,N_soma_segments):
    soma_transition_back.append(sg.Dendritic_segment(parent=soma_transition_back[-1],
                                                     protein_diffusion_constant=0.24/(np.pi*100),
                                                     translation_rate=0,
                                                     name="soma_transition_back_" + str(i+1),
                                                     length=soma_length/N_soma_segments,
                                                     radius=np.sqrt(soma_length**2-(i*soma_length/N_soma_segments)**2)))

primary_branch = [sg.Dendritic_segment(protein_diffusion_constant=0.24/(np.pi*100),
                                       translation_rate=0,
                                       parent=soma_transition_front[-1],
                                       name = "d_1",
                                       length = Dendrite_length/N_dendritic_segments)]
for i in np.arange(1,N_dendritic_segments):
    primary_branch.append(sg.Dendritic_segment(protein_diffusion_constant=0.24/(np.pi*100),
                                               translation_rate=0,
                                               parent=primary_branch[-1],
                                               name="d_" + str(i+1),
                                               length=Dendrite_length/N_dendritic_segments))
  
secondary_branch_1 = [sg.Dendritic_segment(protein_diffusion_constant=0.24/(np.pi*100),
                                           translation_rate=0,
                                           parent=primary_branch[-1],
                                           name="d_1_1",
                                           length=Dendrite_length/N_dendritic_segments,
                                           d_theta=30*np.pi/360,
                                           d_phi=0)]
for i in np.arange(1,N_dendritic_segments):
    secondary_branch_1.append(sg.Dendritic_segment(protein_diffusion_constant=0.24/(np.pi*100),
                                                   translation_rate=0,
                                                   parent=secondary_branch_1[-1],
                                                   name="d_1_" + str(i+1),
                                                   length=Dendrite_length/N_dendritic_segments))

secondary_branch_2 = [sg.Dendritic_segment(protein_diffusion_constant=0.24/(np.pi*100),
                                           translation_rate=0,
                                           parent=primary_branch[-1],
                                           name="d_1_2",
                                           length=Dendrite_length/N_dendritic_segments,
                                           d_theta=-30*np.pi/360,
                                           d_phi=0)]
for i in np.arange(1,N_dendritic_segments):
    secondary_branch_2.append(sg.Dendritic_segment(protein_diffusion_constant=0.24/(np.pi*100),
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

neuron = sg._Neuron(soma, "Test_neuron")

ae = sg.Analytic_engine(neuron)
print("Computing mRNA expectations...")
mRNA_expectations = np.array(ae.stationary_mRNA_expectations())
print("Computing protein expectations...")
prot_expectations = np.array(ae.stationary_protein_expectations())

# print("Computing gene-mRNA correlations...")
# gene_mRNA_covariances = np.array(ae.stationary_gene_mRNA_covariances())
# print("Computing mRNA-mRNA correlations...")
# mRNA_mRNA_covariances = np.array(ae.stationary_mRNA_mRNA_covariances())
# print("Computing gene-protein correlations...")
# gene_prot_covariances = np.array(ae.stationary_gene_protein_covariances())
# print("Computing mRNA-protein correlations...")
# mRNA_prot_covariances = np.array(ae.stationary_mRNA_protein_covariances())
# print("Computing protein-protein correlations...")
# prot_prot_covariances = np.array(ae.stationary_protein_protein_covariances())

# prot_FFs = [prot_prot_covariances[i,i]-prot_expectations[i]**2 for i in range(len(prot_expectations))]/prot_expectations

# Get the neuron segments and their volumes
me = sg.Morphologic_engine(neuron)
segments, volumes = me.segments(), me.volumes()

start_points = [segments[i][0][:3] for i in range(len(segments))]
end_points = [segments[i][1][:3] for i in range(len(segments))]
radii = [segments[i][1][3] for i in range(len(segments))]

prot_concentrations = prot_expectations/volumes
mRNA_concentrations = mRNA_expectations/volumes
# segment_values = np.flip(mRNA_expectations)
segment_values = np.flip(prot_expectations)
# segment_values = np.flip(np.log(prot_expectations))

# Convert the neuron morphology into a mesh for visualization in PyVista
# Make a tube for each segment
tubes = []
all_scalars = [] # Collect scalars for the final mesh
for i in range(0,len(segments)):
    start, end = start_points[i], end_points[i]
    radius = radii[i]  # Radius of the current segment

    # Create a tube (cylinder) between two coordinates
    tube = pv.Line(start, end).tube(radius=radius)
    tubes.append(tube)

    # Extend scalars to match the number of points in the current tube
    all_scalars.append(np.full(tube.n_cells, segment_values[i]))

# Combine all tubes into a single mesh
neuron_mesh = tubes[0] if tubes else None
for tube in tubes[1:]:
    neuron_mesh += tube

# Concatenate scalars for all segments
flat_scalars = np.concatenate(all_scalars)

# Assign scalars to the combined mesh
neuron_mesh.cell_data["Protein Levels"] = flat_scalars

# Visualize the neuron
plotter = pv.Plotter()
plotter.add_mesh(neuron_mesh, scalars="Protein Levels", cmap="coolwarm", show_edges=False)
plotter.show_axes()
plotter.show()

# Pearson correlation coefficient
# (ae.protein_protein_correlation(s_1_1, s_12_1) - ae.protein_expectation(s_1_1)*ae.protein_expectation(s_12_1))/(np.sqrt(ae.protein_protein_correlation(s_1_1, s_1_1)-ae.protein_expectation(s_1_1)**2)*np.sqrt(ae.protein_protein_correlation(s_12_1, s_12_1)-ae.protein_expectation(s_12_1)**2))
