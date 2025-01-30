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

me = sg.Morphologic_engine(neuron)

segments = me.segments()
volumes = me.volumes()

ae = sg.Analytic_engine(neuron)
print("Computing mRNA expectations...")
mRNA_expectations = np.array(ae.stationary_mRNA_expectations())
print("Computing protein expectations...")
prot_expectations = np.array(ae.stationary_protein_expectations())
print("Computing gene-mRNA correlations...")
gene_mRNA_covariances = np.array(ae.stationary_gene_mRNA_covariances())
print("Computing mRNA-mRNA correlations...")
mRNA_mRNA_covariances = np.array(ae.stationary_mRNA_mRNA_covariances())
print("Computing gene-protein correlations...")
gene_prot_covariances = np.array(ae.stationary_gene_protein_covariances())
print("Computing mRNA-protein correlations...")
mRNA_prot_covariances = np.array(ae.stationary_mRNA_protein_covariances())
print("Computing protein-protein correlations...")
prot_prot_covariances = np.array(ae.stationary_protein_protein_covariances())

prot_FFs = [prot_prot_covariances[i,i]-prot_expectations[i]**2 for i in range(len(prot_expectations))]/prot_expectations

# ae.stationary_expectations_and_correlations()

# Extract the neuron segments and nodes
start_points = [segments[i][0][:3] for i in range(len(segments))]
end_points = [segments[i][1][:3] for i in range(len(segments))]
radii = [segments[i][0][3] for i in range(len(segments))]

# prot_expectations = np.genfromtxt("protein_expectations", delimiter='\n')
# prot_expectations = np.genfromtxt("protein_expectations.dat", delimiter='\n')
prot_concentrations = prot_expectations/volumes 
# prot_concentrations = [prot_expectations[i]/volumes[i] for i in range(len(volumes))]
segment_values = np.flip(prot_FFs)
# segment_values = np.log(prot_expectations)

# Convert the neuron morphology into a mesh for visualization in PyVista

# You can create a line or tube representation of the morphology by using segments (radii)
# Make a tube for each segment to represent dendrites
tubes = []
all_scalars = [] # Collect scalars for the final mesh
for i in range(0,len(segments)):
    start, end = start_points[i], end_points[i]
    radius = radii[i]  # Radius of the current segment
    # print("radius=", radius)
    # print("segment_values[i]=", segment_values[i])

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
