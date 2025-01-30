import sys
sys.path.append('../build')

import SGEN_Py as sg

import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

#####################################################################################

Dendrite_length = 5000 #um
N_dendritic_segments = 500

soma = sg.Soma(name="soma", length=Dendrite_length/N_dendritic_segments)

dendritic_segments = [sg.Dendritic_segment(parent=soma, name="d_1-1", length=Dendrite_length/N_dendritic_segments)]

for i in np.arange(1,N_dendritic_segments):
    dendritic_segments.append(sg.Dendritic_segment(parent=dendritic_segments[i-1], name="d_1-" + str(i+1), length=Dendrite_length/N_dendritic_segments))
    
neuron = sg.Neuron(soma)

ae = sg.Analytic_engine(neuron)

# print("Computing mRNA expectations...")
# mRNA_expectations = np.array(ae.stationary_mRNA_expectations())
# print("Computing protein expectations...")
# prot_expectations = np.array(ae.stationary_protein_expectations())
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

print("---------------------------------------------")
print("Computing expectations and correlations...")
# As for now, we need another analytic engine to compute again
# sg.Analytic_engine(neuron.refresh()).stationary_expectations_and_correlations()
# sg.Analytic_engine(neuron).stationary_expectations_and_correlations()
# sg.Analytic_engine(neuron).stationary_expectations_and_correlations(dt_mRNA=0.1/10, dt_prot=0.1, t_fin_mRNA=2314.81*10, t_fin_prot=2314.81)
# sg.Analytic_engine(neuron).stationary_expectations_and_correlations(dt_mRNA=.5, dt_prot=.5, t_fin_mRNA=700, t_fin_prot=700)
# ae.stationary_expectations_and_correlations()

##############################################################################################
##############################################################################################

gene_activation_rate = 1/12/3600
gene_deactivation_rate = 1/12/3600

tau_1 = 1/(gene_activation_rate+gene_deactivation_rate)

lambda_2 = transcription_rate = (3.*200/10000) * .001#*3600

mRNA_decay_rate = 1.2e-5#*3600
tau_2 = 1/mRNA_decay_rate

kappa = translation_rate = 0.021#*3600

protein_decay_rate = 1.21e-6#*3600
tau_3 = 1/protein_decay_rate

D_m = mRNA_diffusion_constant = 3.4e-3
D_p = protein_diffusion_constant = .24

mRNA_forward_trafficking_velocity = .5e-2
mRNA_backward_trafficking_velocity = .1e-2
v_m = mRNA_forward_trafficking_velocity - mRNA_backward_trafficking_velocity

protein_forward_trafficking_velocity = 0 # .5e-2/10
protein_backward_trafficking_velocity = 0
v_p = protein_forward_trafficking_velocity - protein_backward_trafficking_velocity

x_lim = Dendrite_length

lambda_m = (v_m - np.sqrt(v_m**2+4*D_m/tau_2))/(2*D_m)

n_exp = .5

fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(10*1.6*1.1, 3.2*1.9*1.7*4/3*4/5))

### mRNA concentrations
xi = np.linspace(0,x_lim, 100000)

R = (2*lambda_2*n_exp)/(v_m + np.sqrt(v_m**2+4*D_m/tau_2)) * np.exp(lambda_m * xi)
R_discrete = np.genfromtxt("mRNA_expectations_500.csv", delimiter=',')

axs[0,0].plot(xi, R, label="R_analytic")

axs[0,0].plot(np.arange(0,len(R_discrete)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), R_discrete * N_dendritic_segments/Dendrite_length, label="R_discrete")

# axs[0,0].plot(np.arange(0,len(mRNA_expectations)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), mRNA_expectations * N_dendritic_segments/Dendrite_length, label="mRNA_concentrations")

### Protein concentrations
Lambda_m = (v_p - np.sqrt(v_p**2+4*D_p/tau_3))/(2*D_p)
Lambda_p = (v_p + np.sqrt(v_p**2+4*D_p/tau_3))/(2*D_p)

P = (2*lambda_2*n_exp)/(v_m + np.sqrt(v_m**2+4*D_m/tau_2))*kappa/((lambda_m-Lambda_m)*(Lambda_p-lambda_m)*D_p) * (np.exp(lambda_m * xi) + 1/Lambda_p*(lambda_m - v_p/D_p)*np.exp(Lambda_m * xi))

P_discrete = np.genfromtxt("protein_expectations_500.csv", delimiter=',')

axs[1,0].plot(xi, P, label="P_analytic")
axs[1,0].plot(np.arange(0,len(P_discrete)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), P_discrete *N_dendritic_segments/Dendrite_length, label="P_discrete")

# axs[1,0].plot(np.arange(0,len(prot_expectations)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), prot_expectations *N_dendritic_segments/Dendrite_length, label="Prot_concentration")

for ax in axs[:2,0]:
    ax.set_xlim([0,x_lim])
    ax.axhline(0, color='black')
    ax.set_yscale('log')
    ax.legend()

axs[2,0].set_xlabel("Distance from the soma, " + r'$\mu m$', fontsize=17)
axs[0,0].set_ylabel("mRNA concentration, " + r'$\mu m^{-3}$', fontsize=17)
axs[1,0].set_ylabel("Protein concentration, " + r'$\mu m^{-3}$', fontsize=17)
axs[2,0].set_ylabel("Distance from the soma, " + r'$\mu m$', fontsize=17)

### gene-mRNA correlations
lambda_m_tilde = (v_m - np.sqrt(v_m**2+4*D_m*(1/tau_1+1/tau_2)))/(2*D_m)

G2_nm = 2*lambda_2*n_exp*gene_deactivation_rate*tau_1/(v_m + np.sqrt(v_m**2+4*D_m*(1/tau_1+1/tau_2)))*np.exp(lambda_m_tilde*xi) + n_exp*R

G2_nm_discrete = np.genfromtxt("gene_mRNA_covariances_500.csv", delimiter=',')

axs[0,1].plot(xi, G2_nm, label="G2_nm")
axs[0,1].plot(np.arange(0,len(G2_nm_discrete)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), G2_nm_discrete*N_dendritic_segments/Dendrite_length, label="G2_nm_discrete")

# axs[0,1].plot(np.arange(0,len(gene_mRNA_covariances)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), gene_mRNA_covariances *N_dendritic_segments/Dendrite_length, label="G2_nm_discrete_")

axs[0,1].set_ylabel(r"$G^{(2)}_{nm}$ density, " + r'$\mu m^{-3}$', fontsize=17)

############## mRNA-mRNA ##############
mRNA_covariances = np.genfromtxt("mRNA_mRNA_covariances_500.csv", delimiter=',')

mRNA_stds = [np.sqrt(mRNA_covariances[i,i] - R_discrete[i]**2) for i in range(len(R_discrete))]

PCCs = np.zeros(mRNA_covariances.shape)
for i in range(PCCs.shape[0]):
    for j in range(i):
        PCCs[i,j] = PCCs[j,i] = (mRNA_covariances[i,j] - R_discrete[i]*R_discrete[j])/(mRNA_stds[i]*mRNA_stds[j])

# Calculate extent
x_min = 0
x_max = PCCs.shape[1] * Dendrite_length/N_dendritic_segments
y_min = 0
y_max = PCCs.shape[0] * Dendrite_length/N_dendritic_segments

im = axs[2,0].imshow(PCCs, extent=[x_min, x_max, y_min, y_max], origin='lower', aspect='equal')

fig.colorbar(im, ax=axs[2,0])

### mRNA noise
G2_m2 = (2*lambda_2*gene_deactivation_rate*tau_1/(v_m + np.sqrt(v_m**2+4*D_m*(1/tau_1+1/tau_2))) + R[0])*R #(2*lambda_2*n_exp)/(v_m + np.sqrt(v_m**2+4*D_m/tau_2)) * np.exp(2*lambda_m * xi)

axs[1,1].plot(xi, G2_m2, label=r"$G^{(2)}_{m^2}$")

G2_m2_discrete = np.array([mRNA_covariances[i,i] - R_discrete[i] for i in range(len(R_discrete))])

axs[1,1].plot(np.arange(0,len(G2_m2_discrete)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), G2_m2_discrete*(N_dendritic_segments/Dendrite_length)**2, label="G2_discrete")

# G2_m2_discrete_ = np.array([mRNA_mRNA_covariances[i,i] - mRNA_expectations[i] for i in range(len(mRNA_expectations))])
# axs[1,1].plot(np.arange(0,len(G2_m2_discrete_)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), G2_m2_discrete_*(N_dendritic_segments/Dendrite_length)**2, label="G2_discrete_matrix_inversion")


axs[1,1].set_ylabel(r"$G^{(2)}_{m^2}$ density, " + r'$\mu m^{-3}$', fontsize=17)

######### mRNA-prot ###########
prot_covariances = np.genfromtxt("prot_prot_covariances_500.csv", delimiter=',')
mRNA_prot_covariances_ = np.genfromtxt("mRNA_prot_covariances_500.csv", delimiter=',')

prot_stds = [np.sqrt(prot_covariances[i,i] - P_discrete[i]**2) for i in range(len(P_discrete))]

PCCs = np.zeros(prot_covariances.shape)
for i in range(PCCs.shape[0]):
    PCCs[i,i] = (mRNA_prot_covariances_[i,i] - R_discrete[i]*P_discrete[i])/(mRNA_stds[i]*prot_stds[i])
    for j in range(i):
        PCCs[i,j] = PCCs[j,i] = (mRNA_prot_covariances_[i,j] - R_discrete[i]*P_discrete[j])/(mRNA_stds[i]*prot_stds[j])

# Calculate extent
x_min = 0
x_max = PCCs.shape[1] * Dendrite_length/N_dendritic_segments
y_min = 0
y_max = PCCs.shape[0] * Dendrite_length/N_dendritic_segments

im = axs[2,1].imshow(PCCs, extent=[x_min, x_max, y_min, y_max], origin='lower', aspect='equal')

fig.colorbar(im, ax=axs[2,1])

axs[2,1].set_xlabel("Distance from the soma, " + r'$\mu m$', fontsize=17)

### gene-prot correlations
G2_np_discrete = np.genfromtxt("gene_prot_covariances_500.csv", delimiter=',')

axs[0,2].plot(np.arange(0,len(G2_np_discrete)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), G2_np_discrete*N_dendritic_segments/Dendrite_length, label="G2_np_discrete")

# axs[0,2].plot(np.arange(0,len(gene_prot_covariances)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), gene_prot_covariances *N_dendritic_segments/Dendrite_length, label="G2_np_discrete_")

axs[0,2].set_ylabel(r"$G^{(2)}_{np}$ density, " + r'$\mu m^{-3}$', fontsize=17)

### prot-prot noise ###
G2_p2_discrete = np.array([prot_covariances[i,i]-P_discrete[i] for i in range(len(P_discrete))])
axs[1,2].plot(np.arange(0,len(G2_p2_discrete)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), G2_p2_discrete *(N_dendritic_segments/Dendrite_length)**2, label="G2_discrete")

# G2_p2_discrete_ = np.array([prot_prot_covariances[i,i]-prot_expectations[i] for i in range(len(prot_expectations))])

# axs[1,2].plot(np.arange(0,len(G2_p2_discrete_)*Dendrite_length/N_dendritic_segments-.001,Dendrite_length/N_dendritic_segments), G2_p2_discrete_*(N_dendritic_segments/Dendrite_length)**2, label="G2_discrete_matrix_inversion")

axs[1,2].set_ylabel(r"$G^{(2)}_{p^2}$ density, " + r'$\mu m^{-3}$', fontsize=17)

for ax in axs[:2,2]:
    ax.set_xlim([0,x_lim])
    ax.axhline(0, color='black')
    ax.set_yscale('log')
    ax.legend()

########################################################################################
PCCs = np.zeros(prot_covariances.shape)
for i in range(PCCs.shape[0]):
    for j in range(i):
        PCCs[i,j] = PCCs[j,i] = (prot_covariances[i,j] - P_discrete[i]*P_discrete[j])/(prot_stds[i]*prot_stds[j])

# Calculate extent
x_min = 0
x_max = PCCs.shape[1] * Dendrite_length/N_dendritic_segments
y_min = 0
y_max = PCCs.shape[0] * Dendrite_length/N_dendritic_segments

im = axs[2,2].imshow(PCCs, extent=[x_min, x_max, y_min, y_max], origin='lower', aspect='equal')

fig.colorbar(im, ax=axs[2,2])

for ax in axs[:2,1]:
    ax.set_xlim([0,x_lim])
    ax.axhline(0, color='black')
    ax.set_yscale('log')
    ax.legend()

axs[2,2].set_xlabel("Distance from the soma, " + r'$\mu m$', fontsize=17)

plt.tight_layout()
plt.show()
