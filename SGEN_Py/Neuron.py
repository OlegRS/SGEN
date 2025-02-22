import sys
sys.path.append('../build')

import _SGEN_Py as _sg
import numpy as np

class Neuron:
    def __init__(self, soma, name):
        """Wrapper around the C++ Neuron class with additional plotting."""
        self._neuron = _sg._Neuron(soma, name)
        self._analytic_engine = _sg._Analytic_engine(self._neuron)
        self._morphologic_engine = _sg._Morphologic_engine(self._neuron)
        self._gillespie_engine = _sg._Gillespie_engine(self._neuron)

    def __init__(self, file_name, name='no_name'):
        """Wrapper around the C++ Neuron class with additional plotting."""
        self._neuron = _sg._Neuron(file_name, name)
        self._analytic_engine = _sg._Analytic_engine(self._neuron)
        self._morphologic_engine = _sg._Morphologic_engine(self._neuron)
        self._gillespie_engine = _sg._Gillespie_engine(self._neuron)


    # Analytics
    ### Stationary expectations
    def active_gene_expectation(self):
        return self._analytic_engine.active_genes_expectation()
    
    def mRNA_expectations(self, dict_return=False):
        mRNA_expectations = self._analytic_engine.stationary_mRNA_expectations()
        if dict_return:
            return dict(zip(self._analytic_engine.o1_mRNA_names(), mRNA_expectations))
        else:
            return np.array(mRNA_expectations)
    
    def protein_expectations(self, dict_return=False):
        prot_expectations = self._analytic_engine.stationary_protein_expectations()
        if dict_return:
            return dict(zip(self._analytic_engine.o1_prot_names(), prot_expectations))
        else:
            return np.array(prot_expectations)

    def expected_counts(self, dict_return=False):
        self._analytic_engine = _sg._Analytic_engine(self._neuron) # DEBUG (Temporary fix)
        return {'gene': self.active_gene_expectation(),
                'mRNA': self.mRNA_expectations(dict_return=dict_return),
                'prot': self.protein_expectations(dict_return=dict_return)}

    def mRNA_time_scales(self):
        return 1/np.array(self._analytic_engine.mRNA_o1_eigenvalues())

    def protein_time_scales(self):
        return 1/np.array(self._analytic_engine.protein_o1_eigenvalues())


    ### Stationary covariances (memory intense)
    def gene_gene_correlation(self):
        print('Computing gene-gene correlations...')
        return self._analytic_engine.stationary_gene_gene_covariance()
    def gene_mRNA_correlations(self):
        print('Computing gene-mRNA correlations...')
        return np.array(self._analytic_engine.stationary_gene_mRNA_covariances())
    def mRNA_mRNA_correlations(self):
        print('Computing mRNA-mRNA correlations...')
        return np.array(self._analytic_engine.stationary_mRNA_mRNA_covariances())
    def gene_protein_correlations(self):
        print('Computing gene-protein correlations...')
        return np.array(self._analytic_engine.stationary_gene_protein_covariances())
    def mRNA_protein_correlations(self):
        print('Computing mRNA-protein correlations...')
        return np.array(self._analytic_engine.stationary_mRNA_protein_covariances())
    def protein_protein_correlations(self):
        print('Computing protein-protein correlations...')
        return np.array(self._analytic_engine.stationary_protein_protein_covariances())

    def correlations(self):
        return {'gene-gene': self.gene_gene_correlation(),
                'gene-mRNA': self.gene_mRNA_correlations(),
                'mRNA-mRNA': self.mRNA_mRNA_correlations(),
                'gene-prot': self.gene_protein_correlations(),
                'mRNA-prot': self.mRNA_protein_correlations(),
                'prot-prot': self.protein_protein_correlations()}

    def expectations_and_correlations(self):
        return np.array(self._analytic_engine.stationary_expectations_and_correlations())

    def mRNA_mRNA_Pearson_correlation(self, comp1, comp2):
        return (self._analytic_engine.mRNA_mRNA_correlation(comp1,comp2) - self._analytic_engine.mRNA_expectation(comp1)*self._analytic_engine.mRNA_expectation(comp2))/(np.sqrt(self._analytic_engine.mRNA_mRNA_correlation(comp1,comp1)-self._analytic_engine.mRNA_expectation(comp1)**2)*np.sqrt(self._analytic_engine.mRNA_mRNA_correlation(comp2,comp2)-self._analytic_engine.mRNA_expectation(comp2)**2))

    def mRNA_prot_Pearson_correlation(self, comp1, comp2):
        return (self._analytic_engine.mRNA_protein_correlation(comp1,comp2) - self._analytic_engine.mRNA_expectation(comp1)*self._analytic_engine.protein_expectation(comp2))/(np.sqrt(self._analytic_engine.mRNA_mRNA_correlation(comp1,comp1)-self._analytic_engine.mRNA_expectation(comp1)**2)*np.sqrt(self._analytic_engine.protein_protein_correlation(comp2,comp2)-self._analytic_engine.protein_expectation(comp2)**2))

    def prot_prot_covariance(self, comp1, comp2):
        return self._analytic_engine.protein_protein_correlation(comp1,comp2) - self._analytic_engine.protein_expectation(comp1)*self._analytic_engine.protein_expectation(comp2)

    
    def prot_prot_Pearson_correlation(self, comp1, comp2):
        return self.prot_prot_covariance(comp1, comp2)/(np.sqrt(self._analytic_engine.protein_protein_correlation(comp1,comp1)-self._analytic_engine.protein_expectation(comp1)**2)*np.sqrt(self._analytic_engine.protein_protein_correlation(comp2,comp2)-self._analytic_engine.protein_expectation(comp2)**2))

        
    def protein_variance(self, comp):
        return self._analytic_engine.protein_protein_correlation(comp, comp) - self._analytic_engine.protein_expectation(comp)**2

    def protein_standard_deviation(self, comp):
        return np.sqrt(self._analytic_engine.protein_protein_correlation(comp, comp) - self._analytic_engine.protein_expectation(comp)**2)


    # Simulation
    def Gillespie_sim(self, record_times, output_file_name='', n_avrg_trajectories=1, burn_in=0):
        if n_avrg_trajectories < 1:
            raise ValueError("ERROR in run_Gillespie: Number of trajectories should be greater than zero")
        var_names = ['time'] + self._gillespie_engine.variable_names()
        
        print("Running Gillespie...")
        if n_avrg_trajectories == 1:
            var_values = np.array(self._gillespie_engine.run_Gillespie(record_times, output_file_name+'.csv', burn_in))
            var_dict = {var_names[i]: var_values[:,i] for i in range(len(var_names))}
        elif n_avrg_trajectories > 1:
            if output_file_name == '':
                var_values = np.array(self._gillespie_engine.run_Gillespie(record_times=record_times, burn_in=burn_in))[:,1:]/n_avrg_trajectories
                for tn in range(1, n_avrg_trajectories):
                    var_values += np.array(self._gillespie_engine.run_Gillespie(record_times=record_times, burn_in=burn_in))[:,1:]/n_avrg_trajectories
            else:
                var_values = np.array(self._gillespie_engine.run_Gillespie(record_times, output_file_name + '_1.csv', burn_in))[:,1:]/n_avrg_trajectories
                for tn in range(2, n_avrg_trajectories+1):
                    var_values += np.array(self._gillespie_engine.run_Gillespie(record_times, output_file_name + '_' + str(tn) + '.csv', burn_in))[:,1:]/n_avrg_trajectories

            var_dict = {var_names[i]: var_values[:,i-1] for i in range(1,len(var_names))}
            var_dict.update({'time': np.array(record_times)})
        return var_dict

    def stationary_Gillespie_sim(self, record_times, output_file_name='', n_avrg_trajectories=1, burn_in_factor=5):
        gene_timescale = max(1/self.soma().gene_activation_rate(), 1/self.soma().gene_deactivation_rate())
        burn_in = burn_in_factor*max(gene_timescale, max(self.mRNA_time_scales()), max(self.protein_time_scales()))
        print("Gillespie burn-in time: ", burn_in)
        return self.run_Gillespie(record_times, output_file_name, n_avrg_trajectories, burn_in)
    
    def load_Gillespie_sim(self, file_name):
        with open(file_name, "r") as f:
            var_names = f.readline().strip().split(",")  # Read and split header
            var_values = np.loadtxt(f, delimiter=",", dtype=float)  # Load numeric data
        return {var_names[i]: var_values[:,i] for i in range(len(var_names))}


    # Mrphology
    def segments(self):
        return np.array(_sg._Morphologic_engine(self._neuron).segments())
    
    def volumes(self):
        return np.array(_sg._Morphologic_engine(self._neuron).volumes())

    def soma(self):
        return self._neuron.soma()


    # Plotting 
    def draw_3d(self, visualisation_values=[]):
        import pyvista as pv

        segments = self.segments()
        
        start_points = [segments[i][0][:3] for i in range(len(segments))]
        end_points = [segments[i][1][:3] for i in range(len(segments))]
        radii = [segments[i][1][3] for i in range(len(segments))]

        # Make a tube for each segment to represent dendrites
        tubes = []
        if visualisation_values != []:
            visualisation_values = np.flip(visualisation_values)
            all_scalars = [] # Collect scalars for the final mesh
        for i in range(0,len(segments)):
            start, end = start_points[i], end_points[i]
            radius = radii[i]  # Radius of the current segment

            # Create a tube (cylinder) between two coordinates
            tube = pv.Line(start, end).tube(radius=radius)
            tubes.append(tube)

            if visualisation_values != []:
                # Extend scalars to match the number of points in the current tube
                all_scalars.append(np.full(tube.n_cells, visualisation_values[i]))

        # Combine all tubes into a single mesh
        neuron_mesh = tubes[0] if tubes else None
        for tube in tubes[1:]:
            neuron_mesh += tube

        if visualisation_values != []:
            # Concatenate scalars for all segments
            flat_scalars = np.concatenate(all_scalars)
            # Assign scalars to the combined mesh
            neuron_mesh.cell_data["Protein Levels"] = flat_scalars

        # Visualize the neuron
        plotter = pv.Plotter()
        if visualisation_values != []:
            plotter.add_mesh(neuron_mesh, scalars="Protein Levels", cmap="coolwarm", show_edges=False)
        else:
            plotter.add_mesh(neuron_mesh, cmap="coolwarm", show_edges=False)
        plotter.show_axes()
        plotter.show()



# # Prevent users from instantiating Cpp_Neuron directly
# def _raise_error(*args, **kwargs):
#     raise TypeError("Direct instantiation of Cpp_Neuron is not allowed. Use Neuron instead!")

# _sg.Cpp_Neuron.__new__ = _raise_error  # Override constructor
