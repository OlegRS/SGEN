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


    # Simulation
    def run_Gillespie(self, times, output_file_name, time_offset=0):
        return np.array(self._gillespie_engine.run_Gillespie(record_times, output_file_name, time_offset))


    # Mrphology
    def segments(self):
        return np.array(_sg._Morphologic_engine(self._neuron).segments())
    
    def volumes(self):
        return np.array(_sg._Morphologic_engine(self._neuron).volumes())


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
