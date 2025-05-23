import sys
sys.path.append('../build')

import _SGEN_Py as _sg
import numpy as np
import pyvista as pv

class Neuron:
    def __init__(self, input_data, name='no_name'):
        """Wrapper around the C++ Neuron class with additional plotting."""
        if isinstance(input_data, str):
            import os
            if not os.path.exists(input_data):
                raise FileNotFoundError(f"Error: The file '{input_data}' does not exist.")
            self._neuron = _sg._Neuron(input_data, name)
        elif isinstance(input_data, _sg.Soma):
            self._neuron = _sg._Neuron(input_data, name)
        self._analytic_engine = _sg._Analytic_engine(self._neuron)
        self._morphologic_engine = _sg._Morphologic_engine(self._neuron)
        self._gillespie_engine = _sg._Gillespie_engine(self._neuron)
        
    def __str__(self):
        return str(self._neuron)

    # Analytics
    ### Stationary expectations
    def active_gene_expectation(self):
        return self._analytic_engine.active_genes_expectation()
    
    def mRNA_expectations(self, dict_return=False):
        mRNA_expectations = self._analytic_engine.stationary_mRNA_expectations()
        if dict_return:
            return dict(zip(self._analytic_engine.o1_prot_names(), mRNA_expectations))
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

    def mRNA_standard_deviation(self, comp):
        return np.sqrt(self._analytic_engine.mRNA_mRNA_correlation(comp, comp) - self._analytic_engine.mRNA_expectation(comp)**2)


    ### Nonstationary covariances
    def nonstationary_moments(self, record_times, file_name, dt_factor=.3):#, t_max_factor=3):
        ## Determining time step from characteristic time scales of the system
        # gene_long_timescale = max(1/self.soma().gene_activation_rate(), 1/self.soma().gene_deactivation_rate())
        gene_short_timescale = min(1/self.soma().gene_activation_rate(), 1/self.soma().gene_deactivation_rate())
        # max_time = t_max_factor * max(gene_long_timescale, max(self.mRNA_time_scales()), max(self.protein_time_scales()))
        dt = dt_factor * min(gene_short_timescale, min(self.mRNA_time_scales()), min(self.protein_time_scales()))

        print("nonstationary_moments dt=", dt)

        results = self._analytic_engine.nonstationary_moments(record_times, dt, file_name)
        output = {'gene' : np.array([results[t_ind][0][0][0] for t_ind in range(len(record_times))]),
                  'mRNA' : np.array([results[t_ind][1][0] for t_ind in range(len(record_times))]),
                  'prot' : np.array([results[t_ind][2][0] for t_ind in range(len(record_times))]),
                  'gene_gene' : np.array([results[t_ind][3][0][0] for t_ind in range(len(record_times))]),
                  'gene_mRNA' : np.array([results[t_ind][4][0] for t_ind in range(len(record_times))]),
                  'mRNA_mRNA' : np.array([results[t_ind][5] for t_ind in range(len(record_times))]),
                  'gene_prot' : np.array([results[t_ind][6][0] for t_ind in range(len(record_times))]),
                  'mRNA_prot' : np.array([results[t_ind][7] for t_ind in range(len(record_times))]),
                  'prot_prot' : np.array([results[t_ind][8] for t_ind in range(len(record_times))])}
        return output

    # Simulation
    def Gillespie_sim(self, record_times, output_file_name=None, n_avrg_trajectories=1, burn_in=0, reset=False):
        if(reset):
            self.reset()
        if n_avrg_trajectories < 1:
            raise ValueError("ERROR in run_Gillespie: Number of trajectories should be greater than zero")
        var_names = self._gillespie_engine.variable_names()
        
        print("Running Gillespie...")
        if n_avrg_trajectories == 1:
            if output_file_name is None:
                var_values = np.array(self._gillespie_engine.run_Gillespie(record_times=record_times, burn_in=burn_in))
            else:
                var_values = np.array(self._gillespie_engine.run_Gillespie(record_times, output_file_name+'.csv', burn_in))
            var_dict = {var_names[i]: var_values[:,i] for i in range(len(var_names))}
        elif n_avrg_trajectories > 1:
            if output_file_name is None:
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

    def stationary_Gillespie_sim(self, record_times, output_file_name=None, n_avrg_trajectories=1, burn_in_factor=5):
        gene_timescale = max(1/self.soma().gene_activation_rate(), 1/self.soma().gene_deactivation_rate())
        burn_in = burn_in_factor*max(gene_timescale, max(self.mRNA_time_scales()), max(self.protein_time_scales()))
        print("Gillespie burn-in time: ", burn_in)
        if output_file_name is None:
            return self.Gillespie_sim(record_times, n_avrg_trajectories=n_avrg_trajectories, burn_in=burn_in)
        return self.Gillespie_sim(record_times, output_file_name, n_avrg_trajectories, burn_in)
   
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

    def dendritic_segments(self):
        return self._neuron.dendritic_segments()

    def spines(self):
        return self._neuron.spines()

    def reset(self):
        self._gillespie_engine.reset()


    # Plotting
    def draw_3d(self, visualisation_values=None, color='#32CD32', file_name=None, plotter=None, scalar_bar_args=None, cbar_scale="linear", clim=None, mRNA_visualisation=False, mRNA_radius=None, mRNA_data=None, gene_visualisation=False, active_gene_count=None, gene_length=None, show_axes=False, camera_position=None, opacity=1, annotations=None, cmap="viridis"):
        if "google.colab" in sys.modules:
            # Seems that only static plotting is supported by colab at the moment
            pv.global_theme.jupyter_backend = 'static'
            pv.global_theme.notebook = True
            pv.start_xvfb()
            
        segments = self.segments()
        
        start_points = [segments[i][0][:3] for i in range(len(segments))]
        end_points = [segments[i][1][:3] for i in range(len(segments))]
        radii = [segments[i][1][3] for i in range(len(segments))]
        
        tubes = []
        all_scalars = []
        
        if visualisation_values is not None:
            if clim is None:
                clim = [visualisation_values.min(), visualisation_values.max()] if visualisation_values is not None else None
                
            visualisation_values = np.flip(visualisation_values)
            
            if cbar_scale == "uniform":
                expected_concentrations = self.expected_counts()["prot"] / self.volumes()
                expected_concentrations_appended = np.append(expected_concentrations, [clim[0] - 1e-7, clim[1] + 1e-7])
                sorted_ec = np.sort(expected_concentrations_appended)
                rescaled_pc = np.zeros(visualisation_values.shape)
                for i in range(visualisation_values.shape[0]):
                    m = np.searchsorted(sorted_ec, visualisation_values[i])-1 # Find index: sorted_ec[m] < visualisation_values[i] < sorted_ec[m+1]
                    rescaled_pc[i] = m + (visualisation_values[i] - sorted_ec[m])/(sorted_ec[m+1] - sorted_ec[m])

            elif cbar_scale == "log":
                rescaled_pc = np.log(visualisation_values)
            elif cbar_scale == "linear":
                rescaled_pc = visualisation_values
            else:
                raise ValueError("cbar_scale can be either \"uniform\", \"log\" or \"linear\"")

        for i in range(len(segments)):
            start, end = start_points[i], end_points[i]
            radius = radii[i]
            
            tube = pv.Line(start, end).tube(radius=radius)
            tubes.append(tube)
            
            if visualisation_values is not None:
                all_scalars.append(np.full(tube.n_cells, visualisation_values[i]))
                
        if not tubes:
            raise ValueError("No segments to visualise")

        neuron_mesh = tubes[0]
        for tube in tubes[1:]:
            neuron_mesh += tube

        SHOW_PLOTTER = False
        if plotter is None:
            plotter = pv.Plotter()
            SHOW_PLOTTER = True
        else:
            plotter.clear()
            
        if visualisation_values is not None:
            flat_scalars = np.concatenate(all_scalars)
            neuron_mesh.cell_data["Protein Levels"] = flat_scalars

            if scalar_bar_args is None:
                # Define colorbar parameters
                scalar_bar_args = {
                    "title": "Protein Concentration",
                    # "title_font_size": 18,
                    # "label_font_size": 14,
                    # "color": "black",
                    # "position_x": 0.8,   # Horizontal position in range [0,1]
                    # "position_y": 0.05,  # Vertical position in range [0,1]
                    # "width": 0.08,       # Width of the bar in figure coords
                    # "height": 0.9,       # Height of the bar in figure coords
                    "n_labels": 0,       # Number of tick labels (can be overridden below)
                    "fmt": "%.2e",       # Format string for tick labels
                    "vertical": False    # Orientation
                }

                
            if mRNA_visualisation:
                plotter.add_mesh(neuron_mesh, scalars="Protein Levels", cmap=cmap, clim=clim, show_edges=False, opacity=opacity, scalar_bar_args=scalar_bar_args, annotations=annotations)
            else:
                plotter.add_mesh(neuron_mesh, scalars="Protein Levels", cmap=cmap, clim=clim, show_edges=False, scalar_bar_args=scalar_bar_args, annotations=annotations)
        else:
            if mRNA_visualisation:
                plotter.add_mesh(neuron_mesh, color=color, show_edges=False, opacity=opacity, annotations=annotations)
            else:
                plotter.add_mesh(neuron_mesh, color=color, show_edges=False, opacity=opacity, annotations=annotations)

        if gene_visualisation:
            if gene_length is None:
               gene_length = min([self.soma().length(), self.soma().radius()])/2
            end_pos = self.soma().position()
            theta, phi = self.soma().orientation()
            r, l = self.soma().radius(), self.soma().length()
            n = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]) # Normal vector to the transversal plane
            n_gene_copies = self.soma().n_gene_copies()
            gene_radius = gene_length/(2*n_gene_copies)
            center = end_pos - ((l+gene_length)/2-gene_radius)*n
            
            if active_gene_count is None:
                for i in range(n_gene_copies):
                    if i < self.soma().n_active_genes():
                        plotter.add_mesh(pv.Sphere(radius=gene_length/(2*n_gene_copies), center=center), color="green", opacity=.8)
                    else:
                        plotter.add_mesh(pv.Sphere(radius=gene_length/(2*n_gene_copies), center=center), color="black", opacity=.8)
                    center += 2*gene_radius*n
            else:
                if n_gene_copies < active_gene_count:
                    raise ValueError("n_gene_copies < active_gene_count")
                for i in range(n_gene_copies):
                    if i < active_gene_count:
                        plotter.add_mesh(pv.Sphere(radius=gene_length/(2*n_gene_copies), center=center), color="green", opacity=.8)
                    else:
                        plotter.add_mesh(pv.Sphere(radius=gene_length/(2*n_gene_copies), center=center), color="black", opacity=.8)
                    center += 2*gene_radius*n
                    
        # Add mRNA molecules as small spheres
        if mRNA_visualisation:
            dendritic_segments = self._neuron.dendritic_segments()
            if mRNA_radius is None:
                mRNA_radius = min([ds.radius() for ds in dendritic_segments])/5
            if mRNA_data is None:
                end_pos = self.soma().position()
                theta, phi = self.soma().orientation()
                r, l = self.soma().radius(), self.soma().length()
                    
                for i in range(self.soma().mRNA_count()):
                    l_pos = np.random.uniform(0,l) # Axial displacement
                    rd = np.random.rand(3) # Random direction
                    n = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]) # Normal vector to the transversal plane
                    rdip =  rd - (rd@n)*n # Random direction in transversal plane
                    rdip /= np.linalg.norm(rdip)
                    td = np.random.uniform(mRNA_radius-r, r-mRNA_radius) # Transversal displacement
                    center = end_pos - l_pos*n + td*rdip
                    plotter.add_mesh(pv.Sphere(radius=mRNA_radius, center=center), color="black")
                        
                for ds in dendritic_segments:
                    end_pos = ds.position()
                    theta, phi = ds.orientation()
                    r, l = ds.radius(), ds.length()
                    for i in range(ds.mRNA_count()):
                        l_pos = np.random.uniform(0,l) # Axial displacement
                        rd = np.random.rand(3) # Random direction
                        n = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]) # Normal vector to the transversal plane
                        rdip =  rd - (rd@n)*n # Random direction in transversal plane
                        rdip /= np.linalg.norm(rdip)
                        td = np.random.uniform(mRNA_radius-r,r-mRNA_radius) # Transversal displacement
                        center = end_pos - l_pos*n + td*rdip
                        plotter.add_mesh(pv.Sphere(radius=mRNA_radius, center=center), color="black")
            else:
                mRNA_count = next((value for key, value in mRNA_data.items() if self.soma().name() in key), None)
                if mRNA_count is not None:
                    end_pos = self.soma().position()
                    theta, phi = self.soma().orientation()
                    r, l = self.soma().radius(), self.soma().length()
                    
                    for i in range(int(mRNA_count)):
                        l_pos = np.random.uniform(0,l) # Axial displacement
                        rd = np.random.rand(3) # Random direction
                        n = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]) # Normal vector to the transversal plane
                        rdip =  rd - (rd@n)*n # Random direction in transversal plane
                        rdip /= np.linalg.norm(rdip)
                        td = np.random.uniform(mRNA_radius-r,r-mRNA_radius) # Transversal displacement
                        center = end_pos - l_pos*n + td*rdip
                        plotter.add_mesh(pv.Sphere(radius=mRNA_radius, center=center), color="black")
                    
                for ds in dendritic_segments:
                    mRNA_count = next((value for key, value in mRNA_data.items() if ds.name() in key), None)
                    end_pos = ds.position()
                    theta, phi = ds.orientation()
                    r, l = ds.radius(), ds.length()
                    if mRNA_count is not None:
                        for i in range(int(mRNA_count)):
                            l_pos = np.random.uniform(0,l) # Axial displacement
                            rd = np.random.rand(3) # Random direction
                            n = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]) # Normal vector to the transversal plane
                            rdip =  rd - (rd@n)*n # Random direction in transversal plane
                            rdip /= np.linalg.norm(rdip)
                            td = np.random.uniform(mRNA_radius-r,r-mRNA_radius) # Transversal displacement
                            center = end_pos - l_pos*n + td*rdip
                            plotter.add_mesh(pv.Sphere(radius=mRNA_radius, center=center), color="black")
                    else:
                        raise ValueError(f"mRNA count not found for dendritic segment: {ds.name()}")
                    
        if show_axes:
            plotter.show_axes()

        if camera_position is not None:
            plotter.camera_position = camera_position
            
        if file_name is not None:
            plotter.window_size = [1500, 1500]
            plotter.enable_anti_aliasing()
            plotter.save_graphic(file_name)

        if SHOW_PLOTTER:
            plotter.show()
   
# # Prevent users from instantiating Cpp_Neuron directly
# def _raise_error(*args, **kwargs):
#     raise TypeError("Direct instantiation of Cpp_Neuron is not allowed. Use Neuron instead!")

# _sg.Cpp_Neuron.__new__ = _raise_error  # Override constructor
