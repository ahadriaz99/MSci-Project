'''
Analysis tool for visualising M7 FCIQMC output
Author: Marcell Dorian Kovacs
06/12/2021
'''
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

path_to_output = ''

input_file = path_to_output+'M7.stats'
params = {
   'axes.labelsize': 40,
   'font.size': 35,
   'legend.fontsize': 35,
   'lines.linewidth' : 2,
   'lines.markersize' : 35,
   'xtick.labelsize': 35,
   'ytick.labelsize': 35,
   'figure.figsize': [40, 20]
   }
plt.rcParams.update(params)
class Data:
    def __init__(self, N, M, L, path=''):
    
        assert L >= 0 and L <= N*(M-1)
        assert M >= 0 and N >= 0
    
        self.N = N
        self.M = M
        self.L = L
        
        self.path = path
        
        self.cycles = []
        self.timesteps = []
        self.diag_shift = []
        self.nwalkers = []
        self.energy_nom = []
        self.normalisation = []
        
        self.timestep = None
        self.ncycle = None
        self.nwalker = None
        
        self.E0 = None
        self.E0_std_err = None
        self.E0_exact = None
                
    def read_output(self):
        
        file = self.path+'M7.stats'
        with open(file, 'r') as f:
            for line in f:
                split_line = line.split()
                if (split_line[0][0] == '#'):
                    continue
                self.cycles.append(float(split_line[0]))
                self.timesteps.append(float(split_line[1]))
                
                # Keep track of total walker population at each timestep
                self.nwalkers.append(float(split_line[3]))
                
                self.diag_shift.append(float(split_line[2]))
                self.energy_nom.append(float(split_line[7]))
                self.normalisation.append(float(split_line[8]))
                
        self.timestep = self.timesteps[1]
        self.ncycle = len(self.cycles)
        self.nwalker = np.max(self.nwalkers)
                
    def fig_population_evolution(self):
    
        plt.title('Total walker population evolution \n'+\
                  'Disc geometry: No. bosons %d   No. Landau levels %d  Total ang. mom. %d \n'%(self.N, self.M, self.L)+\
                  r'Timestep $\tau$ = %1.7f'%(self.timestep))
        plt.ylabel('Total walker population')
        plt.xlabel(r'Timesteps [$\tau$]')
        plt.plot(self.cycles, self.nwalkers)
        plt.grid()
        plt.savefig('Disc_N%d_M%d_L%d_Walker_Evol.jpeg'%(self.N, self.M, self.L))
        plt.close()
        #for i in range(len(self.cycles)):
        #    print(self.cycles[i], self.nwalkers[i])
        
    def fig_shift_evolution(self):
    
        plt.title('Diagonal shift evolution \n'+\
                  'Disc geometry: No. bosons %d   No. Landau levels %d  Total ang. mom. %d \n'%(self.N, self.M, self.L)+\
                  r'Timestep $\tau$ = %1.7f'%(self.timestep))
        plt.ylabel('Diagonal shift [$V_0$]')
        plt.xlabel(r'Timesteps [$\tau$]')
        plt.plot(self.cycles, self.diag_shift)
        plt.grid()
        plt.savefig('Disc_N%d_M%d_L%d_Shift_Evol.jpeg'%(self.N, self.M, self.L))
        plt.close()
        #for i in range(len(self.cycles)):
        #    print(self.cycles[i], self.nwalkers[i])
        
    def fig_proj_energy_evolution(self, cshift=0):
    
        self.e_proj = np.array(self.energy_nom)/np.array(self.normalisation)
    
        plt.title('Projected energy evolution \n'+\
                  'Disc geometry: No. bosons %d   No. Landau levels %d  Total ang. mom. %d \n'%(self.N, self.M, self.L)+\
                  r'Timestep $\tau$ = %1.7f'%(self.timestep))
        plt.ylabel('Projected energy [$V_0$]')
        plt.xlabel(r'Timesteps [$\tau$]')
        plt.plot(self.cycles[cshift:], self.e_proj[cshift:])
        if (self.E0 is not None):
            plt.plot(self.cycles[cshift:], self.E0*np.ones(len(self.cycles))[cshift:], color='red', linewidth=3.5, label='QMC $E_0 \; [V_0] = %1.7f \pm %1.7f $'%(self.E0, self.E0_std_err))
        if (self.E0_exact is not None):
            plt.plot(self.cycles[cshift:], self.E0_exact*np.ones(len(self.cycles))[cshift:],linewidth=3.5, color='green', label='FCI    $E_0 \; [V_0] = $%1.6f'%self.E0_exact)
        plt.grid()
        plt.legend()
        if (cshift == 0):
            plt.savefig('Disc_N%d_M%d_L%d_Eproj_Evol.jpeg'%(self.N, self.M, self.L))
        else:
            plt.savefig('Disc_N%d_M%d_L%d_Eproj_Evol_zoom.jpeg'%(self.N, self.M, self.L))
        plt.close()
        #for i in range(len(self.cycles)):
        #    print(self.cycles[i], self.nwalkers[i])

    def exact_energy(self, E0_exact):
        self.E0_exact = E0_exact
        
    def energy_analysis(self, E0, E0_std_err):
        self.E0 = E0
        self.E0_std_err = E0_std_err
        

output = Data(5, 5, 15, path="Disc_N5_M5_L15/")
output.read_output()
output.exact_energy(2.10692052)
output.energy_analysis(2.107056984584, 1.025373938055e-04)
output.fig_population_evolution()
output.fig_shift_evolution()
output.fig_proj_energy_evolution()
output.fig_proj_energy_evolution(2000)


output = Data(6, 8, 14, path="Disc_N6_M8_L14/")
output.read_output()
output.exact_energy(2.69433994)
output.energy_analysis(2.691843511287, 4.336893497106e-04)

output.fig_population_evolution()
output.fig_shift_evolution()
output.fig_proj_energy_evolution()
output.fig_proj_energy_evolution(2000)

