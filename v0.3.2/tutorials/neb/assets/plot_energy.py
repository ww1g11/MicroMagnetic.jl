import matplotlib as mpl
mpl.use("Agg")
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 12,
        }

def custom_legend(legend):
    frame = legend.get_frame()
    frame.set_facecolor('0.90')

    # Set the fontsize
    for label in legend.get_texts():
        label.set_fontsize(10)

    for label in legend.get_lines():
        label.set_linewidth(1.5)  # the legend line width

def load_data_simple(name):
    
    energy = np.loadtxt('%s_energy.txt' %(name))
    distance = np.loadtxt('%s_distance.txt' %(name))
    distance[:,0] = 0
    eV = 1.602176565e-19
    energy = (energy - np.min(energy))/eV
    
    return energy[:,1:], np.cumsum(distance[:,:], axis=1)

def plot_energy(name):
    energy_data, distances = load_data_simple(name)

    plt.rc('font', **font)
    fig = plt.figure(figsize=(5,4))

    xnew = np.linspace(distances[-1, 0], distances[-1, -1], num=101, endpoint=True)
    f = interp1d(distances[-1, :], energy_data[-1, :], kind='cubic')
    p = plt.plot(distances[-1, :], energy_data[-1, :], 'o', color='slateblue', label='Minimum energy path')
    plt.plot(xnew, f(xnew), '-', c='slateblue')

    p = plt.legend(loc='best')
    custom_legend(p)
    plt.xlabel("Distances (arb. unit)")
    plt.ylabel("Energy (eV)")
    filename = '%s_energy_2d.png' % name
    plt.tight_layout()
    fig.savefig(filename, dpi=300)

plot_energy(name="skx_fm")
