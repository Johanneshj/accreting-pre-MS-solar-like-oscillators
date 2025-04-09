import os
import matplotlib as plt
import matplotlib.cm as cm
import numpy as np
import seaborn as sns
import pandas as pd

sns.set(style="ticks", palette="muted", 
        rc={"xtick.bottom" : True, "ytick.left" : True})
plt.style.use(r'./matplotlibrc.txt')

def CD_diagram(models_folder,
               model_name, ax, label,
               classII=True,
               CD_folder='CD_diagram'):
    
    def get_frequencies(profile_number, freq_dir):
        path = os.path.join(freq_dir, 'profile' + str(profile_number) + '-freqs.dat')
        return pd.read_table(path, skiprows=5, sep='\s+')

    model_folder = os.path.join(models_folder, model_name)

    hist_model_number = np.genfromtxt(os.path.join(model_folder, f"{model_folder}_value_model_number.txt"))

    if classII==True:
        masses = np.genfromtxt(os.path.join(model_folder, f"{model_folder}_value_star_mass.txt"))
        center_h1 = np.genfromtxt(os.path.join(model_folder, f"{model_folder}_value_center_h1.txt"))
        mass_idx = np.where(masses >= 0.9*max(masses))[0]
        prems_idx = np.where(center_h1 >= 0.95*max(center_h1))[0]
        classII_idx = np.arange(mass_idx[0], prems_idx[-1])
        hist_model_number = hist_model_number[classII_idx]
    
    profile_folder = os.path.join(model_folder, "profile_data")
    index_file = np.genfromtxt(os.path.join(profile_folder, "index_file.txt"))
    model_numbers = np.array(index_file[:, 1], dtype=int)
    profile_numbers = np.array(index_file[:, 0], dtype=int)

    idx = np.isin(model_numbers, hist_model_number)
    profile_numbers = profile_numbers[idx]

    large_sep, small_sep = [], []
    for profile_number in enumerate(profile_numbers):
        freq_dir = os.path.join(f'{CD_folder}', f'{model_folder}')
        freq = get_frequencies(profile_number, freq_dir)

        l0 = freq[np.logical_and(freq.l == 0, freq.n_p > 10)].copy()
        l2 = freq[np.logical_and(freq.l == 2, freq.n_p > 10)].copy()
        l2.n_p += 1 

        freq_merge = pd.merge(l0, l2, on='n_p')
        freq_merge['d02'] = freq_merge['Re(freq)_x'] - freq_merge['Re(freq)_y']
        
        large_sep += [np.mean(np.diff(l0['Re(freq)'].values))]
        small_sep += [np.mean(freq_merge['d02'])]
    
    ax.plot(large_sep, small_sep, label=label)
    ax.set_xlabel(r'$\Delta\nu \ [\mu$Hz]')
    ax.set_ylabel(r'$\delta\nu_{02} \ [\mu$Hz]')

        