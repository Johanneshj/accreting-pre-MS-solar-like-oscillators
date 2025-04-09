### Inspired by
### https://github.com/earlbellinger/mesa-gyre-tutorial-2022/blob/main/tutorial.ipynb

import os
import matplotlib as plt
import numpy as np

def propagation_diagram(model_name, profile_number, ax,
                        brunt_color='k', lamb_color='red', lamb_ls='-', 
                        brunt_filled=True, brunt_hatched=False,
                        plot_log_scale=True):

    profile_folder = os.path.join(model_name, "profile_data")

    N2_folder = os.path.join(profile_folder, "brunt_N2")
    logR_folder = os.path.join(profile_folder, "logR")
    lamb_Sl1_folder = os.path.join(profile_folder, "lamb_Sl1")

    # Load Brunt and Lamb profiles and the radius
    brunt_N2 = np.genfromtxt(os.path.join(N2_folder, f"brunt_N2-{profile_number}.txt"))
    lamb_Sl1 = np.genfromtxt(os.path.join(lamb_Sl1_folder, f"lamb_Sl1-{profile_number}.txt"))
    logR = np.genfromtxt(os.path.join(logR_folder, f"logR-{profile_number}.txt"))
    radius = 10**logR

    # ------ PLOT ------ #
    # Buoyancy
    ax.plot(radius, np.sqrt(abs(brunt_N2))/(2*np.pi)*1e6, 
            ls='-', color=brunt_color)
    
    if brunt_filled:
        ax.fill_between(radius, 
                        np.zeros(len(np.sqrt(abs(brunt_N2))/(2*np.pi)*1e6)), 
                        np.sqrt(abs(brunt_N2))/(2*np.pi)*1e6, 
                        color=brunt_color, alpha=0.6, zorder=-100)
    elif brunt_hatched:
        ax.fill_between(radius, 
                        np.zeros(len(np.sqrt(abs(brunt_N2))/(2*np.pi)*1e6)), 
                        np.sqrt(abs(brunt_N2))/(2*np.pi)*1e6, color='none', 
                        edgecolor='k', hatch='///', zorder=-100)
    
    # Lamb l1
    ax.plot(radius, lamb_Sl1, 
            ls=lamb_ls, color=lamb_color, zorder=-150)
    
    # Labels and limits
    ax.set_xlabel(r'radius [R$_\odot$]')
    ax.set_ylabel(r'frequency $\nu \ [\mu\rm{Hz}]$')

    if plot_log_scale:   
        ax.set_yscale('log')
        ax.set_xscale('log')                     
        ax.set_xlim(0.01, np.max(radius))
        ax.set_ylim(2, 1e5)
        
        # Ticks
        xticks = [1e-2, 1e-1, 0.5, 1, np.round(np.max(radius), 1)]
        ax.set_xticks(xticks)
        ax.set_xticklabels(['0.01', '0.1', '0.5', '1', str(np.round(np.max(radius), 1))])
    else:
        ax.set_xlim(0.0, np.max(radius))
        ax.set_yscale('log')
        ax.set_ylim(2, 1e5)
    
