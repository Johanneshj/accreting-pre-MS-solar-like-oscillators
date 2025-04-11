import os
import matplotlib as plt
import numpy as np
import seaborn as sns
sns.set(style="ticks", palette="muted", 
        rc={"xtick.bottom" : True, "ytick.left" : True})
plt.style.use(r'./matplotlibrc.txt')

def HR_diagram(models_folder,
               model_name, ax,
               ls='-', c='k', lw=2,
               label='model'):
    
    model_folder = os.path.join(models_folder, model_name)
    
    Teff = np.genfromtxt(os.path.join(model_folder, f"{model_name}_value_log_Teff.txt"))
    Lum = np.genfromtxt(os.path.join(model_folder, f"{model_name}_value_log_L.txt"))
    center_h1 = np.genfromtxt(os.path.join(model_folder, f"{model_name}_value_center_h1.txt"))

    pre_ms_mask = (center_h1 > 0.999*max(center_h1))

    Teff = Teff[pre_ms_mask]
    Lum = Lum[pre_ms_mask]

    if label is not None:
        ax.plot(Teff, Lum, 
                color=c, linewidth=lw, linestyle=ls,
                label=label)
    else:
        ax.plot(Teff, Lum, 
                color=c, linewidth=lw, linestyle=ls)
    
    ax.set_xlabel(r'$\log_{10}{T_{\rm{eff}} \ [\rm{K}]}$')
    ax.set_ylabel(r'$\log_{10}{\rm{luminosity} \ [L_\odot]}$')
    ax.invert_xaxis()