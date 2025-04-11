import os
import matplotlib as plt
import numpy as np
import seaborn as sns
sns.set(style="ticks", palette="muted", 
        rc={"xtick.bottom" : True, "ytick.left" : True})
plt.style.use(r'./matplotlibrc.txt')

def kiel_diagram(models_folder,
                 model_name, ax,
                 ls='-', c='k', lw=2,
                 label='model'):
        
    model_folder = os.path.join(models_folder, model_name)
    
    Teff = np.genfromtxt(os.path.join(model_folder, f"{model_name}_value_log_Teff.txt"))
    logg = np.genfromtxt(os.path.join(model_folder, f"{model_name}_value_log_g.txt"))
    center_h1 = np.genfromtxt(os.path.join(model_folder, f"{model_name}_value_center_h1.txt"))

    pre_ms_mask = (center_h1 > 0.999*max(center_h1))

    Teff = Teff[pre_ms_mask]
    logg = logg[pre_ms_mask]

    if label is not None:
        ax.plot(Teff, logg, 
                color=c, linewidth=lw, linestyle=ls,
                label=label)
    else:
        ax.plot(Teff, logg, 
                color=c, linewidth=lw, linestyle=ls)
    
    ax.set_xlabel(r'$\log_{10}{T_{\rm{eff}} \ [\rm{K}]}$')
    ax.set_ylabel(r'$\log_{10}{\rm{g} \ [dex]}$')
    ax.invert_xaxis()