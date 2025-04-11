import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set(style="ticks", palette="muted", 
        rc={"xtick.bottom" : True, "ytick.left" : True})
plt.style.use(r'./matplotlibrc.txt')
from HR_diagram import HR_diagram
from kiel_diagram_plot import kiel_diagram


colors = 'k'

fig, axs = plt.subplots(1,2,figsize=(10, 6))

for mass in [1.0]:
    for acchist in [30]:
        model_name = f'mass-{mass}Msun_acchist-{acchist}'
        
        HR_diagram(models_folder='data_resolution_test',
                    model_name=model_name, ax=axs[0],
                    ls='-', c=colors, lw=1.5,
                    label='')
        kiel_diagram(models_folder='data_resolution_test',
                    model_name=model_name, ax=axs[1],
                    ls='-', c=colors, lw=1.5,
                    label='')
            
axs[0].set_xlim(3.5, 3.85)
axs[0].set_ylim(-.5, 1.5)
axs[0].invert_xaxis()
axs[1].invert_xaxis()
axs[1].invert_yaxis()
fig.savefig('figures/resolution_test')
