import numpy as np
import os

keys = ["log_Teff", "log_g", "log_L", "log_R", "model_number", "star_age", "star_mass", 
        "acoustic_cutoff", "nu_max", "center_h1"]

for i in range(1, 21):
    string1 = "mix_qtop_"+str(i)
    keys.append(string1)
    string2 = "mix_type_"+str(i)
    keys.append(string2)
    string3 = "burn_qtop_"+str(i)
    keys.append(string3)
    string4 = "burn_type_"+str(i)
    keys.append(string4)
    string5 = "mix_relr_top_"+str(i)
    keys.append(string5)
    string6 = "mix_relr_type_"+str(i)
    keys.append(string6)
    string7 = "burn_relr_top_"+str(i)
    keys.append(string7)
    string8 = "burn_relr_type_"+str(i)
    keys.append(string8)

keys = ["log_Teff", "log_g", "log_L", "center_h1"]

def collect_history_data(mass, acchist,
                        grid_folder='',
                        dest_folder=''):
    folder = grid_folder
    hist = np.genfromtxt(f"{folder}/mass-{mass}Msun_acchist-{acchist}/LOGS/history.data", 
                        skip_header = 5, 
                        names = True
                        )
    model_number =  hist["model_number"]
    sort=np.argsort(model_number)

    for key in keys:
        data = hist[key]
        np.savetxt(f"{dest_folder}/mass-{mass}Msun_acchist-{acchist}_value_{key}.txt", data[sort].T)
        print(f"mass-{mass}Msun_acchist-{acchist}_value_{key}")

# ------------------------------------------------ #
'''prefix = ''
for mass in [1.0]:
    for acchist in [30]:
        
        folder = '../resolution_test'
        data_folder = 'data_resolution_test'
        if not os.path.exists(data_folder):
            os.mkdir(data_folder)

        name = f"{prefix}mass-{mass}Msun_acchist-{acchist}"

        if not os.path.exists(os.path.join(data_folder, name)):
            os.mkdir(os.path.join(data_folder, name))
            dest_folder = os.path.join(data_folder, name)
        else:
            dest_folder = os.path.join(data_folder, name)
        
        collect_history_data(mass, acchist,
                            grid_folder=folder,
                            dest_folder = dest_folder)'''