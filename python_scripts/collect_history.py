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

mass = float('<<mass>>')
acchist = int('<<acchist>>')

def collect_history_data(mass, acchist):
    folder = 'grid'
    hist = np.genfromtxt(f"../../{folder}/mass-{mass}Msun_acchist-{acchist}/LOGS/history.data", 
                        skip_header = 5, 
                        names = True
                        )
    model_number =  hist["model_number"]
    sort=np.argsort(model_number)
    for key in keys:
        data = hist[key]
        np.savetxt(f"mass-{mass}Msun_acchist-{acchist}_value_{key}.txt", data[sort].T)
        print(f"mass-{mass}Msun_acchist-{acchist}_value_{key}")

collect_history_data(mass=mass, acchist=acchist)
