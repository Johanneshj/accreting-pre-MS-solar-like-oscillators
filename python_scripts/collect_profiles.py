import numpy as np
import os

profile_keys = ['brunt_N2', 'lamb_Sl1', 'mass', 'logR']

mass = float('<<mass>>')
acchist = int('<<acchist>>')

def collect_profile_data(mass, acchist):
    folder = 'grid'

    index_table = np.genfromtxt(os.path.join(f"../../../{folder}/mass-{mass}Msun_acchist-{acchist}/LOGS/", 'profiles.index'), 
                skip_header = 1, names = ['model_number', 'priority', 'profile_number'], dtype=int)

    profile_numbers = index_table['profile_number']
    model_numbers = index_table['model_number']

    index_data = np.array([profile_numbers, model_numbers], dtype=int)
    np.savetxt(f"index_file.txt", index_data.T)

    for key in profile_keys:
        if not os.path.exists(f"{key}"):
            os.mkdir(f"{key}")

        for profile_number in index_table['profile_number']:
            profile = np.genfromtxt(f"../../../{folder}/mass-{mass}Msun_acchist-{acchist}/LOGS/profile{profile_number}.data",
                skip_header = 5, names = True)
            
            data = profile[f"{key}"]
            zone = profile['zone']
            sort = np.argsort(zone)

            np.savetxt(f"{key}/{key}-{profile_number}.txt", data[sort].T)
            print(f"{key}/{key}-{profile_number}.txt")
        
collect_profile_data(mass=mass, acchist=acchist)