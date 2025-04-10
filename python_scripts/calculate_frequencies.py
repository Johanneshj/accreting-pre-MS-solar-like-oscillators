### Script for generating GYRE inlists and calculating frequencies
### gyre6freqs.sh file adopted from Earl Bellinger's MESA/GYRE summer school 2022
### https://github.com/earlbellinger/mesa-gyre-tutorial-2022/

import numpy as np
import os
import pandas as pd
import subprocess

def write_bash_script(acchist, mass, profile_number, gyre6freqs_folder,
                      lower_val, upper_val, destination_folder):
        bash_script_content = f"""
        #!/bin/bash -i
        source ~/.bashrc
        gyre_file=profile{profile_number}.data.GYRE
        {gyre6freqs_folder}/gyre6freqs.sh -l {lower_val} -u {upper_val} -D {destination_folder} -i $gyre_file
        """
        with open(f"temp_script-mass_{mass}Msun-acchist_{acchist}.sh", "w") as f:
            f.write(bash_script_content)

def calculate_frequencies(acchist, mass, lower_bound=True, classII=True,
                          gyre6freqs_folder='',
                          path_to_data='',
                          path_to_destination_folder=''):
    
    full_path_to_data = path_to_data
    model = f"mass-{mass}Msun_acchist-{acchist}"
    full_path_to_model = os.path.join(full_path_to_data, model)
    
    # Load history files
    hist_model_number = np.genfromtxt(os.path.join(full_path_to_model, f"{model}_value_model_number.txt"))
    masses = np.genfromtxt(os.path.join(full_path_to_model, f"{model}_value_star_mass.txt"))
    center_h1 = np.genfromtxt(os.path.join(full_path_to_model, f"{model}_value_center_h1.txt"))
    numax = np.genfromtxt(os.path.join(full_path_to_model, f"{model}_value_nu_max.txt"))
    
    # Find indexes corresponding to class II objects
    if classII:
        mass_indexes = np.where(masses >= 0.9*mass)[0]
        center_h1_indexes = np.where( (center_h1 >= 0.95*max(center_h1)) )[0]
        classII_indexes = np.arange(mass_indexes[0], center_h1_indexes[-1])
        hist_model_number = hist_model_number[classII_indexes]
    
    # Load MESA index file
    profile_folder = os.path.join(full_path_to_model, "profile_data")
    index_file = np.genfromtxt(os.path.join(profile_folder, "index_file.txt"))
    index_model_number = np.array(index_file[:, 1], dtype=int)
    profile_numbers = np.array(index_file[:, 0], dtype=int)

    # Find profiles corresponding to class II
    
    idx = np.isin(index_model_number, hist_model_number)
    profile_numbers = profile_numbers[idx]
    indexes_for_numax = np.isin(hist_model_number, index_model_number[idx])

    # Numax values for frequency calculation
    numax = numax[indexes_for_numax]
    
    # Define lower and upper bounds and destination folder
    # for GYRE calculations/outputs
    for i, profile_number in enumerate(profile_numbers):
        best_numax = numax[i]
        
        if lower_bound:
            lower = best_numax + 2*(0.66*best_numax**0.88) 
        else:
            lower = 1
        upper = best_numax + 2*(0.66*best_numax**0.88)

        path_to_destination_folder = path_to_destination_folder
        destination_folder = os.path.join(path_to_destination_folder)     
        if not os.path.exists(destination_folder):
            os.mkdir(destination_folder)          
        
        # Write gyre inlist file
        write_bash_script(mass, acchist, profile_number=profile_number, 
                          gyre6freqs_folder=gyre6freqs_folder,
                          lower_val=lower, upper_val=upper,
                          destination_folder=destination_folder)
        
        # Run gyre inlist file
        subprocess.run(["bash", f"temp_script-mass_{mass}Msun-acchist_{acchist}.sh"], check=True)

# For parallelization
mass = float('<<mass>>')
acchist = int('<<acchist>>')

calculate_frequencies(acchist=acchist, mass=mass)
