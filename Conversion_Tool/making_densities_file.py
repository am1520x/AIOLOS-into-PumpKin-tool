"""
Gets the number densities and Temperature at a specific radius for all species
"""
import numpy as np
import os

def get_essential_data(directory, sim, timestep, species, rplanet=-1):
    # Load first species to get dimensions
    data_h0 = np.loadtxt(f"{directory}output_{sim}_{species[0][0]}_t{timestep}.dat", skiprows=2, usecols=[0, 1, 2, 4, 10, 11, 12, 19], max_rows=410)
    max_rr = len(data_h0)
    
    if rplanet < 0:
        rp = data_h0[2, 0]  # Planet radius from data if not provided
    else:
        rp = rplanet

    # Load all species data
    all_species_datas = []
    for s, spc in enumerate(species):
        file = f"{directory}output_{sim}_{spc[0]}_t{timestep}.dat"
        if os.path.isfile(file):
            tmpdata = np.loadtxt(file, skiprows=2, usecols=[0, 1, 4, 10, 12], max_rows=max_rr-1)
        else:
            tmpdata = np.zeros_like(data_h0)
        all_species_datas.append(tmpdata)

    # Calculate averaged quantities
    r = all_species_datas[0][:, 0] / rp
    total_n = np.zeros_like(r)
    avg_T = np.zeros_like(r)
    
    for s, spc_data in enumerate(all_species_datas):
        total_n += spc_data[:, 2]  # Number density
        avg_T += spc_data[:, 2] * spc_data[:, 4]  # n * T
    
    avg_T /= total_n  # Temperature weighted by number density

    return avg_T, r, all_species_datas

def process_timesteps(directory, sim, timesteps, species, rplanet=1.37e8, index=10, output_file="species_densities.txt", output_file_2="qt_conditions.txt"):
    """
    Process multiple timesteps and write number densities at a specific index to a file.
    Returns a nested dictionary with the data and list of average temperatures.
    
    Args:
        directory: Path to simulation directory
        sim: Simulation name prefix
        timesteps: List of timestep numbers to process
        species: List of species specifications
        rplanet: Planet radius (default 1.37e8 cm)
        index: Radial index to extract data from (default 10)
        output_file: Name of output text file
        
    Returns:
        tuple: (avg_T_list, data_dict) where:
            - avg_T_list: List of average temperatures at the specified index
            - data_dict: Nested dictionary {'timestep': {'species': num_den}}
    """
    
    avg_T_list = []
    data_dict = {}
    header = "Time_s\t" #+ "\t".join([spc[0] for spc in species])
    with open(output_file, 'w') as f_densities, open(output_file_2, 'w') as f_temps:
        # Write headers
        f_densities.write(header + "\n")
        f_temps.write(header + "\n")
            
        for t in timesteps:
            try:
                # Get data for this timestep
                avg_T, r, all_species_datas = get_essential_data(
                    directory=directory,
                    sim=sim,
                    timestep=t,
                    species=species,
                    rplanet=rplanet
                )
                
                # Store avg_T at our index
                avg_T_list.append(avg_T[index])
                f_temps.write(f"{t}\t{avg_T[index]:.4f}\n")
                # Initialize dictionary entry for this timestep
                data_dict[t] = {}
                
                # Collect number densities at our index
                densities = []
                for s, spc_data in enumerate(all_species_datas):
                    num_den = spc_data[:, 2]  # Number density column
                    species_name = species[s][0]
                    densities.append(f"{num_den[index]:.4E}")
                    # Store in dictionary
                    data_dict[t][species_name] = num_den[index]
                
                # Write to file: timestep then all densities
                f_densities.write(f"{t}\t" + "\t".join(densities) + "\n")
                
            except Exception as e:
                print(f"Error processing timestep {t}: {str(e)}")
                # Add placeholder entry in dictionary for failed timesteps
                data_dict[t] = {spc[0]: np.nan for spc in species}
                continue
    return avg_T_list, data_dict
