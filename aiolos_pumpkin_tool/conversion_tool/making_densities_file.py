"""
Module for extracting number densities and temperature profiles of species \nfrom AIOLOS outputs.

This script supports:
- Loading and aggregating data across species.
- Extracting data at specific radial indices or over all radial bins.
- Writing outputs to file for later analysis or plotting.

Functions:
- get_essential_data: Load and average simulation data for selected species.
- process_timesteps: Process multiple timesteps and record values at a specific radius.
- process_radial_profile: Extract full radial profile of number densities and temperature.
"""

import numpy as np
import os


def get_essential_data(directory, sim, timestep, species, rplanet=-1):
    """
    Load number density and temperature data for a set of species at a specific timestep.

    This function reads AIOLOS simulation output files for specified species and
    extracts key data including number densities, temperatures, and radial coordinates.
    It handles both individual species names and species data arrays.

    Args:
        directory (str): Path to directory containing output files. Should end with '/'.
        sim (str): Simulation name prefix used in filenames.
        timestep (int): Simulation timestep number.
        species (list): List of species names or species data arrays to load.
            Each element can be either a string (species name) or a list where
            the first element [0] is the species name.
        rplanet (float, optional): Planet radius in cm for normalizing radius column.
            If negative, auto-detects from the first data file. Defaults to -1.

    Returns:
        tuple: A 3-element tuple containing:
            - avg_T (np.ndarray): Array of mass-weighted average temperatures at each radial bin (K).
            - r (np.ndarray): Array of normalized radius values (dimensionless, r/rplanet).
            - all_species_datas (list): List of numpy arrays, one per species, with shape (n_radial, 5).
                Each array contains columns: [radius, ?, number_density, ?, temperature].

    Raises:
        FileNotFoundError: If the first species file cannot be found.
        ValueError: If data files have inconsistent formats or invalid data.

    Examples:
        >>> # Load data for simple species list
        >>> avg_T, r, datas = get_essential_data("../data/", "sim1", 10, ["H", "He"])
        >>> len(datas)  # Should equal number of species
        2
        >>> avg_T.shape == r.shape  # Both arrays should have same length
        True

        >>> # Load data with species data arrays (typical format from main pipeline)
        >>> species_arrays = [["H", 1, 1, "$\\\\rm H$"], ["He", 4, 2, "$\\\\rm He$"]]
        >>> avg_T, r, datas = get_essential_data("../data/", "sim1", 10, species_arrays)
        >>> isinstance(avg_T, np.ndarray) and len(avg_T) > 0
        True
    """
    # Load first species to get dimensions
    first_species_name = species[0][0] if isinstance(species[0], list) else species[0]
    first_file_path = f"{directory}output_{sim}_{first_species_name}_t{timestep}.dat"

    data_h0 = np.genfromtxt(
        first_file_path,
        skip_header=2,
        usecols=[0, 1, 2, 4, 10, 11, 12, 19],
        max_rows=410,
        invalid_raise=False,
        filling_values=np.nan,  # Replace invalid entries like '-' with NaN
    )

    max_rr = data_h0.shape[0]

    if rplanet < 0:
        rp = data_h0[2, 0]  # Planet radius from data if not provided
    else:
        rp = rplanet

    # Load all species data
    all_species_datas = []
    for s, spc in enumerate(species):
        species_name = spc[0] if isinstance(spc, list) else spc
        file = f"{directory}output_{sim}_{species_name}_t{timestep}.dat"
        if os.path.isfile(file):
            tmpdata = np.genfromtxt(
                file,
                skip_header=2,
                usecols=[0, 1, 4, 10, 12],
                max_rows=max_rr,
                invalid_raise=False,
                filling_values=np.nan,
            )
        else:
            tmpdata = np.zeros_like(data_h0)
        all_species_datas.append(tmpdata)

    # Calculate averaged quantities
    r = all_species_datas[0][:, 0] / rp
    total_n = np.zeros_like(r)
    avg_T = np.zeros_like(r)

    for spc_data in all_species_datas:
        # Replace NaNs with 0 before using data
        n_vals = np.nan_to_num(spc_data[:, 2])
        t_vals = np.nan_to_num(spc_data[:, 4])
        total_n += n_vals
        avg_T += n_vals * t_vals

    # Avoid division by zero
    with np.errstate(divide="ignore", invalid="ignore"):
        avg_T = np.where(total_n > 0, avg_T / total_n, 0)
    return avg_T, r, all_species_datas


def process_timesteps(
    directory,
    sim,
    timesteps,
    species,
    rplanet=1.37e8,
    index=10,
    output_file="qt_densities.txt",
    output_file_2="qt_conditions.txt",
):
    """
    Process multiple timesteps and extract number densities at a specific radial location.

    This function processes a time series of AIOLOS simulation outputs, extracting
    number densities and temperatures at a fixed radial index across multiple timesteps.
    Results are written to output files and returned as structured data.

    Args:
        directory (str): Path to simulation output directory. Should end with '/'.
        sim (str): Simulation name prefix used in output filenames.
        timesteps (list[int]): List of timestep indices to process sequentially.
        species (list): List of species names or species data arrays.
            Each element can be either a string or a list where [0] is the species name.
        rplanet (float, optional): Planet radius in cm for normalization. Defaults to 1.37e8 cm.
        index (int, optional): Radial bin index to extract data from. Defaults to 10.
        output_file (str, optional): Output filename for species densities. Defaults to "qt_densities.txt".
        output_file_2 (str, optional): Output filename for temperature data. Defaults to "qt_conditions.txt".

    Returns:
        tuple: A 2-element tuple containing:
            - avg_T_list (list[float]): Average temperatures at specified radial index for each timestep (K).
            - data_dict (dict): Nested dictionary with structure {timestep: {species_name: number_density}}.
                Failed timesteps will have NaN values for all species.

    Raises:
        IOError: If output files cannot be written.
        Exception: Individual timestep failures are caught and logged, with NaN placeholders added.

    Examples:
        >>> # Process time evolution at radial index 10
        >>> temps, data = process_timesteps("../data/", "sim1", [0, 5, 10, 15], ["H", "He"])
        >>> len(temps) == 4  # Should match number of timesteps
        True
        >>> isinstance(data[0], dict)  # Each timestep should be a dictionary
        True
        >>> "H" in data[0]  # Species should be present in data
        True

        >>> # Check data structure
        >>> isinstance(data[0]["H"], (int, float))  # Number density should be numeric
        True
    """

    avg_T_list = []
    data_dict = {}
    header = "Time_s\t"
    with open(output_file, "w") as f_densities, open(output_file_2, "w") as f_temps:
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
                    rplanet=rplanet,
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


def process_radial_profile(
    directory,
    sim,
    timestep,
    species,
    rplanet=1.37e8,
    output_file="qt_densities.txt",
    output_file_2="qt_conditions.txt",
):
    """
    Extract full radial profiles of number densities and temperatures at a single timestep.

    This function processes a single timestep from AIOLOS simulation outputs to extract
    the complete radial distribution of species number densities and temperatures across
    all radial bins from the planet surface to the simulation boundary.

    Args:
        directory (str): Path to simulation output directory. Should end with '/'.
        sim (str): Simulation name prefix used in output filenames.
        timestep (int): Specific timestep index to process.
        species (list): List of species names or species data arrays.
            Each element can be either a string or a list where [0] is the species name.
        rplanet (float, optional): Planet radius in cm for normalization. Defaults to 1.37e8 cm.
        output_file (str, optional): Output filename for radial number density profiles.
            Defaults to "qt_densities.txt".
        output_file_2 (str, optional): Output filename for radial temperature profile.
            Defaults to "qt_conditions.txt".

    Returns:
        tuple: A 2-element tuple containing:
            - avg_T_list (list[float]): Average temperatures at each radial bin (K).
                Length equals the number of radial bins in the simulation.
            - data_dict (dict): Nested dictionary with structure {radius_index: {species_name: number_density}}.
                Each radius_index corresponds to a radial bin from 0 (planet surface) to max_radius.

    Raises:
        FileNotFoundError: If simulation output files for the specified timestep don't exist.
        IOError: If output files cannot be written.
        ValueError: If data files have inconsistent formats.

    Examples:
        >>> # Extract radial profile at timestep 85
        >>> temps, data = process_radial_profile("../data/", "sim1", 85, ["H", "He", "H2"])
        >>> len(temps) > 100  # Typical simulation has >100 radial bins
        True
        >>> len(data) == len(temps)  # Should have data for each radial bin
        True
        >>> isinstance(data[0]["H"], (int, float))  # Number densities should be numeric
        True

        >>> # Check radial structure
        >>> all(isinstance(idx, int) for idx in data.keys())  # Keys should be integer indices
        True
        >>> "H" in data[10] and "He" in data[10]  # All species should be present
        True
    """
    avg_T_list = []
    data_dict = {}

    header = "Cell_Number\t"

    with open(output_file, "w") as f_densities, open(output_file_2, "w") as f_temps:
        # Write headers
        f_densities.write(header + "\n")
        f_temps.write("Radius_index\tAvg_Temperature\n")

        # Load all data
        avg_T, r, all_species_datas = get_essential_data(
            directory=directory,
            sim=sim,
            timestep=timestep,
            species=species,
            rplanet=rplanet,
        )

        n_radial = len(avg_T)

        for i in range(n_radial):
            data_dict[i] = {}
            densities = []

            for s, spc_data in enumerate(all_species_datas):
                num_den = spc_data[i, 2]
                species_name = species[s][0]
                densities.append(f"{num_den:.4E}")
                data_dict[i][species_name] = num_den

            avg_T_list.append(avg_T[i])
            f_temps.write(f"{i}\t{avg_T[i]:.4f}\n")
            f_densities.write(f"{i}\t" + "\t".join(densities) + "\n")

    return avg_T_list, data_dict
