from calculating_rates import parse_reactions, calc_rate_coeff, calculate_reaction_rates 

def make_rates(output_file, avg_T_at_index, number_densities, timesteps):
    # 1. Parse the reactions from file
    reactions = parse_reactions('test.reac')
    header = "Time_s\t" #+ "\t".join([spc[0] for spc in species])   
    with open(output_file, 'w') as f:
        # Write header
        f.write(header + "\n")

        for i, t in enumerate(timesteps):
            rates = calculate_reaction_rates(reactions, number_densities[t], calc_rate_coeff, avg_T_at_index[i])
            # Write to file: timestep then all densities
            f.write(f"{t}\t" + "\t".join(rates) + "\n")
                
       
    