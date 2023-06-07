# main script to convert mechanism file to KMC input
import os
import sys
import mech_reader
import kmc_writer


output_format = 'zacros'
# output_format = 'montecoffee'

###########################################################
# Read

# Check for correct usage. Must specify input file.
if len(sys.argv) < 2:
    ValueError('Need to specify mechanism input file')

mechanism_file = sys.argv[1]
if mechanism_file.endswith('.cti') or mechanism_file.endswith('.yaml'):
    print("Using Cantera Reader")
    reader = mech_reader.CanteraMechanismReader()
    species_list, reaction_list = reader.read(mechanism_file)

elif mechanism_file.endswith('.inp'):
    print("Using Chemkin Reader")
    
    gas_mech_file = sys.argv[1]
    surface_mech_file = sys.argv[2]
    if len(sys.argv) > 3:
        species_dict_file = sys.argv[3]
    else:
        species_dict_file = None

    reader = mech_reader.ChemkinMechanismReader()
    species_list, reaction_list = reader.read(gas_mech_file, surface_mech_file, species_dict_file)

    # get the site density from the file
    with open(surface_mech_file, 'r') as f:
        for line in f:
            if 'SDEN' in line and ('mol/cm2' in line or 'mol/cm^2' in line):
                site_density = float(line.split('/')[1]) * 100.0 * 100.0  # convert to mol/m2
                break
        else:
            site_density = 2.72E-5  # default value
            raise ValueError('Could not find site density in surface mechanism file')
        
    # specify the initial concentrations
    starting_gas_conc = ''
    for sp in species_list:
        if sp.contains_surface_site():
            continue
        elif sp.label == 'Ar':
            starting_gas_conc += ' 0.9'
        elif sp.label == 'CO':
            starting_gas_conc += ' 0.1'
        else:
            starting_gas_conc += ' 0.0'

else:
    raise ValueError("Mechanism file extension not recognized")


###########################################################
# Write
if output_format == 'montecoffee':
    writer = kmc_writer.MonteCoffeeWriter()
    output_dir = os.path.join(os.path.dirname(mechanism_file), 'montecoffee')
elif output_format == 'zacros':
    writer = kmc_writer.ZacrosWriter()
    output_dir = os.path.join(os.path.dirname(mechanism_file), 'zacros')
else:
    raise ValueError(f'Output format {output_format} not recognized')



writer.write(output_dir, species_list, reaction_list, starting_gas_conc, T=1000, site_density=site_density)
