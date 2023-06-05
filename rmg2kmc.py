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

else:
    raise ValueError("Mechanism file extension not recognized")


###########################################################
# Write
if output_format == 'montecoffee':
    writer = kmc_writer.MonteCoffeeWriter()
    output_dir = os.path.dirname(os.path.join(mechanism_file, 'montecoffee'))
elif output_format == 'zacros':
    writer = kmc_writer.ZacrosWriter()
    output_dir = os.path.dirname(os.path.join(mechanism_file, 'zacros'))
else:
    raise ValueError(f'Output format {output_format} not recognized')

writer.write(output_dir, species_list, reaction_list)
