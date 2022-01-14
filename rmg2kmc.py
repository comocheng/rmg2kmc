# main script to convert mechanism file to KMC input
import os
import sys
import mech_reader
import kmc_writer


# output_format = 'zacros'
output_format = 'montecoffee'

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
    raise NotImplementedError()

else:
    raise ValueError("Mechanism file extension not recognized")


###########################################################
# Write
if output_format == 'montecoffee':
    writer = kmc_writer.MonteCoffeeWriter()
elif output_format == 'zacros':
    writer = kmc_writer.ZacrosWriter()
else:
    raise ValueError(f'Output format {output_format} not recognized')

output_dir = os.path.dirname(mechanism_file)
writer.write(output_dir, species_list, reaction_list)
