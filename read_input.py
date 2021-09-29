# script to read in the mechanism file
from rmgpy.chemkin import load_chemkin_file


def import_mechanism():
    print("reading mechanism")


my_cti = '/home/moon/methanol/perturb_5000/run_0000/cantera/chem_annotated.cti'
my_chem = '/home/moon/methanol/perturb_5000/run_0000/chemkin/chem_annotated-surface.inp'
my_dict = '/home/moon/methanol/perturb_5000/run_0000/chemkin/species_dictionary.txt'

species_list, reaction_list = load_chemkin_file(my_chem, dictionary_path=my_dict)
print(species_list)


output_dir = "/home/moon/autokmc/task1/"
