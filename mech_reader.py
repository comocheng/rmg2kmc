# Class to read in mechanism files
# from rmgpy.chemkin import load_chemkin_file


# def import_mechanism():
#     print("reading mechanism")


# my_cti = '/home/moon/methanol/perturb_5000/run_0000/cantera/chem_annotated.cti'
# my_chem = '/home/moon/methanol/perturb_5000/run_0000/chemkin/chem_annotated-surface.inp'
# my_dict = '/home/moon/methanol/perturb_5000/run_0000/chemkin/species_dictionary.txt'

# species_list, reaction_list = load_chemkin_file(my_chem, dictionary_path=my_dict)
# print(species_list)


# output_dir = "/home/moon/autokmc/task1/"


import abc


class MechanismReader(abc.ABC):
    """
    abstract class to read in mechanism files and create a list of species and reactions
    Takes in mechanism files
    """
    @abc.abstractmethod
    def read(self):
        raise NotImplementedError


class CanteraMechanismReader(MechanismReader):
    """
    Read in the Cantera mechanism
    """
    def read(self, mech_file, species_dict=None):
        import cantera as ct

        gas = ct.Solution(mech_file, "gas")
        surf = ct.Interface(mech_file, "surface1", [gas])

        species_list = surf.species()
        reaction_list = surf.reactions()
        # print(surf.species_list)
        return species_list, reaction_list


class ChemkinMechanismReader(MechanismReader):
    def read(self):
        return super().read()
