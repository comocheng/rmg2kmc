# Class definitions for KMC writers

import os
import numpy as np
import shutil
import datetime
import abc
import cantera as ct


def rename_species(name):
    """function to rename species names to strings that can be used as Python variables"""
    name = name.replace('(', '_')
    name = name.replace(')', '')
    return name


def detect_reaction_type(reaction):
    """
    function to detect the forward reaction type as surface reaction adsorption or desorption
    convention is that the first reactant is the site, and the second is other_site
    """
    reaction_types = ['adsorption', 'desorption', 'surface_reaction']
    # adsorption/desorption has a gas phase species
    return 'adsorption'


class KMCWriter(abc.ABC):
    @abc.abstractmethod
    def write(self):
        raise NotImplementedError


class MonteCoffeeWriter(KMCWriter):
    def write(self, output_dir, species_list, reaction_list):
        # copy the template files to the output_dir
        template_dir = os.path.join(os.path.dirname(__file__), 'montecoffee', 'template_files')
        shutil.copy(os.path.join(template_dir, 'run.py'), output_dir)
        shutil.copy(os.path.join(template_dir, 'user_kmc.py'), output_dir)
        shutil.copy(os.path.join(template_dir, 'user_system.py'), output_dir)

        # Make a dictionary mapping species variable names to species

        # need to write the user_sites.py classes
        user_sites_path = os.path.join(output_dir, 'user_sites.py')
        user_sites_lines = self.generate_user_sites(species_list)
        with open(user_sites_path, 'w') as f:
            f.writelines(user_sites_lines)

        # need to write the user_events.py classes
        user_events_path = os.path.join(output_dir, 'user_events.py')
        user_events_lines = self.generate_user_events(reaction_list)
        with open(user_events_path, 'w') as f:
            f.writelines(user_events_lines)

    def generate_user_sites(self, species_list):
        lines = [f'# User sites autogenerate by rmg2kmc {datetime.datetime.now()}\n']
        lines.append('import base.sites\n\n\n')

        # for now, this only handles 4 site types for 111 facets
        lines.append('# Define constant site types\n')
        lines.append('SITE_FCC111_TOP = 0\n')
        lines.append('SITE_FCC111_BRIDGE = 1\n')
        lines.append('SITE_FCC111_FCC_HOLLOW = 2\n')
        lines.append('SITE_FCC111_HCP_HOLLOW = 3\n')
        lines.append('\n')
        lines.append('# Define species constants\n')

        # convert cantera names to python variable names
        # TODO check that no species was lost in the conversion
        # TODO reserve zero for vacant site
        for i, species in enumerate(species_list):
            lines.append(f'SPECIES_{rename_species(species.name)} = {i}\n')
        lines.append('\n')

        # add a basic site
        lines.append('class Site(base.sites.SiteBase):\n')
        lines.append('    def __init__(\n')
        lines.append('        self,\n')
        lines.append('        stype=SITE_FCC111_TOP,\n')
        lines.append('        covered=0,      # covered is the species index\n')
        lines.append('        ind=[],         # ase indices of site-related atoms\n')
        lines.append('        lattice_pos=None\n')
        lines.append('    ):\n')
        lines.append('        base.sites.SiteBase.__init__(\n')
        lines.append('            self,\n')
        lines.append('            stype=stype,\n')
        lines.append('            covered=covered,\n')
        lines.append('            ind=ind,\n')
        lines.append('            lattice_pos=lattice_pos\n')
        lines.append('        )\n')

        return lines

    def generate_user_events(self, reaction_list):
        # TODO remove dependence of MonteCoffee simulations on Cantera
        lines = [f'# User events autogenerate by rmg2kmc {datetime.datetime.now()}\n']
        lines.append('import cantera as ct\n')
        lines.append('import base.events\n')
        lines.append('import user_sites\n\n\n')

        # need to convert reaction into a class name
        # TODO add Diffusion reactions
        for i, reaction in enumerate(reaction_list):
            if len(reaction.reactants) > 2:
                raise NotImplementedError('Trimolecular reactions not yet supported')

            # figure out if it's a surface reaction or adsorption
            reaction_type = detect_reaction_type(reaction)

            # define forward reaction
            lines.append(f'class Reaction{i}Fwd(base.events.EventBase):\n')
            lines.append(f'    # {reaction.equation}\n')
            lines.append('    def __init__(self, params):\n')
            lines.append('        base.events.EventBase.__init__(self, params, name={reaction.equation})\n\n')

            lines.append('    def possible(self, system, site, other_site):\n')
            # TODO actually fill this out
            lines.append('        return True\n\n')

            # Define the reaction rate
            lines.append('    def get_rate(self, system, site, other_site):\n')
            # TODO add temperature dependence
            if not isinstance(reaction.rate, ct.Arrhenius):
                raise NotImplementedError('Kinetics only implemented for type Arrhenius')
            lines.append('        T = 1000.0\n')
            A = reaction.rate.pre_exponential_factor
            Ea = reaction.rate.activation_energy
            b = reaction.rate.temperature_exponent
            lines.append(f'        kinetics = ct.Arrhenius(A={A}, b={b}, E={Ea})\n')
            lines.append('        return kinetics(T)\n\n')

            # TODO fill out actual event actions
            lines.append('    def do_event(self, system, site, other_site):\n')
            lines.append('        pass\n\n')

            lines.append('    def get_involve_other(self):\n')
            lines.append('        return False\n\n\n')

            if reaction.reversible:
                # TODO define reverse reaction
                pass

        return lines



def get_species_name(species, species_list):
    for sp in species_list:
        if sp.is_isomorphic(species):
            return sp.label

class ZacrosWriter(KMCWriter):
    def write_lattice_file(self, output_dir):
        
        lattice_path = os.path.join(output_dir, 'lattice_input.dat')
        lines = []
        lines.append(f'# Lattice input autogenerate by rmg2kmc {datetime.datetime.now()}\n')
        lines.append('lattice  periodic_cell\n\n')
        lines.append('# R = sqrt(2) a / 4\n')
        lines.append('# a = 3.9239A\n')
        lines.append('# R = 1.3873 A\n')
        lines.append('# (2R, 0)\n')
        lines.append('# (R, sqrt(3)R)\n')

        lines.append('cell_vectors       # in row format (Angstroms)\n')
        lines.append('   2.774616298697894   0.000000000000000\n')
        lines.append('   1.387308149348947   2.402888200426728\n\n')

        lines.append('repeat_cell       20   20\n\n')

        lines.append('n_site_types      1\n')
        lines.append('site_type_names   top\n\n')

        lines.append('n_cell_sites      1\n')
        lines.append('site_types        top\n\n')

        lines.append('site_coordinates   # fractional coordinates (x,y) in row format\n')
        lines.append('   0.000000000000000   0.000000000000000\n\n')
        
        lines.append('neighboring_structure\n')
        lines.append('   1-1  north\n')
        lines.append('   1-1  east\n')
        lines.append('   1-1  southeast\n\n')
        lines.append('end_neighboring_structure\n\n')

        lines.append('end_lattice\n\n')

        with open(lattice_path, 'w') as f:
            f.writelines(lines)

    def write_energetics_file(self, output_dir, species_list):
        energetics_path = os.path.join(output_dir, 'energetics_input.dat')
        lines = []
        lines.append(f'# Energetics input autogenerate by rmg2kmc {datetime.datetime.now()}\n')
        lines.append('energetics\n\n')

        for species in species_list:
            if not species.contains_surface_site():
                continue
            if species.is_surface_site():
                continue
            lines.append(f'cluster {species.label}_top\n')  # TODO other sites?
            lines.append('  sites 1\n')
            lines.append('  lattice_state\n')
            lines.append(f'    1 {species.label}   1\n')
            lines.append('  site_types top\n')
            lines.append('  cluster_eng 0.0\n\n')
            lines.append(f'end_cluster {species.label}_top\n')  # TODO other sites?

        lines.append('end_energetics\n\n')
        with open(energetics_path, 'w') as f:
            f.writelines(lines)


    def write(self, output_dir, species_list, reaction_list, T):
        # make the mechanism file
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

        # WRITE THE GENERAL SIMULATION INPUT FILE
        simulation_path = os.path.join(output_dir, 'simulation_input.dat')
        lines = []
        lines.append(f'# Simulation input file autogenerate by rmg2kmc {datetime.datetime.now()}\n')
        random_seed = 8949321
        lines.append(f'random_seed\t\t{random_seed}\n\n')
        lines.append(f'temperature\t\t{float(T)}  # K\n')
        lines.append(f'pressure\t\t1.01325  # bar\n\n')

        # get num gas species
        num_gas_species = 0
        num_surf_species = 0
        gas_species = []
        surf_species = []

        for species in species_list:
            if species.contains_surface_site():
                surf_species.append(species)
                num_surf_species += 1
            else:
                gas_species.append(species)
                num_gas_species += 1
        

        lines.append(f'n_gas_species\t\t{num_gas_species}\n')
        lines.append(f'gas_specs_names\t\t{" ".join([str(sp.label) for sp in gas_species])}\n')
        lines.append(f'gas_energies\t\t{" ".join([str(0.00) for sp in gas_species])}  # eV\n')
        lines.append(f'gas_molec_weights\t\t{" ".join([str(np.round(sp.molecular_weight.value, 3)) for sp in gas_species])}\n')
        
        # TODO add something that makes sense
        lines.append(f'gas_molar_fracs\t\t{" ".join([str(np.round(1.0 / num_gas_species, 2)) for sp in gas_species])}\n\n')

        lines.append(f'n_surf_species\t\t{num_surf_species}\n')
        lines.append(f'surf_specs_names\t\t{" ".join([str(sp.label) for sp in surf_species])}\n')
        # TODO add bidentate species
        lines.append(f'surf_specs_dent\t\t{" ".join(["1" for sp in surf_species])}\n\n')

        # TODO get timeframe from Cantera simulation?
        lines.append(f'snapshots\t\ton time 1e-8\n')
        lines.append(f'process_statistics\t\ton time 1e-8\n')
        lines.append(f'species_numbers\t\ton time 1e-8\n\n')

        lines.append(f'# event_report\t\ton\n\n')

        lines.append(f'max_steps\t\tinfinity  # seconds\n')
        lines.append(f'max_time\t\t5e-6  # seconds\n\n')

        lines.append(f'wall_time\t\t3000  # seconds\n\n')

        lines.append(f'no_restart\t\t\n\n')

        lines.append(f'# debug_report_processes\n')
        lines.append(f'# debug_report_global_energetics\n')
        lines.append(f'# debug_check_processes\n')
        lines.append(f'# debug_check_lattice\n\n')

        lines.append(f'finish\n\n')



        with open(simulation_path, 'w') as f:
            f.writelines(lines)


        # WRITE THE MECHANISM FILE
        mechanism_path = os.path.join(output_dir, 'mechanism_input.dat')
        lines = []
        lines.append(f'# Mechanism file autogenerate by rmg2kmc {datetime.datetime.now()}\n')
        lines.append(f'# Temperature: {T}K\n')
        lines.append('\n')
        lines.append('\n')
        lines.append('mechanism\n')
        lines.append('\n')

        # write each reaction as an irreversible step
        for reaction in reaction_list:

            step_name = ''.join(str(reaction).split())
            lines.append(f'step {step_name}\n\n')

            gas_reacs_prods_string = ''
            for sp in reaction.reactants + reaction.products:
                if sp.contains_surface_site():
                    continue

                gas_reacs_prods_string += f' {get_species_name(sp, species_list)} {reaction.get_stoichiometric_coefficient(sp)}'

            lines.append(f'  gas_reacs_prods {gas_reacs_prods_string}\n')
            # TODO read the sites from the adjacency list

            # TODO add more site types
            lines.append('  sites 1\n')
            dentate = 1
            if any([sp.is_surface_site() for sp in reaction.reactants]):
                lines.append('  initial\n')
                lines.append(f'    1  *  {dentate}\n')
                lines.append('  final\n')
                # This needs to be the new surface species
                lines.append(f'    1  {reaction.products[0].label}  {dentate}\n')
            else:
                lines.append('  initial\n')
                lines.append(f'    1  {reaction.reactants[0].label}  {dentate}\n')
                lines.append('  final\n')
                lines.append(f'    1  *  {dentate}\n')
            lines.append('\n')

            # TODO add more site types
            lines.append('  variant top\n')
            lines.append('    site_types\ttop\n')
            lines.append(f'    pre_expon\t{reaction.kinetics.A.value_si}\n')
            lines.append(f'    activ_eng\t0.00  # eV\n')
            lines.append('  end_variant\n\n')
            lines.append('end_step\n\n\n')
        
        lines.append('end_mechanism\n')

        with open(mechanism_path, 'w') as f:
            f.writelines(lines)


        # WRITE THE LATTICE FILE
        self.write_lattice_file(output_dir)

        # WRITE THE ENERGETICS FILE
        self.write_energetics_file(output_dir, species_list)
