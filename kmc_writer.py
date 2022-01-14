# Class definitions for KMC writers

import os
import shutil
import datetime
import abc
import cantera as ct


def rename_species(name):
    """function to rename species names to strings that can be used as Python variables"""
    name = name.replace('(', '_')
    name = name.replace(')', '')
    return name


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


class ZacrosWriter(KMCWriter):
    def write(self, output_dir, species_list, reaction_list):
        pass
