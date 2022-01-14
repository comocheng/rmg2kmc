# Class definitions for KMC writers

import os
import shutil
import datetime
import abc
import cantera as ct


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

        # need to write the user_sites.py classes
        user_sites_path = os.path.join(output_dir, 'user_sites.py')
        user_sites_lines = self.generate_user_sites(species_list)
        with open(user_sites_path, 'w') as f:
            f.writelines(user_sites_lines)

        # need to write the user_events.py classes
        user_events_path = os.path.join(output_dir, 'user_events.py')
        user_events_lines = self.generate_user_sites(reaction_list)
        with open(user_events_path, 'w') as f:
            f.writelines(user_events_lines)

    def generate_user_sites(self, species_list):
        lines = [f'User sites autogenerate by rmg2kmc {datetime.datetime.now()}\n']
        return lines

    def generate_user_events(self, reaction_list):
        lines = [f'User events autogenerate by rmg2kmc {datetime.datetime.now()}\n']
        return lines


class ZacrosWriter(KMCWriter):
    def write(self, output_dir, species_list, reaction_list):
        pass
