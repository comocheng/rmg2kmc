# abstract class definitions for converter
# from abc import ABC, abstractmethod
import abc


class KMCWriter(abc.ABC):
    @abc.abstractmethod
    def write(self):
        raise NotImplementedError


class MonteCoffeeWriter(KMCWriter):
    def write(self, output_dir, species_list, reaction_list):
        pass


class ZacrosWriter(KMCWriter):
    def write(self, output_dir, species_list, reaction_list):
        pass
