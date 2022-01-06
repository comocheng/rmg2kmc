import numpy as np
from ase.build import fcc111


# sparse tensor to represent the system's atoms
# How to populate the list of reactions?
# start with sites
# Use ase to populate the system?


class Position:
    def __init__(self, x, y, z):
        self.position = [x, y, z]


class Atom:
    def __init__(self, pos=Position(0, 0, 0), element='Cu'):
        self.position = pos
        self.element = 'Cu'


class Lattice:
    def __init__(self):
        self.N = 10
        self.system_grid = np.zeros((self.N, self.N, self.N))
        self.vertices = {}
        p1 = Position(0, 0, 0)
        a1 = Atom(p1, 'Cu')
        p2 = Position(0, 1, 0)
        a2 = Atom(p2, 'Cu')
        p3 = Position(1, 0, 0)
        a3 = Atom(p3, 'Cu')
        p4 = Position(1, 1, 0)
        a4 = Atom(p4, 'Cu')

        self.vertices[p1] = a1
        self.vertices[p2] = a2
        self.vertices[p3] = a3
        self.vertices[p4] = a4

        self.edges = {}
    grid_spacing = 1.0

    atom_list = []
    # use a dictionary to represent a graph- one for edges and one for vertices
    # need a graph for the actual atoms


if __name__ == '__main__':
    lattice = Lattice()
    A = Atom()
    B = Atom()
    C = Atom()
    lattice.atom_list.append(A)
    lattice.atom_list.append(B)
    lattice.atom_list.append(C)
    print(lattice)

    # slab = fcc111('Cu', size=(10, 10, 1), vacuum=10.0)
