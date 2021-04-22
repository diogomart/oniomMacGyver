#!/usr/bin/env python

# python modules
from omg import geom
import numpy as np
from openbabel import openbabel as ob

class PeriodicTable:
    table = {1: 'H', 2: 'He', 3: 'Li', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 11: 'Na', 12: 'Mg',
             15: 'P', 16: 'S', 17: 'Cl', 19: 'K', 20: 'Ca', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni',
             29: 'Cu', 30: 'Zn', 34: 'Se', 35: 'Br', 53: 'I'}
    def GetSymbol(self, atomic_number):
        return self.table[atomic_number]

PERIODIC_TABLE = PeriodicTable()


class Atom(ob.OBAtom):
    """A container for all the information relative to an atom

    Attributes:
        type (str)
        x (float)
        y (float)
        z (float)
        
    
    Other info is kept in its own object inside atom:
        - mm      (molecular mechanics)
        - oniom   (oniom info)
        - resinfo (information about the residue the atom belongs to)
        - pdbinfo (information for pdb files)
    """
    
    def __init__(self, atype, xyz):
        """
        
        Args:
            atype (one letter str or entire name): Element of this atom
            xyz (list): a list with the x, y, z coordinates
       
        Example:
            Atom('H', (1.0, 2.0, 3.0))
        """
        ob.OBAtom.__init__(self)
        self.SetType(atype) 
        self.SetVector(float(xyz[0]), float(xyz[1]), float(xyz[2]))
        self.mm = None
        self.oniom = None
        self.resinfo = None
        self.pdbinfo = None
        self.fragment = None
        
    def __repr__(self):
        return self.GetType()

    def set_mm(self, mm_obj):
        self.mm = mm_obj

    def set_oniom(self, oniom_obj):
        self.oniom = oniom_obj

    def set_resinfo(self, resinfo_obj):
        self.resinfo = resinfo_obj

    def set_pdbinfo(self, pdbinfo_obj):
        self.pdbinfo = pdbinfo_obj

    def GetDihedral(self, b, c, d):
        a = (self.GetX(), self.GetY(), self.GetZ())
        b = (b.GetX(), b.GetY(), b.GetZ())
        c = (c.GetX(), c.GetY(), c.GetZ())
        d = (d.GetX(), d.GetY(), d.GetZ())
        return geom.dihedral(a, b, c, d)

class MM(object):
    def __init__(self, atype, charge):
        self.atype = atype
        self.charge = float(charge)

class Oniom(object):
    def __init__(self, mask, layer):
        self.mask = int(mask)
        self.layer = layer
        self.has_link = False
        self.link_atom = None
        self.link_bound_to = None
        self.link_scale1 = None

    def set_link(self, atom, bound_to, scale1):
        self.link_atom = atom
        self.link_bound_to = int(bound_to)
        self.link_scale1 = float(scale1)
        self.has_link = True

class RESinfo(object):
    """ pdb info in g09 Oniom"""
    def __init__(self, name, resname, resnum, chain):
        self.name = name
        self.resname = resname
        self.resnum = int(resnum)
        self.chain = chain

class PDBinfo(object):
    def __init__(self, keyword, serial, 
        occupancy = .0, bfact = .0, altloc = "", icode = "",
        formalcharge = ""):
        self.keyword = keyword
        self.serial = serial
        self.occupancy      = occupancy
        self.bfact          = bfact
        self.altloc         = altloc
        self.icode          = icode
        self.formalcharge   = formalcharge # exists in pybel

