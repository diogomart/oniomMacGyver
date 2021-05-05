#!/usr/bin/env python

import sys
import argparse
from omg.gaussian.gaussian import GaussianLog
from rdkit import Chem
from rdkit.Geometry import Point3D

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--template_sdf', required=True, help="same atom order as in .log file")
    parser.add_argument('-l', '--log', help='Gaussian .log file', required=True)
    parser.add_argument('-o', '--output_fname', help='output .sdf file')
    return parser.parse_args()

if __name__ == "__main__":
    periodic_table = Chem.GetPeriodicTable()
    args = get_args()
    gaulog = GaussianLog(args.log)
    byte = gaulog.grep_bytes['orientation:'][-1][-1]
    coords = gaulog.read_coordinates('all', byte)
    if not (args.template_sdf.endswith('.sdf') or args.template_sdf.endswith('.mol')):
        print('template_sdf must end with SDF or MOL', file=sys.stderr)
        sys.exit(2)
    mol = Chem.MolFromMolFile(args.template_sdf, removeHs=False)
    atoms = mol.GetAtoms()
    assert(len(coords) == len(atoms))
    conformer = mol.GetConformer()
    i = 0
    for (x, y, z) in coords:
        conformer.SetAtomPosition(i, Point3D(x, y, z))
        element_gau = gaulog.atoms_list[i].GetType() # SetType(element) in oniomMacGyver/atoms.py
        element_sdf = periodic_table.GetElementSymbol(atoms[i].GetAtomicNum()) 
        assert(element_gau == element_sdf)
        i += 1

    output_string = Chem.MolToMolBlock(mol)
    if args.output_fname:
        print(output_string, file=open(args.output_fname, 'w'))
    else:
        print(output_string)
