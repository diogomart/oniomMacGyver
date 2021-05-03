#!/usr/bin/env python

import os
import sys
import argparse
import subprocess
from omg.gaussian.gaussian import GaussianLog as GAULOG
from omg.gaussian.gaussian import GaussianCom as GAUCOM
from omg import atoms
from omg import iolines
from omg import geom
from openbabel import openbabel as ob

def pts_to_int(pts):
    try:
        pts = int(pts) - 1
        return int(pts)
    except:
        if pts =='all':
            return pts
        else:
            usage()
            sys.stderr.write('bad input:\n')
            sys.stderr.write('  expected int or all, but got %s:\n' % pts)
            sys.exit(2)

def get_args():
    "Parse arguments of gau_makepdb"
    parser = argparse.ArgumentParser(
        description="""
    Extract structures from gaussian .log or .com files
    and make .pdb files for visualization in VMD or Pymol.

    Typical usage:
        gau_makesdf.py scan_1.log scan_2.log 


""", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-t', '--template_sdf', required=True, help="same atom order as in .log files")
    parser.add_argument('-l', '--logs', help='Gaussian .log files', nargs='+')
    parser.add_argument('-m', '--monomers', help='Gaussian .log files', nargs='*', default=[])
    parser.add_argument('-s', '--scan_step',
                        help='Scan step ("all" or integer) default="all"',
                        default='all')
    parser.add_argument('-o', '--opt_step',
                    help='Scan step ("all" or integer) default = last = 0"',
                        default='0')
    parser.add_argument('--out', help='output .sdf file', required=True)
    args = parser.parse_args()

    # parse input extension
    for log in args.logs:
        basename, input_extension = os.path.splitext(log)
        if input_extension not in ['.log']:
            sys.stderr.write('Aborting!\n')
            sys.stderr.write('  inputname suffix is: %s\n' % input_extension)
            sys.stderr.write('  must be ".com" or ".log"\n')
            sys.exit(2)

    # scan and opt steps
    args.scan_step = pts_to_int(args.scan_step)
    args.opt_step  = pts_to_int(args.opt_step)

    return args

 
def extract_coordinates(gaulog, scan_pts, opt_pts, args):

    # get the bytes
    chosen_bytes = []
    if scan_pts == 'all':
        for sbytes in gaulog.grep_bytes['orientation:']:
            if opt_pts == 'all':
                for byte in sbytes:
                    chosen_bytes.append(byte) # this is crazy
            else:
                chosen_bytes.append(sbytes[opt_pts])
    else:
        sbytes = gaulog.grep_bytes['orientation:'][scan_pts]
        if opt_pts == 'all':
            for byte in sbytes:
                chosen_bytes.append(byte) # this is reasonable 
        else:
            chosen_bytes.append(sbytes[opt_pts])

    model_number = 0
    xyz_list = []
    for byte in chosen_bytes:
        model_number += 1
        xyz = gaulog.read_coordinates('all', byte)
        counter = 0 
        xyz_list.append(xyz)

    return xyz_list


def cat_modreds(GaussianLogList):
    """get studied distances, angles and dihedrals from GaussianLogList"""

    atomids = set() # tuples of 1-indexed atom indexes 
    for gl in GaussianLogList:
        if gl.modreds:
            for modred in gl.modreds:
                atomids.add(tuple(modred.atomids))
    return atomids


def get_level_of_theory(route_section):
    route_section = route_section.lower()
    isCounterpoise = False
    if 'counterpoise' in route_section:
        isCounterpoise = True
    if 'b3lyp/cc-pvdz' in route_section:
        level_of_theory = 'B3LYP_CC-PVDZ'
    elif 'b3lyp/6-31g*' in route_section:
        level_of_theory = 'B3LYP_6-31G*'
    elif 'b3lyp/aug-cc-pvdz' in route_section:
        level_of_theory = 'B3LYP_AUG-CC-PVDZ'
    elif 'wb97xd/aug-cc-pvdz' in route_section:
        level_of_theory = 'WB97XD_AUG-CC-PVDZ'
    elif 'b3lyp/def2svp' in route_section:
        level_of_theory = 'B3LYP_DEF2SVP'
    else:
        sys.stderr.write('Please add that level of theory\n')
        raise RuntimeError
    if isCounterpoise:
        level_of_theory += '_CP'
    return level_of_theory


if __name__ == '__main__':

    args = get_args()

    # build OBMol that openbabel will use to write output .sdf
    obconv = ob.OBConversion()
    obconv.SetInFormat('sdf')
    obmol = ob.OBMol()
    obconv.ReadFile(obmol, args.template_sdf)

    # get smiles of this molecule to write to the output .sdf
    obconv = ob.OBConversion()
    obconv.SetOutFormat('smi')
    smiles = obconv.WriteString(obmol).split()[0]

    # verify that all logs use the same level of theory
    level_of_theory = set()
    gaulogs = [GAULOG(log_filename) for log_filename in args.logs]
    for gaulog in gaulogs:
        level_of_theory.add(get_level_of_theory(gaulog.route_section))
    if len(level_of_theory) != 1: raise RuntimeError
    level_of_theory = level_of_theory.pop() # becomes string
    if len(args.monomers):
        monomer_level_of_theory = set()
        monomer_gaulogs = [GAULOG(log_filename) for log_filename in args.monomers]
        for gaulog in monomer_gaulogs:
            monomer_level_of_theory.add(get_level_of_theory(gaulog.route_section))
        if len(monomer_level_of_theory) != 1: raise RuntimeError
        monomer_level_of_theory = monomer_level_of_theory.pop() # becomes string
        if level_of_theory.replace('_CP', '') != monomer_level_of_theory:
            raise RuntimeError
        energy_apo = 0
        for gaulog in monomer_gaulogs:
            if len(gaulog.energies[energy_keyword_monomer]) != 1:
                raise RuntimeError('monomer does not seem to be an optimization - implement single point parsing?')
            energy_apo += gaulog.energies[energy_keyword_monomer][-1][-1]

    # indexes of scanned/constrained variables
    atomids = cat_modreds(gaulogs)
    if len(atomids) > 1: raise RuntimeError
    is_modred = len(atomids) > 0
    if is_modred:
        atomids = atomids.pop() # will fail if modred not in .log
        scan_atoms = ' '.join(['%d' % i for i in atomids])

    # gather coordinates and energies from multiple Gaussian .log files
    if level_of_theory.endswith('_CP'):
        energy_keyword = 'counterpoise'
        energy_keyword_monomer = 'SCF_energy'
    elif level_of_theory.startswith('B3LYP') or level_of_theory.startswith('WB97XD'):
        energy_keyword = 'SCF_energy'
        energy_keyword_monomer = 'SCF_energy'
    else:
        raise NotImplementedError

    xyz_list = []
    energies = []
    deltaE = []
    for gaulog in gaulogs:
        coords = extract_coordinates(gaulog, args.scan_step, args.opt_step, args)
        xyz_list.extend(coords)
        current_energies = [opt[-1] for opt in gaulog.energies[energy_keyword]]
        energies.extend(current_energies)
        if len(args.monomers):
            deltaE.extend([627.509 * (e - energy_apo) for e in current_energies])

    output_format = os.path.splitext(args.out)[1][1:]
    if output_format not in ['mol', 'sdf']:
        raise RuntimeError('output must be MOL/SDF')
    obconv = ob.OBConversion()
    obconv.SetOutFormat(output_format)
    mol_blocks = []
    sdf_fields = []
    n_frames = len(xyz_list)
    for i in range(n_frames):

        # mol block (analogous to .mol file)
        atom_index = 0
        for atom in ob.OBMolAtomIter(obmol):
            x, y, z = xyz_list[i][atom_index]
            atom.SetVector(x, y, z)
            atom_index += 1
        mol_block = obconv.WriteString(obmol)
        mol_blocks.append(mol_block)

        if output_format == 'sdf':
            fields  = ""
            fields += ">  <SMILES>\n%s\n\n" % smiles
            fields += ">  <MinMethod>\n%s\n\n" % level_of_theory
            fields += ">  <Energy>\n%.7f\n\n" % energies[i]
            if len(args.monomers):
                fields += ">  <deltaE>\n%.7f\n\n" % deltaE[i] # TODO subtract monomers / g16 counterpoise
            if is_modred:
                current_xyz = xyz_list[i]
                tmpcoords = [current_xyz[index_atom-1] for index_atom in atomids]
                scan_var = geom.anymetric(tmpcoords)
                fields += ">  <ScanAtoms_1>\n%s\n\n" % scan_atoms
                fields += ">  <ScanVar_1>\n%.7f\n\n" % scan_var
            sdf_fields.append(fields)
    
    # write output .sdf
    out_text = ""
    for i in range(n_frames):
        if output_format == 'sdf':
            out_text += mol_blocks[i][:-5] # remove "$$$$\n"
            out_text += sdf_fields[i]
            out_text += mol_blocks[i][-5:] #    add "$$$$\n"
        elif output_format == 'mol':
            out_text += mol_blocks[i]

    with open(args.out, 'w') as f:
        f.write(out_text)
