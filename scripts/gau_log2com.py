#!/usr/bin/env python
"""
Reads a structure from a Gaussian Log file and other parameters from a
Gaussian Com file and creates a new Gaussian Com.

"""

import argparse
import sys
import os

# qt_scripts modules
from omg.gaussian.gaussian import GaussianLog
from omg.gaussian.gaussian import GaussianCom
from omg.misc import increment_filename

def get_args():
    "Parse arguments of gau_log2com.py"
    parser = argparse.ArgumentParser(
        description="""
            Reads a structure from a Gaussian Log file and other parameters from a  
            Gaussian Com file and creates a new Gaussian Com.""",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('log', help='gaussian log filename')

    parser.add_argument('-s', '--scan_step',
                        help='scan step (default = last = 0)',
                        default=0,
                        type=int)
    parser.add_argument('-o', '--opt_step',
                        help='opt step  (default = last = 0, -1 for AUTO)',
                        default=0)
    parser.add_argument('--new_com',
                        help='new com name (default = increment template com)',
                        default='')
    parser.add_argument('--template_com',
                        help='(default = same name as .log)',
                        default='')

    args = parser.parse_args()
    if args.template_com == '':
        args.template_com = os.path.splitext(args.log)[0] + '.com'
    if args.new_com == '':
        args.new_com = increment_filename(args.template_com)
    if os.path.exists(args.new_com):
        sys.stderr.write(
            '{0} already exists. Aborting.\n'.format(args.new_com))
        sys.exit(2)
    args.scan_step -= 1
    try:
        args.opt_step = int(args.opt_step)
        args.opt_step -= 1
    except:
        if args.opt_step != 'AUTO':
            sys.stderr.write('--opt_step must be integer or "AUTO"\n')
            sys.exit()
            
    return args

def find_most_converged(gaulog, scan_step):
    """ find opt step with lowest maxforce"""
    lowest_maxforce = float('inf')
    index_best = None
    index = 0
    for byte in gaulog.bytedict['Converged?'][scan_step]:
        (labels, values, thresh) = gaulog.read_converged(byte)
        maxforce, avgforce, maxdisp, avgdisp = values
        if maxforce < lowest_maxforce:
            lowest_maxforce = maxforce
            index_best = index
        index += 1
    return index_best

def main():
    """
    Reads a structure from a Gaussian Log file and other parameters from a
    Gaussian Com file and creates a new Gaussian Com.
    """
    args = get_args()
    gaulog = GaussianLog(args.log)
    gaucom = GaussianCom(args.template_com)
    if args.opt_step == 'AUTO':
        args.opt_step = find_most_converged(gaulog, args.scan_step)
    atoms_log_coords = gaulog.read_geometry(args.opt_step, args.scan_step)
    for no, atom in enumerate(gaucom.atoms_list):
        atom.SetVector(atoms_log_coords[no].GetVector())
    gaucom.write_to_file(args.new_com)

if __name__ == "__main__":
    main()


