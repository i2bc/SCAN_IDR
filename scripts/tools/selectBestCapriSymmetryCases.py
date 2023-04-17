import os,re,sys
from pathlib import Path
import shutil
from glob import glob
import tempfile
import argparse

#
# This program calculates a consensus Capri analysis when required
# Mainly developed for 58_6GP7 but could have larger usage if required
#


#
# Specify with run types should cut the models for CAPRI evaluation 
#
list_of_possible_protocols = [
    'mixed_ali-delim-delim',
    'mixed_ali-fl-fl',
    'mixed_ali-delim-fl',
    'mixed_ali-delim-100',
    'mixed_ali-delim-200',
    'unpaired_ali-delim-delim',
    'unpaired_ali-delim-fl',
    'single_pep-delim-delim',
    'single_pep-delim-100',
    'single_pep-delim-200',
]
d_protocol = dict()
for proto in list_of_possible_protocols:
    d_protocol[proto] = dict()


def selectBest(caprieval_dir, entry):

    pdb = entry.split('_')[1]
    dir_output = caprieval_dir

    for rtype in d_protocol:
        l_capri_file = glob(f"{dir_output}/CAPRI_eval_ch*_{entry}_{rtype}.out")
        if len(l_capri_file) == 0:
            print(f"!ERROR: {entry} in capri evaluation has not the expected capri evaluation files in {caprieval_dir} for protocol {rtype}")
            continue
        fout = open(f"{dir_output}/CAPRI_eval_{entry}_{rtype}.out","w")
        fout.write("# file RMSAvsA RMSBvsB RMSABvsAB LRMSA LRMSB FRNAT FRNNAT FRIR IRMS CAPRIrank DockQ DockQrank ref_index\n")
        if len(l_capri_file) == 0:
            print(f"Nothing to do with {entry}")
            continue
        d_pdb = {}
        for fcapri in l_capri_file:
            for l in open(fcapri).readlines():
                if l[0] == "#":
                    continue
                s = l.split()
                if s[0] not in d_pdb:
                    d_pdb[s[0]] = []
                d_pdb[s[0]].append([eval(s[11]), l])
        for elt in d_pdb:
            d_pdb[elt].sort(reverse=True)
            line_best_capri = d_pdb[elt][0][1]
            fout.write(line_best_capri)
    fout.close()

def Main():
    """Main function to select the best capri scores obtained with several references
    """
    USAGE = ""

    parser = argparse.ArgumentParser(usage=USAGE)
    parser.add_argument('-e','--entry',action="store",dest='entry',type='string',
                      help="Code of the pdb example <index>_<pdbcode>"+\
                      "the best rank and irms (according to diff ref) will be reported")
    parser.add_argument('-c','--caprieval_dir',action="store",dest='caprieval_dir',type='string',
                      help="the directory where the capri evaluation file are stored ")

    if len(sys.argv) == 1:
        print(USAGE)
        print("type -h or --help for more help information")
        sys.exit(1)
    (options, args) = parser.parse_args(sys.argv)

    selectBest(options.caprieval_dir, options.entry)


if __name__ == "__main__":
    # preparation
    if len(sys.argv) == 1:
        print("type -h or --help for usage.")
        sys.exit(1)

    Main()
