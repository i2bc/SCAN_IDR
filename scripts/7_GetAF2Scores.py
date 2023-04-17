import os,re,sys
from pathlib import Path
import shutil
from glob import glob
import tempfile
import copy
import configparser
import argparse

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)

config = configparser.ConfigParser()
config.read(f'{script_dir}/config.ini')

WORKING_DIR = config['DEFAULT']['WORKING_DIR']
SCRIPTS_DIR = os.path.abspath(config['DEFAULT']['SCRIPTS_DIR'])
DATA_DIR = os.path.abspath(os.path.join(WORKING_DIR, 'data', ))
AF2RUNS_DIR = os.path.abspath(os.path.join(DATA_DIR, 'af2_runs'))
FASTAMSA_DIR = os.path.abspath(os.path.join(DATA_DIR, 'fasta_msa'))

CUTMODELS_DIR = os.path.abspath(os.path.join(DATA_DIR, config['DEFAULT']['CUTMODELS_DIR']))
CAPRIEVAL_DIR = os.path.abspath(os.path.join(DATA_DIR, config['DEFAULT']['CAPRIEVAL_DIR']))

pdblist_formatted = os.path.abspath(os.path.join(DATA_DIR, config['DEFAULT']['INPUT_PDB_LIST']))

RESULTS_DIR = os.path.abspath(os.path.join(DATA_DIR, config['DEFAULT']['RESULTS_DIR']))
if not os.path.isdir(RESULTS_DIR):
    os.system(f"mkdir {RESULTS_DIR}")
    
fscorepath = os.path.abspath(os.path.join(RESULTS_DIR, config['DEFAULT']['OUTPUT_AF2SCORES']))

n_iterations = int(config['DEFAULT']['N_ITERATIONS'])

###
###  Recovers AF2 scores in the af2_runs directory.
###

def parse_reffile(fpath, ini_index=1):
    """
    Parse the reference file Curated indicated the uniprot and the delimitations of every pdb input
    """
    dref = dict()
    count = ini_index - 1
    pat = re.compile("CHAIN:(\w+).+?UNIPROT:(\S+).+?START:(\d+).+?STOP:(\d+)")
    for l in open(fpath).readlines():
        if l[0] == "#":
            pdb = l.split()[2]
            count += 1
            dref[f"{count}_{pdb}"] = []
        else:
            m = pat.search(l)
            if not m:
                print(f"ERROR: Problem in {fpath} of delimitation in pdb {pdb}")
                sys.exit()
            dref[f"{count}_{pdb}"].append({'CHAIN':m.group(1), 'UNIPROT':m.group(2), 'START':m.group(3), 'STOP':m.group(4)})
    return count, dref

f_data = [pdblist_formatted]

#
# Creation of a dictionary listing the protocols
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

# We recover in d_pdbref, the START and STOP index in the reference PDB for CAPRI
# --> d_pdbref[pdbcode]['START'] and d_pdbref[pdbcode]['STOP']
d_pdbref = dict()
index = 1
for f in f_data:
    # Can be used to process multiple input files
    fpath = os.path.join(f)
    last_index, d_tmp = parse_reffile(fpath, ini_index=index)
    d_pdbref.update(d_tmp)
    index = last_index + 1

# Get a sorted list of all index_pdb directories in which the model structures are stored
l_af2_entrypath = [os.path.join(AF2RUNS_DIR, k) for k in d_pdbref]
l_entries = [[int(os.path.basename(x).split('_')[0]), x] for x in l_af2_entrypath]
l_entries = [x[1] for x in sorted(l_entries)]



#pat = re.compile("model_(\d+) took .+? recycles. with pLDDT (\d+.?\d*) and ptmscore (\d+.\d+)")
pat = re.compile("model_(\d+) took .+? recycles. with pLDDT (\d+.?\d*).+?ptmscore (\d+.\d+).+?iptm\w*? (\d+.\d+)")
pat2 = re.compile("length (\d+)")

fout = open(fscorepath, "w")
fout.write('#entry\tprotocol\tversion\tNmodel\tpLDDT\tpTMScore\tipTMScore\tCombinedScore\tLength\n')
for entry_index, dpath in enumerate(l_entries):
    entry = os.path.basename(dpath)
    pdb = entry.split('_')[1]
    index_pdb = int(entry.split('_')[0])
    

    ali_nb_positions = None
    #set_types = set(copy.deepcopy(list(d_runtype.keys())))
    for version in range(n_iterations):
        set_types = set(copy.deepcopy(list(d_protocol.keys())))
        for protocol in d_protocol:
            #print(rtype)
            flog_score_path = os.path.join(dpath,f'af2_{entry}_{protocol}',f'af2_predictions_v{version+1}/log.txt')
            #print(flog_score_path)
            #flog_score_path_bckup = os.path.join(dpath,f'af2_{entry}_{protocol}',f'af2_predictions_*v{version+1}/af2_predictions/log.txt')
            flog_score = glob(flog_score_path)
            with open(flog_score[0]) as f:
                fin = f.readlines()
                Length = 0
                for l in fin:
                    m2 = pat2.search(l)
                    if m2:
                        Length = m2.group(1)
                    m = pat.search(l)
                    if m:
                        set_types.discard(protocol)
                        combined_score = float(m.group(3)) * 0.2 + float(m.group(4)) * 0.8 # combined score is ptmscore * 0.2 + iptmscore * 0.8
                        fout.write(f"{entry}\t{protocol}\t{version+1}\t{m.group(1)}\t{m.group(2)}\t{m.group(3)}\t{m.group(4)}\t{combined_score:.3f}\t{Length}\n")
                        

        if len(set_types) > 0:
            print("Some scores could not be recorded")
            print(pdb, version, set_types)
        else:
            print(f"{entry}, iteration {version}, all scores recovered")
fout.close()
