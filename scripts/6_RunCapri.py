import os,re,sys
from pathlib import Path
import shutil
from glob import glob
import tempfile
import configparser

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)

config = configparser.ConfigParser()
config.read(f'{script_dir}/config.ini')

WORKING_DIR = config['DEFAULT']['WORKING_DIR']
SCRIPTS_DIR = os.path.abspath(config['DEFAULT']['SCRIPTS_DIR'])
DATA_DIR = os.path.abspath(os.path.join(WORKING_DIR, 'data', ))
AF2RUNS_DIR = os.path.abspath(os.path.join(DATA_DIR, 'af2_runs'))
FASTAMSA_DIR = os.path.abspath(os.path.join(DATA_DIR, 'fasta_msa'))
CUTMODELS_DIR = os.path.join(DATA_DIR, config['DEFAULT']['CUTMODELS_DIR'])
REFERENCE_STRUCT_DIR = os.path.join(DATA_DIR, config['DEFAULT']['REFERENCE_DIR'])

CAPRIEVAL_DIR = os.path.join(DATA_DIR, config['DEFAULT']['CAPRIEVAL_DIR'])
if not os.path.isdir(CAPRIEVAL_DIR):
    os.system(f"mkdir {CAPRIEVAL_DIR}")

CMD_CAPRI = os.path.join(SCRIPTS_DIR, "tools", "CAPRI_wrapper.py")
CMD_CAPRI_BEST_FOR_SYMMETRY = os.path.join(SCRIPTS_DIR, "tools", "selectBestCapriSymmetryCases.py")
CMD_PYTHON = config['PROGRAMS']['PYTHON']

cmdfile_runcapri = os.path.join(DATA_DIR, "cmd_run_capri_all_models.sh")
cmdfile = open(cmdfile_runcapri, "w")
cmdfile.close()

pdblist_formatted = config['DEFAULT']['INPUT_PDB_LIST']

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

f_data = [pdblist_formatted]

# We recover in d_pdbref, the START and STOP index in the reference PDB for CAPRI
# --> d_pdbref[index_pdbcode]['START'] and d_pdbref[index_pdbcode]['STOP']
d_pdbref = dict()
index = 1
for f in f_data:
    # Can be used to process multiple input files
    fpath = os.path.join(DATA_DIR, f)
    last_index, d_tmp = parse_reffile(fpath, ini_index=index)
    d_pdbref.update(d_tmp)
    index = last_index + 1

# Get a sorted list of all index_pdb directories in which the model structures are stored
l_af2_entrypath = [os.path.join(AF2RUNS_DIR, k) for k in d_pdbref]
l_entries = [[int(os.path.basename(x).split('_')[0]), x] for x in l_af2_entrypath]
l_entries = [x[1] for x in sorted(l_entries)]

for entry_index, dpath in enumerate(l_entries):
    entry = os.path.basename(dpath)
    pdb = entry.split('_')[1]
    if not os.path.isdir(os.path.join(CAPRIEVAL_DIR,entry)):
        os.system(f"mkdir {os.path.join(CAPRIEVAL_DIR,entry)}")
    print(f"# Running on pdb {pdb}")
    dir_output = os.path.join(CAPRIEVAL_DIR, entry)
    if not os.path.isdir(dir_output):
        os.system(f"mkdir {dir_output}")

    cmdfile = open(cmdfile_runcapri, "a")
    cmdfile.write(f"# Moving to directory {entry}/\ncd {os.path.join(CUTMODELS_DIR, entry)}\n")
    for protocol in list_of_possible_protocols:
        # Get the list of models to work on
        dir_models = os.path.join(CUTMODELS_DIR, entry, f"af2_{entry}_{protocol}")
        l_models = glob(os.path.join(CUTMODELS_DIR, entry, f"af2_{entry}_{protocol}","*.pdb"))
        str_models = os.path.join(CUTMODELS_DIR, entry, f"af2_{entry}_{protocol}", "*.pdb")
        is_REFopt = False
        REFmono = glob(os.path.join(REFERENCE_STRUCT_DIR, f"{pdb}_?-?.pdb"))
        if len(REFmono) > 0:
            REF = REFmono[0]
            chainRef = "".join(REF.split('/')[-1][-7:-4].split('-'))
        else:
            REFmulti = glob(os.path.join(REFERENCE_STRUCT_DIR, f"{pdb}_?-?-ref*.pdb"))
            if len(REFmulti) > 0:
                REF = ",".join(REFmulti)
                chainRef = "".join(REFmulti[0].split('/')[-1][-12:-9].split('-'))
            else:
                is_REFopt = True
                # Successive runs of CAPRI are required
                REFsucc = glob(os.path.join(REFERENCE_STRUCT_DIR, f"{pdb}_?-?-opt*.pdb"))
                l_chainRef = []
                for opt in REFsucc:
                    chainRef="".join(opt.split('/')[-1][-12:-9].split('-'))
                    l_chainRef.append(chainRef)
                 
        print(chainRef)
        #os.chdir(os.path.join(CUTMODELS_DIR, entry))

        if not is_REFopt:
            # Parallel run on 5 cpus (-p 5)
            # Important : add the -a option to prevent that failures fo files with very long filename absolute path (PROFIT related issue)
            cmd = f"{CMD_PYTHON} {CMD_CAPRI} -r {REF} -d af2_{entry}_{protocol} -c {chainRef}:AB " \
                  f"-o {dir_output}/CAPRI_eval_{entry}_{protocol}.out -a "

            cmdfile.write(cmd+"\n")

            #os.system(cmd)
        else:
            # Treatment of the symetry for 6GP7
            for ii, chainRef in enumerate(l_chainRef):
                REF = REFsucc[ii]
                for chainMob in ['B','C']:
                    # We consider the cas where the chains were not changes
                    # AF2 builds chains B,C,D and so on. Here it is tuned for 6GP7 case with 3 chains.
                    # Adapt the script for larger complexes
                    cmd = f"{CMD_PYTHON} {CMD_CAPRI} -r {REF} -d af2_{entry}_{protocol} -c {chainRef}:{chainMob}D " \
                          f"-o {dir_output}/CAPRI_eval_ch{chainRef}-{chainMob}D_{entry}_{protocol}.out -a "

                    cmdfile.write(cmd + "\n")
                    #os.system(cmd)
            
            cmd = f"{CMD_PYTHON} {CMD_CAPRI_BEST_FOR_SYMMETRY} -e {entry} -c {dir_output} " 
            cmdfile.write(cmd + "\n")   
    cmdfile.close()