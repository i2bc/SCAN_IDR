import os,re,sys
from pathlib import Path
import shutil
from glob import glob
import tempfile
import importlib
import tools_pdb as tp
importlib.reload(tp)

import configparser

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)

config = configparser.ConfigParser()
config.read(f'{script_dir}/config.ini')

WORKING_DIR = config['DEFAULT']['WORKING_DIR']
SCRIPTS_DIR = os.path.abspath(config['DEFAULT']['SCRIPTS_DIR'])
DATA_DIR = os.path.abspath(os.path.join(WORKING_DIR, 'data', ))
AF2RUNS_DIR = os.path.abspath(os.path.join(WORKING_DIR, 'data', 'af2_runs'))
FASTAMSA_DIR = os.path.abspath(os.path.join(WORKING_DIR, 'data', 'fasta_msa'))

pdblist_formatted = config['DEFAULT']['INPUT_PDB_LIST']

CUTMODELS_DIR = os.path.join(DATA_DIR, "cutmodels_for_caprieval")
if not os.path.isdir(CUTMODELS_DIR):
    os.system(f"mkdir {CUTMODELS_DIR}")

logerrpath = os.path.join(DATA_DIR, "CutModelsForCAPRI.logerr")
logerrfile = open(logerrpath, "w")
logerrfile.close()

l_exclude_change_chain = [] # In case some entries should not be processed add cound
f_data = [pdblist_formatted]

"""
###
For every type of runs, defined by the protocol extension, we need to define in a reference dictionary the rules that should be applied on the delimitations 
The models created by AF2 are always numbered from residue 1 to end for every chains
To properly assess the prediction, we need first to adjust the delimitations of the modelled chains so that they match with the reference PDB structure
###

1 - We recover the PDB delimitations from the file cmd_xxx_mixed_ali-delim-delim.sh in af2_runs directory because they match the original PDB
      ==> They will be applied to recover a pdb model that can be compared directly to CAPRI

2 - If MODEL is full-length we can use directly the delimitations above to cut the pdb 

3 - If the MODEL is extended for the IDP (+100, +200 or more ), we need to remove only the residues added to create the extension
       Since they are numbered from 1 to N, we need first to renumber the PDBchains from DELIMstart to DELIMstop and only after apply the PDB cutter with the correct index

4 - In case of multiple chains for the receptor :  We fuse and renumber so they can be compared with the CAPRI reference prepared in the same manner.
       If homo dimer we need to define 2 different CAPRI references.
"""

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

# 0 means do no cut. 1 means cut from 1. 2 means renumber pdb first and cut.
# The first index is for receptor, the second for the ligand
d_protocol['mixed_ali-delim-delim']['cutRecIDP'] = [0, 0]
d_protocol['single_pep-delim-delim']['cutRecIDP'] = [0, 0]
d_protocol['mixed_ali-fl-fl']['cutRecIDP'] = [1, 1]
d_protocol['mixed_ali-delim-fl']['cutRecIDP'] = [0, 1]
d_protocol['unpaired_ali-delim-delim']['cutRecIDP'] = [0, 0]
d_protocol['unpaired_ali-delim-fl']['cutRecIDP'] = [0, 1]
d_protocol['mixed_ali-delim-100']['cutRecIDP'] = [0, 2]
d_protocol['mixed_ali-delim-200']['cutRecIDP'] = [0, 2]
d_protocol['single_pep-delim-100']['cutRecIDP'] = [0, 2]
d_protocol['single_pep-delim-200']['cutRecIDP'] = [0, 2]

for proto in list_of_possible_protocols:
    d_protocol[proto]['cmdfile'] = proto

def extract_delim(cmdfile):
    """
    Get the start and stop values of DELIM1, DELIM2, ...
    """
    pat = re.compile("DELIM(\d+)='?(\d*?)_'?(\d*?)[\s\n]")
    fin = open(cmdfile)
    l = pat.findall(fin.read())
    return l

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

#
# Loop over all the directories and models PDB
#

#possible chains in AF2 models (last one is IDP).
# TODO: Checking. Usually AF2 models start at B but to be changed if it's not the case anymore
chain_models="BCDEFGHIJKLMNOPQRSTUVWXYZ"

for entry_index, dpath in enumerate(l_entries):
    # dpath is the full path of the af2_run directory
    entry = os.path.basename(dpath)
    pdb = entry.split('_')[1]
    if not os.path.isdir(os.path.join(CUTMODELS_DIR, entry)):
        os.mkdir(os.path.join(CUTMODELS_DIR, entry))

    print(f"# Running on pdb {entry}")
    # Get the list of dict containing chainid, uniprotid, starts and stops of every chains in CAPRI ref
    l_data_ref = d_pdbref[entry]
    l_starts_capri = [x['START'] for x in l_data_ref]
    l_stops_capri = [x['STOP'] for x in l_data_ref]
    n_chain = len(l_data_ref)
    l_chain_receptor = list(chain_models[:n_chain-1])
    l_chain_IDP = list(chain_models[n_chain-1])

    for protocol in d_protocol:

        print(f"## Running in protocol {protocol}")
        
        # Get the list of models to work on
        pat_revision = re.compile("af2_predictions.*?_v(\d+)")
        l_models = glob(os.path.join(dpath, f"af2_{entry}_{protocol}", "af2_predictions*_v*", "*.pdb"))

        cmdfile = os.path.join(dpath, f"cmd_{entry}_{d_protocol[protocol]['cmdfile']}.sh")
        print(cmdfile)
        
        # Define the path where the cut models will be stored
        protocol_model_cut_dir = os.path.join(CUTMODELS_DIR, entry, f"af2_{entry}_{protocol}")

        # Get the protocol index to be applied for cutting the chains in the model
        cutRec = d_protocol[protocol]['cutRecIDP'][0]
        cutIDP = d_protocol[protocol]['cutRecIDP'][1]
        
        #if cutIDP == 2:
        ## Extracts the delimitations used for the run reading the command file
        l_delim_in_cmdfile = extract_delim(cmdfile)
        #print(l_delim_in_cmdfile)
        
        # Cleanup the previous models cut for capri
        if os.path.isdir(protocol_model_cut_dir):
            shutil.rmtree(protocol_model_cut_dir)
        os.mkdir(protocol_model_cut_dir)
        if not os.path.isdir(protocol_model_cut_dir):
            break
        
        for model in l_models:
            print(f"### Running on model {os.path.basename(model)}")

            # Inside every directory, it is not clear whether the pdb are indexed by the version dir in which they were created.
            # We change the name of the pdb so that it incorporates the version number.
            m = pat_revision.search(model)
            if m:
                version = m.group(1)
            else:
                print(f"Could not get the version of the model in {model}")
                logerrfile = open(logerrpath, "a")
                logerrfile.write(f"ERROR: model version NOT FOUND {entry}\t{protocol}\t{model}\n")
                logerrfile.close()
            
            model_basename = os.path.basename(model)
            #root, ext = os.path.splitext(model_basename)
            #out_model_basename = f'{root}_v{version}{ext}'
            out_model_basename = model_basename
            fout_model_name = os.path.join(protocol_model_cut_dir, out_model_basename)
            
            # Generate a tmp file for every chain
            l_tmpchain, _ = tp.splitPDB(model, l_chain_receptor+l_chain_IDP)

            # Process on both receptor and IDP
            for ii, cut_protocol in enumerate([cutRec,cutIDP]):
                
                if ii == 0:
                    # For first chain 
                    list_chains = l_tmpchain[:-1]
                    new_chain = "A"
                    l_start = l_starts_capri[:-1]
                    l_stop = l_stops_capri[:-1]
                    current_chain = l_chain_receptor
                    fpdb_receptor = l_tmpchain[0]
                    l_limit_in_cmdfile = l_delim_in_cmdfile[:-1]
                    
                if ii == 1:
                    # For second chain
                    #print("#### Running the IDP")
                    list_chains = [l_tmpchain[-1]]
                    new_chain = "B"
                    l_start = [l_starts_capri[-1]]
                    l_stop = [l_stops_capri[-1]]
                    current_chain = l_chain_IDP
                    fpdb_ligand = l_tmpchain[-1]
                    l_limit_in_cmdfile = l_delim_in_cmdfile[-1]
                    
                # Renumbering step    
                if cut_protocol == 2:
                    model_start = l_limit_in_cmdfile[1]
                    model_stop = l_limit_in_cmdfile[2]
                    for jj, fch in enumerate(list_chains): # list_chains is a list of files

                        cmd_status = tp.renumberPDB(fch, current_chain[jj], model_start)
                        if cmd_status != 0:
                            logerrfile = open(logerrpath, "a")
                            logerrfile.write(f"ERROR: Problem with ChainRenumber {entry}\t{protocol}\t{model}\n")
                            logerrfile.close()
                        last_res_index = tp.getLastResidueIndex(fch)
                        if last_res_index != model_stop:
                            print("ERROR: Inconsistent values between the DELIM in cmd file and in the AF2 model")
                            print(f"Model is expected to start at resi {model_start} and end at resi {model_stop}")
                            print(f"Renumbering the chain for CAPRI evaluation we found it ends at {last_res_index}")
                            print("Potentially highly detrimental to the interpretation. Please correct first !!")
                            if entry_index != 4:
                                sys.exit()
                            else:
                                print("Entry 5_6DO3 has a selenocystein 'U' residue in position. "
                                      "That's why there is a shift by one residue. "
                                      "Check it is not too detrimental for CAPRI eval.")
                
                # Once renumbered we cut in a consistent manner        
                if cut_protocol == 1 or cut_protocol == 2:
                    # residues index in the PDB are consistent with FL protein => cut using the capri delimitations
                    for jj,fch in enumerate(list_chains):
                        # cut every chain file
                        cmd_status = tp.cutPDB(fch, l_start[jj], l_stop[jj], current_chain[jj])
                
                # Change change name
                for jj,fch in enumerate(list_chains): # list_chains is a list of files
                    if entry not in l_exclude_change_chain: 
                        cmd_status = tp.alterchainPDB(fch, current_chain[jj], new_chain)
                        #print(cmd_status)
                        if cmd_status != 0:
                            logerrfile = open(logerrpath, "a")
                            logerrfile.write(f"ERROR: Inconsistent Chains {entry}\t{protocol}\t{model}\n")
                            logerrfile.close()

                if len(list_chains) > 1:
                    first_chain_in_PDB = list_chains[0]
                    last_res_index = tp.getLastResidueIndex(first_chain_in_PDB)
                    for jj,fch in enumerate(list_chains[1:]):
                        # We concatenate the chains into a single pdb
                        if entry not in l_exclude_change_chain: # If we don't change chain names no need to renumber
                            cmd_status = tp.renumberPDB(fch, new_chain, str(int(last_res_index)+1))
                            if cmd_status != 0:
                                logerrfile = open(logerrpath, "a")
                                logerrfile.write(f"ERROR: Problem with ChainRenumber {entry}\t{protocol}\t{model}\n")
                                logerrfile.close()
                        tp.fusePDB([first_chain_in_PDB, fch], fused_filename=first_chain_in_PDB)
                        if entry not in l_exclude_change_chain: # If we don't change chain names no need to renumber
                            last_res_index = tp.getLastResidueIndex(first_chain_in_PDB)
                            
            # Concatenate receptor and IDP ligand file
            _ = tp.fusePDB([fpdb_receptor, fpdb_ligand], fout_model_name, add_ter=True, add_end=True) 
            
            # Cleanup everything
            for elt in l_tmpchain: os.remove(elt)
logerrfile.close()
