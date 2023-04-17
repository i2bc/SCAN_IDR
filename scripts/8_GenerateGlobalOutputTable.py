import os,re,sys
from pathlib import Path
import shutil
from glob import glob
import tempfile
import copy
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

CUTMODELS_DIR = os.path.abspath(os.path.join(DATA_DIR, config['DEFAULT']['CUTMODELS_DIR']))
CAPRIEVAL_DIR = os.path.abspath(os.path.join(DATA_DIR, config['DEFAULT']['CAPRIEVAL_DIR']))

pdblist_formatted = config['DEFAULT']['INPUT_PDB_LIST']


pdblist_formatted = os.path.abspath(os.path.join(DATA_DIR, config['DEFAULT']['INPUT_PDB_LIST']))

RESULTS_DIR = os.path.abspath(os.path.join(DATA_DIR, config['DEFAULT']['RESULTS_DIR']))
fscorepath = os.path.abspath(os.path.join(RESULTS_DIR, config['DEFAULT']['OUTPUT_AF2SCORES']))

n_iterations = int(config['DEFAULT']['N_ITERATIONS'])

fglobaloutput_path = os.path.abspath(os.path.join(RESULTS_DIR, config['DEFAULT']['OUTPUT_GLOBAL']))

"""

Scripts which recovers the data describing the length of the inputs and the AF2 scores
        parses the CAPRI evaluation files
        generates a global table with all the information

"""


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

def capri_rank_for_peptide(fnat, lrmsd, irmsd):
    if fnat >= 80.0 and (lrmsd <= 1.0 or irmsd <= 0.5):
        rank = 'High'
    elif ((fnat >= 50.0 and fnat < 80.0) and (lrmsd <= 2.0 or irmsd <= 1.0)) or ((fnat >= 80.0) and (lrmsd > 1.0 and irmsd > 0.5)) :
        rank = 'Medium'
    elif ((fnat >= 20.0 and fnat < 50.0) and (lrmsd <= 4.0 or irmsd <= 2.0)) or ((fnat >= 50.0) and (lrmsd > 2.0 and irmsd > 1)) :
        rank = 'Acceptable'
    else:
        rank = 'Incorrect'
    return rank



f_data = [pdblist_formatted]
d_global_res = dict()

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
    fpath = os.path.join(DATA_DIR, f)
    last_index, d_tmp = parse_reffile(fpath, ini_index=index)
    d_pdbref.update(d_tmp)
    index = last_index + 1
   
# Get a sorted list of all index_pdb directories in which the model structures are stored
l_af2_entrypath = [os.path.join(AF2RUNS_DIR, k) for k in d_pdbref]
l_entries = [[int(os.path.basename(x).split('_')[0]), x] for x in l_af2_entrypath]
l_entries = [x[1] for x in sorted(l_entries)]

# Parse the file listing the AF2 scores
d_score = dict()
fscore = open(fscorepath, 'r').readlines()
counter = 0

for l in fscore:
    if l[0] == "#":
        continue
    s = l.split()
    entry = s[0]
    protocol = s[1]
    version = int(s[2])
    model_index = s[3]
    pLDDT = s[4]
    pTMSCORE = s[5]
    ipTMSCORE = s[6]
    CombinedSCORE = s[7]
    Length_in_scorefile = s[8]
    if entry not in d_score:
        d_score[entry] = {}
    if protocol not in d_score[entry]:
        d_score[entry][protocol] = {}
    if version not in d_score[entry][protocol]:
        d_score[entry][protocol][version] = {}
    d_score[entry][protocol][version][model_index] = {}
    d_score[entry][protocol][version][model_index]['pLDDT'] = pLDDT
    d_score[entry][protocol][version][model_index]['pTMSCORE'] = pTMSCORE
    d_score[entry][protocol][version][model_index]['ipTMSCORE'] = ipTMSCORE
    d_score[entry][protocol][version][model_index]['CombinedSCORE'] = CombinedSCORE
    d_score[entry][protocol][version][model_index]['Length_in_scorefile'] = Length_in_scorefile
    
pat = re.compile("rank_(\d+)_model_(\d+)_v(\d+)")
fout = open(fglobaloutput_path, "w")
fout.write('# $1:index_pdb\n# $2:protocol\n# $3:SizeModel\n# $4:version\n# $5:modelIndex\n# $6:AF2Rank\n# $7:pLDDT\n# $8:pTMSCORE\n# $9:ipTMSCORE\n# $10:CombinedSCORE\n# $11:DockQ\n# $12:DockQrank\n# $13:CAPRIrank\n# $14:CAPRIrankPeptide\n# $15:IRMS\n# $16:RMSAvsA\n# $17:RMSBvsB\n# $18:RMSABvsAB\n# $19:LRMSA\n# $20:LRMSB\n# $21:FRNAT\n# $22:FRNNAT\n# $23:FRIR\n')
for entry_index,dpath in enumerate(l_entries):
    entry = os.path.basename(dpath)
    d_global_res[entry] = dict()
    pdb = entry.split('_')[1]
    index_pdb = int(entry.split('_')[0])
    ali_nb_positions = None
    #set_types = set(copy.deepcopy(list(d_protocol.keys())))
    #for version in range(nmodels):
    
    set_types = set(copy.deepcopy(list(d_protocol.keys())))
    for protocol in d_protocol:
        #d_global_res[entry][version][protocol] = dict()
        
        #print(protocol)
        fcapri_eval_path = os.path.join(CAPRIEVAL_DIR, f'{entry}', f'CAPRI_eval_{entry}_{protocol}.out')
        with open(fcapri_eval_path) as f:
            fin = f.readlines()                
            for l in fin:
                if l[0] == "#":
                    continue
                s = l.split()
                m = pat.search(l)
                if m:
                    version = int(m.group(3))
                    if version not in d_global_res[entry]:
                        d_global_res[entry][version] = dict()
                    if protocol not in d_global_res[entry][version]:
                        d_global_res[entry][version][protocol] = dict()
                    Rank        = m.group(1)
                    d_global_res[entry][version][protocol][Rank] = dict()
                    dgr = d_global_res[entry][version][protocol][Rank]
                    dgr['model_index'] = m.group(2)
                    model_index = dgr['model_index']
                    dgr['RMSAvsA']     = s[1]
                    dgr['RMSBvsB']     = s[2]
                    dgr['RMSABvsAB']   = s[3]
                    dgr['LRMSA']       = s[4] 
                    dgr['LRMSB']       = s[5]
                    dgr['FRNAT']       = s[6]
                    dgr['FRNNAT']      = s[7]
                    dgr['FRIR']        = s[8]
                    dgr['IRMS']        = s[9]
                    dgr['CAPRIrank']   = s[10]
                    dgr['DockQ']       = s[11]
                    dgr['DockQrank']   = s[12]
                    
                    CAPRIrank_peptide = capri_rank_for_peptide(eval(dgr['FRNAT']), eval(dgr['LRMSA']), eval(dgr['IRMS']))
                    dgr['CAPRIrankPeptide'] = CAPRIrank_peptide

                    dgr['pLDDT']       = d_score[entry][protocol][version][model_index]['pLDDT']
                    dgr['pTMSCORE']    = d_score[entry][protocol][version][model_index]['pTMSCORE']
                    dgr['ipTMSCORE']        = d_score[entry][protocol][version][model_index]['ipTMSCORE']
                    dgr['CombinedSCORE']    = d_score[entry][protocol][version][model_index]['CombinedSCORE']
                    dgr['SizeModel']   = d_score[entry][protocol][version][model_index]['Length_in_scorefile']
                    fout.write(f"{entry}\t{protocol}\t{dgr['SizeModel']}\t{version}\t{dgr['model_index']}\t{Rank}\t"
                               f"{dgr['pLDDT']}\t{dgr['pTMSCORE']}\t{dgr['ipTMSCORE']}\t{dgr['CombinedSCORE']}\t"
                               f"{dgr['DockQ']}\t{dgr['DockQrank']}\t{dgr['CAPRIrank']}\t{dgr['CAPRIrankPeptide']}\t"
                               f"{dgr['IRMS']}\t{dgr['RMSAvsA']}\t{dgr['RMSBvsB']}\t{dgr['RMSABvsAB']}\t"
                               f"{dgr['LRMSA']}\t{dgr['LRMSB']}\t{dgr['FRNAT']}\t{dgr['FRNNAT']}\t{dgr['FRIR']}\n")
                else:
                    print(f"Problem capri eval line could not be parsed :\n{l}")
                    sys.exit()
                        
fout.close()
