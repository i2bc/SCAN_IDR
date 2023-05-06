import os,re,sys
from pathlib import Path
import shutil
import configparser
import argparse
import tempfile

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)

config = configparser.ConfigParser()
config.read(f'{script_dir}/config.ini')

WORKING_DIR = os.path.abspath(config['DEFAULT']['WORKING_DIR'])
SCRIPTS_DIR = os.path.abspath(config['DEFAULT']['SCRIPTS_DIR'])
PYTHON_CMD = config['PROGRAMS']['PYTHON']
COLABFOLD_CMD = config['PROGRAMS']['COLABFOLD_SIF']
SINGULARITY_CMD = config['PROGRAMS']['SINGULARITY']

DATA_DIR = os.path.abspath(os.path.join(WORKING_DIR, 'data', ))
AF2RUNS_DIR = os.path.abspath(os.path.join(WORKING_DIR, 'data', 'af2_runs'))
FASTAMSA_DIR = os.path.abspath(os.path.join(WORKING_DIR, 'data', 'fasta_msa'))

pdblist_formatted = config['DEFAULT']['INPUT_PDB_LIST']




if not os.path.isdir(AF2RUNS_DIR):
    os.system(f"mkdir {AF2RUNS_DIR}")

def assignProtocolInformation(protocol_name):
    d_proto = dict()
    split_dash = protocol_name.split('-')
    d_proto['tag'] = protocol_name
    d_proto['ali'] = split_dash[0]
    d_proto['rec_ext'] = split_dash[1]
    d_proto['lig_ext'] = split_dash[2]
    return d_proto

def get_protein_size(path2fasta):
    """
    Get the full length of the protein from its fasta
    """
    fin = open(path2fasta).readlines()
    size = 0
    for ll in fin[1:]:
        s = ll.split()[0]
        size += len(s)
    return size

def fillGapLigand(start, stop, fin, fout):
    """
    Add Gaps in the end of the sequences of an alignment. Only works for ligand gapping.
    """
    nbgap=int(stop) -int(start)
    #print(nbgap)
    # read file
    with open(fin, 'r') as fp:
        # read an store all lines into list
        lines = fp.readlines()
    # Write file
    with open(fout, 'w') as fp:
        # iterate each line
        lineprec = ''
        for number, line in enumerate(lines):
            if number > 2 and line[0] != '>' and '>10' not in lineprec and "\t10" not in lineprec:
                line_B = len(line)
                line = line.replace('\n', '')
                line = line[0:-1*nbgap-1]+(nbgap+1)*'-'+'\n' # add the gaps at the end of the sequence
                if len(line) != line_B:
                    print('ERROR len line', len(line), line_B )
            fp.write(line)
            lineprec = line

def get_new_delimitations(start_ini, stop_ini, extension, protein_size):
    """
    Extends the start and stop delimitation for a given extension

    """
    index_start_IDP = int(start_ini)
    index_stop_IDP = int(stop_ini)
    half_extension = int(extension/2)

    if index_start_IDP - half_extension < 1:
        length_to_be_added_above = half_extension - index_start_IDP
        new_start = 1
    else:
        # new_stop = min(protein_size, index_stop_IDP + half_IDP_length)
        new_start = index_start_IDP - half_extension
        length_to_be_added_above = 0

    if index_stop_IDP + half_extension > protein_size:
        length_to_be_added_below = half_extension - (protein_size - index_stop_IDP)
        new_stop = protein_size
    else:
        length_to_be_added_below = 0
        new_stop = index_stop_IDP + half_extension

    new_start = max(1, new_start - length_to_be_added_below)
    new_stop = min(protein_size, new_stop + length_to_be_added_above)
    return new_start, new_stop

def change_singu_colab_cmd(runcolab_template, tmp_runcolab, colabfoldpath, singularitypath):
    """
    Changes the path to the colabfold and singularity command following the content of the config.ini
    """
    with open(runcolab_template) as fin:
        with open(tmp_runcolab, 'w') as fout:
            for ll in fin.readlines():
                newll = None
                sp1 = ll.split('=')
                sp2 = ll.split()
                if len(sp1) > 0 and sp1[0] == "COLAB_VERSION":
                    newll = f"COLAB_VERSION={colabfoldpath}\n"
                elif len(sp2) > 1:
                    if sp2[1] == "exec":
                        singu_execline = ' '.join(sp2[1:])
                        newll = singularitypath + " " + singu_execline + "\n"
                    else:
                        newll = ll
                else:
                    newll = ll
                fout.write(newll)
    return



def changeCmdLine_Npartners(filename, l_start, l_stop, l_uniprot, number_of_partners,
                            co_msa_mode='mixed_ali', rec_extension='delim', lig_extension='delim'):
    """
    Method used in case two chains vs 1 chain
    co_msa_mode can be : 'mixed_ali', 'unpaired_ali', 'single_pep'
    rec_extension can be : 'delim', 'fl', '100' or '200'
    lig_extension can be : 'delim', 'fl', '100' or '200'
    """
    l_applied_delimitations = [] # list that will be returned and which will contain the applied delimitations
    list_indexes = [str(x) for x in list(range(number_of_partners))]
    
    list_partners_indexes = ",".join(list_indexes)
    
    # Variable used to write the characters such as  "${delim_list[0]},${delim_list[1]},${delim_list[2]}" in line 
    list_delims_tmp = ["delim_list[{index}]".format(index=ii) for ii in list_indexes]
    list_delims = ",".join(["${"+strx+"}" for strx in list_delims_tmp]) # written as this, will refer to the values in the DELIM lines

    print(l_start, l_stop)
    # read file
    with open(filename, 'r') as fp:
        # read an store all lines into list
        lines = fp.readlines()
    # Write file
    with open(filename, 'w') as fp:
        # iterate each line
        for number, line in enumerate(lines):
            if line.startswith("DELIM"):
                partnernb = eval(line.split('=')[0][5:])
                is_ligand = False
                if partnernb == len(l_start):
                    is_ligand = True
                print("partnernb", partnernb)
                if rec_extension in ['fl', '100', '200'] or lig_extension in ['fl', '100', '200']:
                    geneid = l_uniprot[partnernb - 1]
                    path2fasta = os.path.join(FASTAMSA_DIR, geneid+".fasta")
                    protein_size = get_protein_size(path2fasta)

                if is_ligand:
                    extension_type = lig_extension
                else:
                    extension_type = rec_extension

                # We rewrite the DELIM<index> line according to the delimitations rules specified
                if extension_type == 'delim':
                    start_in_delim = l_start[partnernb - 1]
                    stop_in_delim = l_stop[partnernb - 1]
                elif extension_type == 'fl':
                    start_in_delim = 1
                    stop_in_delim = protein_size
                else:
                    try:
                        int_extension = int(extension_type)
                        new_start, new_stop = get_new_delimitations(l_start[partnernb-1], l_stop[partnernb-1], int_extension, protein_size)
                        start_in_delim = new_start
                        stop_in_delim = new_stop
                    except:
                        print("Problem with the options used in changeCmdLine_Npartners. "
                              "Integer is expected to define the extension length")
                        sys.exit()
                line = f"DELIM{partnernb}={start_in_delim}_{stop_in_delim}\n"
                l_applied_delimitations.append([start_in_delim, stop_in_delim])

            if line.startswith("DO_MATCH_SPECIES"):
                # Set to true only once, the first time, supposed to be for the mixed_ali-deli-delim protocol.
                # Definition of common species for pairing is the same for all the protocols so we don't need to repeat
                if co_msa_mode == 'mixed_ali' and rec_extension == 'delim' and lig_extension == 'delim':
                    line = "DO_MATCH_SPECIES=true\n"
                else:
                    if os.path.isfile("json_step1.json"):
                        # Normally only for the first protocol, the json containing the matched species should be calculated
                        line = "DO_MATCH_SPECIES=false\n"
                    else:
                        # In case the first protocol was not run, the json_step1 should be generated
                        line = "DO_MATCH_SPECIES=true\n"

            if "${PYTHON} ${ALPHASCRIPTS}/GetAlphaFoldInputs.py" in line and "mixed" in line:
                line = "  ${PYTHON} ${ALPHASCRIPTS}/GetAlphaFoldInputs.py -msa "
                line += list_partners_indexes
                line += " --delimitation "
                line += list_delims
                line += " --jsonname json_step1 --jsonoutput json_step2 --colab -st ${STOICHIOMETRY} -o mixed\n"

            if "${PYTHON} ${ALPHASCRIPTS}/GetAlphaFoldInputs.py" in line and "all_unpaired" in line:
                line = "  ${PYTHON} ${ALPHASCRIPTS}/GetAlphaFoldInputs.py -msa "
                line += list_partners_indexes
                line += " --delimitation "
                line += list_delims
                line += " --jsonname json_step1 --jsonoutput json_step2 --colab -st ${STOICHIOMETRY} -o all_unpaired\n"

            if "RUN_PAIRED_ALI_OPTION=" in line:
                line = "RUN_PAIRED_ALI_OPTION=false\n"
            if "RUN_MIXED_ALI_OPTION=" in line:
                if co_msa_mode in ['mixed_ali', 'single_pep']:
                    line = "RUN_MIXED_ALI_OPTION=true\n"
                elif co_msa_mode in ['unpaired_ali']:
                    line = "RUN_MIXED_ALI_OPTION=false\n"
            if "RUN_ALLUNPAIRED_ALI_OPTION=" in line:
                if co_msa_mode in ['mixed_ali', 'single_pep']:
                    line = "RUN_ALLUNPAIRED_ALI_OPTION=false\n"
                elif co_msa_mode in ['unpaired_ali']:
                    line = "RUN_ALLUNPAIRED_ALI_OPTION=true\n"
            fp.write(line)
    return l_applied_delimitations


def createCmdColab(SCRIPTS_DIR, FASTAMSA_DIR,
                   count, pdb, tag, l_uniprot, l_start, l_stop,
                   co_msa_mode='mixed_ali', rec_extension='delim', lig_extension='delim'):
    """
    Prepare the directory where af2 will be run. The concatenated alignment and the command to run the colab
    """
    for uniprot in l_uniprot:
        fasta = f"{uniprot}.fasta"
        msadir = f"msa_{uniprot}"
        if os.path.exists(fasta):
            os.remove(fasta)
        os.symlink(os.path.join(FASTAMSA_DIR, fasta), fasta)
        if os.path.exists(msadir):
            os.remove(msadir)
        os.symlink(os.path.join(FASTAMSA_DIR, msadir), msadir)

    str_uniprots = ','.join(l_uniprot)
    cmd = f"{PYTHON_CMD} {SCRIPTS_DIR}/createCmdFile.py -l {str_uniprots} -o cmd_{count}_{pdb}_{tag}.sh -s final -a {SCRIPTS_DIR} -k {co_msa_mode}"
    print(cmd)
    os.system(cmd)
    l_applied_delimitations = changeCmdLine_Npartners(f"cmd_{count}_{pdb}_{tag}.sh", l_start, l_stop, l_uniprot, len(l_start),
                            co_msa_mode=co_msa_mode, rec_extension=rec_extension, lig_extension=lig_extension)
    os.system(f"sh cmd_{count}_{pdb}_{tag}.sh")
    af2_dir = f"af2_{count}_{pdb}_{tag}"
    if not os.path.isdir(af2_dir):
        os.mkdir(af2_dir)
    if co_msa_mode in ['mixed_ali']:
        shutil.copy(f"output/MSAforAlphafold_mixed.a3m", os.path.join(af2_dir, f"MSA4af2_{count}_{pdb}_{tag}.a3m"))
    elif co_msa_mode in ['unpaired_ali']:
        shutil.copy(f"output/MSAforAlphafold_all_unpaired.a3m", os.path.join(af2_dir, f"MSA4af2_{count}_{pdb}_{tag}.a3m"))
    elif co_msa_mode in ['single_pep']:
        # Add gaps in the mixed MSA for all sequences in the region of the ligand
        ligand_start, ligand_stop = l_applied_delimitations[-1]
        fillGapLigand(ligand_start, ligand_stop, "output/MSAforAlphafold_mixed.a3m", "output/MSAforAlphafold_single.a3m")
        shutil.copy("output/MSAforAlphafold_single.a3m", os.path.join(af2_dir, f"MSA4af2_{count}_{pdb}_{tag}.a3m"))

    fd, tmp_runcolab = tempfile.mkstemp()
    os.close(fd)

    change_singu_colab_cmd(f"{os.path.join(SCRIPTS_DIR, 'run_colab_template.sh')}", tmp_runcolab, COLABFOLD_CMD, SINGULARITY_CMD)

    shutil.copy(tmp_runcolab, os.path.join(af2_dir, 'run_colab.sh'))
    os.remove(tmp_runcolab)

if __name__ == "__main__":
    usage = " 4_PrepareDirectoryforAF2.py -i <protocol name> " \
            "[-c <name of the cmd_file listing all the calls to run af2>]\n"

    parser = argparse.ArgumentParser(usage=usage)

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
    parser.add_argument(
        '-i',
        '--input_protocol',
        type=str,
        default='all',
        choices=list_of_possible_protocols + ['all'],
        help="Provide the name of the protocols to be prepared. If no argument, all protocols will be prepared",
    )
    parser.add_argument(
        '-c',
        '--cmd_af2run',
        type=str,
        default='cmd_global_af2runs.sh',
        help="Provide the name of the file in which the command lines to run AF2 are stored.",
    )
    options = parser.parse_args()
    if options.input_protocol == 'all':
        l_protocols = list_of_possible_protocols
    else:
        l_protocols = [options.input_protocol]

    fcmd_global_runaf2 = os.path.join(AF2RUNS_DIR, options.cmd_af2run)
    fout_cmd_global_runaf2 = open(fcmd_global_runaf2, "w")

    for protocol in l_protocols:

        print("Preparing concatenated alignment for protocol:",protocol)
        d_protocol = assignProtocolInformation(protocol)

        ## Variables to change for a new run ##
        run_index = 1

        # Option to use to run on a set of PDB only
        is_restrict_run = False
        set_restricted = set(['6A30'])

        working_dir = WORKING_DIR
        os.chdir(working_dir)

        if run_index == 1:
            # Additional runs can be added by incrementing the run_index
            alis_dir = DATA_DIR
            pdblist_formatted = pdblist_formatted
            count = 0
            count_ini = 0

        fileinput = os.path.join(alis_dir, pdblist_formatted)
        tag = d_protocol['tag']
        ali_mode = d_protocol['ali']
        rec_ext = d_protocol['rec_ext']
        lig_ext = d_protocol['lig_ext']

        with open(fileinput,'r') as fin:
            #parse File_Listing_Uniprot_Inputs
            for ll in fin:
                do_run = True
                stab = ll.split("\t")
                s = ll.split()
                if ll[0] != "#":
                    uniprot = stab[1].split(':')[1]
                    start = stab[2].split(':')[1].replace('\n', '')
                    stop = stab[3].split(':')[1].replace('\n', '')
                    l_start.append(start)
                    l_stop.append(stop)
                    l_uniprot.append(uniprot)
                else:
                    # Because of the structure of the Input file, we first get the data and only in the second pass we enter the 'if'
                    if count > count_ini:
                        # It is not the first entry so we can run the previous pdb
                        if is_restrict_run:
                            if not pdb in set_restricted:
                                do_run = False
                        if do_run:
                            createCmdColab(SCRIPTS_DIR, FASTAMSA_DIR,
                                           count, pdb, tag, l_uniprot, l_start, l_stop,
                                           co_msa_mode=ali_mode, rec_extension=rec_ext, lig_extension=lig_ext)
                            fullpath_af2_dir = os.path.join(AF2RUNS_DIR, f"{count}_{pdb}", f"af2_{count}_{pdb}_{tag}")
                            a3m_fname = f"MSA4af2_{count}_{pdb}_{tag}.a3m"
                            fout_cmd_global_runaf2.write(f"cd {fullpath_af2_dir} ; sh run_colab.sh {a3m_fname}\n")
                    count += 1
                    pdb = s[2]
                    print(f"STEP {count} Running {pdb} ")
                    l_uniprot = []
                    l_stop = []
                    l_start = []
                    wpath_dir = os.path.join(AF2RUNS_DIR, f"{count}_{pdb}")
                    if not os.path.isdir(wpath_dir):
                        os.mkdir(wpath_dir)
                    os.chdir(wpath_dir)
                    print(f"Working in {wpath_dir}")
            if is_restrict_run:
                if not pdb in set_restricted:
                    do_run = False
            if do_run:
                createCmdColab(SCRIPTS_DIR, FASTAMSA_DIR,
                               count, pdb, tag, l_uniprot, l_start, l_stop,
                               co_msa_mode=ali_mode, rec_extension=rec_ext, lig_extension=lig_ext)
                fullpath_af2_dir = os.path.join(AF2RUNS_DIR, f"{count}_{pdb}", f"af2_{count}_{pdb}_{tag}")
                a3m_fname = f"MSA4af2_{count}_{pdb}_{tag}.a3m"
                fout_cmd_global_runaf2.write(f"cd {fullpath_af2_dir} ; sh run_colab.sh {a3m_fname}\n")
    fout_cmd_global_runaf2.close()
