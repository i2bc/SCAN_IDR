import os, sys
import argparse
import glob
import re

def create_qid_cov_lines(number_of_chains, number_of_terms=3, default_qid='20.', default_cov='50.'):
    str_qid_cov = ""
    for ii in range(number_of_chains):
        str_qid_cov += "QID{index}={qid}\nCOV{index}={cov}\nDELIM{index}='_'\n".format(index=ii+1,
                                                                                     qid=default_qid,
                                                                                     cov=default_cov)
        str_qid_cov += "\n"
    str_qid_cov += "\n"
    for jj in range(number_of_terms):
        str_qid_cov += "TERM{index}='>'\n".format(index=jj+1)
    str_qid_cov += "\nIT"
    return str_qid_cov

def create_root_fasta_lines(number_of_chains, fasta_list):
    str_new = ""
    for ii in range(number_of_chains):
        root_name, ext = os.path.splitext(os.path.basename(fasta_list[ii]))
        str_new += "ROOT_NAME{index}={rootname}\n".format(index=ii+1,rootname=root_name)
    str_new += "\n"
    for ii in range(number_of_chains):
        str_new += "PROTEIN{index}=${{ROOT_NAME{index}}}.fasta\n".format(index=ii+1)
    stoichio = ','.join(["1" for ii in range(number_of_chains)])
    str_new += "\nSTOICHIOMETRY={}\n".format(stoichio)

    str_new += "\nDO"
    return str_new

def create_runalpha_line(number_of_chains):
    """
    Change the line:
    "${{PYTHON}} ${{ALPHASCRIPTS}}/GetAlphaFoldInputs.py -msa 0,1
    --delimitation ${delim_list[0]},${delim_list[1]} -o mixed --jsonname json_step1 --jsonoutput json_step2
    --colab -st ${STOICHIOMETRY}"
    """
    str_new = "${PYTHON} ${ALPHASCRIPTS}/GetAlphaFoldInputs.py -msa "
    str_msa = ",".join([str(ii) for ii in range(number_of_chains)])
    str_new += str_msa
    str_new += " --delimitation "
    str_delim = ",".join(["${{delim_list[{}]}}".format(ii) for ii in range(number_of_chains)])
    str_new += str_delim
    str_new += " --jsonname json_step1 --jsonoutput json_step2 --colab -st ${STOICHIOMETRY} -o "
    return str_new

def create_generate_line(number_of_chains):
    """
    Change the line:
    "${PYTHON} ${ALPHASCRIPTS}/GenerateMSA4EachPartner.py -i ${WORKING_DIRECTORY} -l ${PROTEIN1},${PROTEIN2} -m mmseqs_profile"

    """
    str_new = "${PYTHON} ${ALPHASCRIPTS}/GenerateMSA4EachPartner.py -i ${WORKING_DIRECTORY} -l "
    str_chains = ",".join(["${{PROTEIN{}}}".format(ii+1) for ii in range(number_of_chains)])
    str_new += str_chains
    str_new += " -m mmseqs_profile"

    return str_new

def create_list_lines(number_of_chains):

    str_new = ""
    str_new += "index_partners=({})\n".format(" ".join([str(x) for x in range(number_of_chains)]))
    str_new += "input_list=({protlist})\n".format(protlist=' '.join(["${{PROTEIN{}}}".format(ii) for ii in range(1,number_of_chains+1)]))
    str_new += "rootname_list=({rootlist})\n".format(rootlist=' '.join(["${{ROOT_NAME{}}}".format(ii) for ii in range(1, number_of_chains + 1)]))
    str_new += "qid_list=({qidlist})\n".format(qidlist=' '.join(["${{QID{}}}".format(ii) for ii in range(1,number_of_chains+1)]))
    str_new += "cov_list=({covlist})\n".format(covlist=' '.join(["${{COV{}}}".format(ii) for ii in range(1,number_of_chains+1)]))
    str_new += "delim_list=({delimlist})\n".format(delimlist=' '.join(["${{DELIM{}}}".format(ii) for ii in range(1,number_of_chains+1)]))
    str_new += "\n\nif"

    return str_new

def create_select_options(status_file):

    if status_file == "initialize":
        str_new = """DO_BUILD_MSA=true
DO_RUN_WITH_MMSEQS=true
DO_RUN_WITH_HHBLITS=false
DO_FILTER_ALI=true
DO_FETCH_FL_SEQS=false
DO_FILTER_MMSEQS_WITH_FULL_LENGTH_SEQS=true
DO_RUN_MAFFT=true
DO_MATCH_SPECIES=true
DO_PAIR_SEQ_DELIM=true

FULL_LENGTH_OPTION"""

    if status_file == "initialize_onlymsa":
        str_new = """DO_BUILD_MSA=true
DO_RUN_WITH_MMSEQS=true
DO_RUN_WITH_HHBLITS=false
DO_FILTER_ALI=true
DO_FETCH_FL_SEQS=false
DO_FILTER_MMSEQS_WITH_FULL_LENGTH_SEQS=true
DO_RUN_MAFFT=true
DO_MATCH_SPECIES=false
DO_PAIR_SEQ_DELIM=false

FULL_LENGTH_OPTION"""

    elif status_file == "final":
        str_new = """DO_BUILD_MSA=false
DO_RUN_WITH_MMSEQS=true
DO_RUN_WITH_HHBLITS=false
DO_FILTER_ALI=false
DO_FETCH_FL_SEQS=false
DO_FILTER_MMSEQS_WITH_FULL_LENGTH_SEQS=false
DO_RUN_MAFFT=false
DO_MATCH_SPECIES=true
DO_PAIR_SEQ_DELIM=true

FULL_LENGTH_OPTION"""

    return str_new

def create_singucall(running_modes, line=None):

    if line == "hhblits":
        str_new = "# 1 - Run HHblits\n"
        if running_modes == "i2bc":
            str_new += "\t${HHBLITS_BIN}/hhblits -i sequence_${i}.fasta -o sequence_${i}.hhr -oa3m sequence_${i}.a3m -n ${IT} -d ${HHBLITS_DB} -cpu ${NCPUS} -all -B ${MAX_ALI} -Z ${MAX_ALI} -maxseq ${MAX_SEQ_HHBLITS} -maxfilt ${MAX_SEQ_HHBLITS} -realign_max ${MAX_SEQ_HHBLITS} -maxres ${MAXRES}\n"
            str_new += "\n\tfi\n"
    if line == "mafft":
        str_new = "# 4 - Run mafft\n"
        if running_modes == "i2bc":
            str_new += "\t\tmafft --randomseed 1 --thread -1 sequence_${i}_filtered_FL_top${MAX_NSEQ_IN_MAFFT}.fasta > sequence_${i}_filtered_FL_mafali.fasta\n"
            str_new += "\n\t\tif"
    return str_new

def create_comsa(comsa_mode):


    if comsa_mode == "all":
        mixed, paired, allunpaired = ['true', 'true', 'true']
    elif comsa_mode == "none":
        mixed, paired, allunpaired = ['false', 'false', 'false']
    elif comsa_mode == "paired":
        mixed, paired, allunpaired = ['false', 'true', 'false']
    elif comsa_mode in ["mixed", "mixed_ali"]:
        mixed, paired, allunpaired = ['true', 'false', 'false']
    elif comsa_mode in ["allunpaired", "unpaired_ali"]:
        mixed, paired, allunpaired = ['false', 'false', 'true']
    elif comsa_mode == "single_pep": # We create a mixed and later in the PrepareDirforAF2 we replace pep ali by gaps.
        mixed, paired, allunpaired = ['true', 'false', 'false']
    else:
        print(f"Unknown option for requested comsa_mode {comsa_mode} called with -k option. Please reformulate")
        sys.exit()
    str_new = f"""RUN_MIXED_ALI_OPTION={mixed}
RUN_PAIRED_ALI_OPTION={paired}
RUN_ALLUNPAIRED_ALI_OPTION={allunpaired}

QID1"""

    return str_new

def create_scripts_line(ali_scripts_dir):
    str_new= f"SCRIPTS={ali_scripts_dir}\nALPHASCRIPTS"
    return str_new

def create_workingdir_line(workingdir):

    str_new = ""
    str_new += "WORKING_DIRECTORY={}\n".format(workingdir)

    return str_new

def create_cmd(input_dir_fasta, list_files='', name_filter='',
               output_cmd_file="./cmd_new.sh",
               status_file='',
               running_mode='',
               comsa_mode='',
               initial_qid='20.',
               initial_cov='50.',
               ali_scripts_dir='',
               working_dir=None):
    """

    Args:
        input_dir_fasta:
        list_files:
        name_filter:
        output_cmd_file:
        status_file:
        running_mode:

    Returns:

    """
    if len(list_files) == 0:
        pattern = "*{}*.fas*".format(name_filter)
        list_fasta_files = glob.glob(pattern)
    else:
        list_fasta_files = list_files.split(",")
    number_of_targets = len(list_fasta_files)
    os.chdir(input_dir_fasta)
    if not working_dir:
        working_dir = os.path.abspath(os.getcwd())

    print("Working dir is :", working_dir)

    with open(output_cmd_file, "w") as fout_cmd:
        with open(os.path.join(os.path.dirname(__file__), "cmd_template.sh"), "r") as fin_cmd:
            newcmd = fin_cmd.read()

            pat0 = re.compile("ROOT_.+?DO", re.DOTALL)
            str_fasta = create_root_fasta_lines(number_of_targets, list_fasta_files)
            newcmd = re.sub(pat0, str_fasta, newcmd)

            pat1 = re.compile("QID1=.+?IT", re.DOTALL)
            str_qid_cov = create_qid_cov_lines(number_of_targets, default_qid=initial_qid, default_cov=initial_cov)
            newcmd = re.sub(pat1, str_qid_cov, newcmd)

            pat2 = re.compile("index_partners.+?if", re.DOTALL)
            str_list_lines = create_list_lines(number_of_targets)
            newcmd = re.sub(pat2, str_list_lines, newcmd, count=1)

            pat3 = re.compile("..PYTHON. ..ALPHASCRIPTS./GetAlphaFoldInputs.+?-o ")
            str_alpha_line = create_runalpha_line(number_of_targets)
            newcmd = re.sub(pat3, str_alpha_line, newcmd)

            pat4 = re.compile("..PYTHON. ..ALPHASCRIPTS./GenerateMSA4EachPartner.+?-m mmseqs_profile")
            str_generate_line = create_generate_line(number_of_targets)
            newcmd = re.sub(pat4, str_generate_line, newcmd)

            pat5 = re.compile("WORKING_DIRECTORY=.+?\n")
            str_working_dir = create_workingdir_line(working_dir)
            if os.name =='nt':
                str_working_dir = '/'.join(str_working_dir.split('\\'))
            newcmd = re.sub(pat5, str_working_dir, newcmd)

            pat6 = re.compile("DO_BUILD_MSA=.+?FULL_LENGTH_OPTION", re.DOTALL)
            str_options_line = create_select_options(status_file)
            newcmd = re.sub(pat6, str_options_line, newcmd)

            pat7 = re.compile("# 4 - Run mafft.+?if", re.DOTALL)
            str_runmaff_line = create_singucall(running_mode, line='mafft')
            newcmd = re.sub(pat7, str_runmaff_line, newcmd)

            pat8 = re.compile("# 1 - Run HHblits.+?fi\n", re.DOTALL)
            str_runhhbl_line = create_singucall(running_mode, line='hhblits')
            newcmd = re.sub(pat8, str_runhhbl_line, newcmd)

            pat9 = re.compile("RUN_MIXED.+?QID1", re.DOTALL)
            str_comsa_line = create_comsa(comsa_mode)
            newcmd = re.sub(pat9, str_comsa_line, newcmd)

            pat10 = re.compile("SCRIPTS=.+?ALPHASCRIPTS", re.DOTALL)
            str_aliscripts_line = create_scripts_line(ali_scripts_dir)
            newcmd = re.sub(pat10, str_aliscripts_line, newcmd)

        fout_cmd.write(newcmd)
    fout_cmd.close()
    print(newcmd)


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i',
        '--input_dir_fasta',
        type=str,
        default="./",
        help="Working directory where all the .fas or .fasta input files are stored. "
             "Filters or list of files can be added to limit the set of files to consider."
    )
    parser.add_argument(
        '-f',
        '--file_name_filter',
        type=str,
        default="",
        help="Filter string used with wild-cards to recover all the fasta files to consider as input."
    )
    parser.add_argument(
        '-l',
        '--list_of_fasta_files',
        type=str,
        default="toto.fasta,titi.fasta,tata.fasta",
        help="List of files written in a comma separated string."
    )
    parser.add_argument(
        '-o',
        '--output_cmd_file',
        type=str,
        default="./cmd_new.sh",
        help="Filename of the file in which we list all the commands to run.",

    )
    parser.add_argument(
        '-s',
        '--status_file',
        type=str,
        default="initialize",
        choices=["initialize", "final", ],
        help="flag to indicate the status of the command file and the options to turn false or true depending on the prior generation of the alignment",
    )
    parser.add_argument(
        '-m',
        '--running_mode',
        type=str,
        default="i2bc",
        choices=["i2bc",  ],
        help="In which mode the script is going to be run.",
    )
    parser.add_argument(
        '-k',
        '--comsa_mode',
        type=str,
        default="all",
        choices=["none", "all", 'mixed', 'allunpaired', 'paired', 'mixed_ali', 'unpaired_ali', 'single_pep' ],
        help="What types of co-msa should be generated. 'mixed_ali'='mixed' and 'allunpaired'='unpaired_ali'",
    )
    parser.add_argument(
        '-w',
        '--working_dir',
        default=None,
        help="What is the directory where the msa and fasta are located and from where the cmd script will be run",
    )
    parser.add_argument(
        '-a',
        '--ali_scripts',
        type=str,
        default=".",
        help="Full path of the alignment scripts",
    )
    parser.add_argument(
        '-q',
        '--query_identity',
        type=str,
        default="25.",
        help="Generic identity to the query sequence threshold applied by default for all the inputs",
    )
    parser.add_argument(
        '-c',
        '--coverage',
        type=str,
        default="50.",
        help="Generic coverage threshold applied by default for all the inputs",
    )

    options = parser.parse_args()
    create_cmd(options.input_dir_fasta,
               list_files=options.list_of_fasta_files,
               name_filter=options.file_name_filter,
               output_cmd_file=options.output_cmd_file,
               status_file=options.status_file,
               running_mode=options.running_mode,
               comsa_mode=options.comsa_mode,
               initial_qid=options.query_identity,
               initial_cov=options.coverage,
               ali_scripts_dir=options.ali_scripts,
               working_dir=options.working_dir,
               )