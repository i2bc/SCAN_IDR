import os,re,sys
import configparser

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)

config = configparser.ConfigParser()
config.read(f'{script_dir}/config.ini')


WORKING_DIR = config['DEFAULT']['WORKING_DIR']
SCRIPTS_DIR = os.path.abspath(config['DEFAULT']['SCRIPTS_DIR'])
PYTHON_CMD = config['PROGRAMS']['PYTHON']

FASTAMSA_DIR = os.path.abspath(os.path.join(WORKING_DIR, 'data', 'fasta_msa'))
if not os.path.isdir(FASTAMSA_DIR):
    os.system(f"mkdir {FASTAMSA_DIR}")

pdblist_formatted = config['DEFAULT']['INPUT_PDB_LIST']

# This file created in data/ will be run to fetch the fasta sequences if required and prepare the command file which will generate the alignments
fout_filename = "cmd_get_fasta_msa.sh"
# This file will be created in data/ and will run the alignment generation process for every entries independently
cmd_file = os.path.join(FASTAMSA_DIR, "cmd_create_alignments.sh")

FLAG_RUN_GetFullSequence = True

fileinput = os.path.join(WORKING_DIR, 'data', pdblist_formatted)
fileoutput = os.path.join(WORKING_DIR, 'data', fout_filename)

l_uniprot = []
count = 0
count_entity = 0

with open(fileinput,'r') as fin:
    with open(fileoutput, 'w') as fout:
        for ll in fin:
            s = ll.split("\t")
            if ll[0] != "#":
                uniprot = s[1].split(':')[1]
                l_uniprot.append(uniprot.strip())
                count_entity += 1
                if FLAG_RUN_GetFullSequence:
                	fout.write(f"echo {count}_ Recovering fasta sequence for {uniprot}\n")
                	fout.write(f"{PYTHON_CMD} {os.path.join(SCRIPTS_DIR,'getFullSequence.py')} -i {uniprot} -d {FASTAMSA_DIR}\n")
            else:
                count += 1
        str_uniprots = ','.join(["{}.fasta".format(x) for x in l_uniprot])
        fout.write(f"echo Generating the command file to create all the individual alignments\n")
        fout.write(f"{PYTHON_CMD} {os.path.join(SCRIPTS_DIR,'createCmdFile.py')} -l {str_uniprots} -k none -a {SCRIPTS_DIR} -w {FASTAMSA_DIR} -o {cmd_file} -s initialize -m i2bc -q 25. -c 50.")
        fout.close()
os.system(f"sh {fileoutput}")
