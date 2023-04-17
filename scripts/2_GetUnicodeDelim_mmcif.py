import os,re,sys
import glob
from gemmi import cif
import configparser

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)

config = configparser.ConfigParser()
config.read(f'{script_dir}/config.ini')

WORKING_DIR = config['DEFAULT']['WORKING_DIR']
pdblist_formatted = config['DEFAULT']['INPUT_PDB_LIST']

mmcif_path = os.path.join(WORKING_DIR, 'data', 'cif', '*.cif*')
pat = re.compile('(\w\w\w\w)\.cif')
list_files = glob.glob(mmcif_path)
#
list_files.sort()
#
fout = open(os.path.join(WORKING_DIR, 'data', pdblist_formatted), "w")
l_full_data = []

for ii,cif_file in enumerate(list_files):
    # Parse every pdb file in cif format
    doc = cif.read(cif_file)
    block = doc.sole_block()
    m = pat.search(cif_file)
    if m:
        pdbcode = m.group(1)
    else:
        print("Problem, could not parse the pdb code ! Please check {}".format(cif_file))
    print("# Running with {}".format(cif_file))
    fout.write("# PDBCODE: {}\n".format(pdbcode))
    category = block.find_mmcif_category('_struct_ref_seq.')
    #print(list(category.tags))
    for row in category:
        # Parse every chain
        d_data = dict(zip(category.tags,row))
        uniprot_id = d_data['_struct_ref_seq.pdbx_db_accession']
        index_chain = d_data['_struct_ref_seq.align_id']
        beg_chain_in_uniprot = d_data['_struct_ref_seq.db_align_beg']
        end_chain_in_uniprot = d_data['_struct_ref_seq.db_align_end']
        pdb_chain = d_data['_struct_ref_seq.pdbx_strand_id']
        if len(uniprot_id) < 6:
            uniprot_id = '???'
        l_info = [pdbcode, uniprot_id, beg_chain_in_uniprot, end_chain_in_uniprot]
        if not l_info in l_full_data:
            l_full_data.append(l_info)
            print(pdbcode,index_chain,pdb_chain,uniprot_id,beg_chain_in_uniprot,end_chain_in_uniprot)
            fout.write("CHAIN:{}\tUNIPROT:{}\tSTART:{}\tSTOP:{}\n".format(pdb_chain,uniprot_id,beg_chain_in_uniprot,end_chain_in_uniprot))
fout.close()

