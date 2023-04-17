#
# A set of tools for PDB manipulations
# 

import os
import re
import tempfile
import shutil
import configparser

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)

config = configparser.ConfigParser()
config.read(f'{script_dir}/config.ini')

SCRIPTS_DIR = os.path.abspath(config['DEFAULT']['SCRIPTS_DIR'])

CMD_RENUMBER = f"python {os.path.join(SCRIPTS_DIR, 'tools', 'RenumberPDB.py')} "
CMD_PDBCHAIN = f"python {os.path.join(SCRIPTS_DIR, 'tools', 'PDBChain.py')} "

def renumberPDB(pdbin, chain, index_res1, pdbout=None):
    """
    Renumber the pdb, better using a single chain
    If pdbout not defined, pdbin is replaced
    """
    fd, tmpout = tempfile.mkstemp()
    os.close(fd)
    
    cmd = CMD_RENUMBER
    cmd += f" -s {pdbin} "
    cmd += f" -r1 {chain}:{index_res1} "
    cmd += f" -o {tmpout} "
    cmd_status = os.system(cmd)
    #print(cmd)
    if pdbout:
        os.system(f"mv {tmpout} {pdbout}")
    else:
        os.system(f"mv {tmpout} {pdbin}")
    return cmd_status    

def cutPDB(pdbin, start, stop, chain, pdbout=None, is_noter=True, is_noend=True, log=False):
    """
    Extract a segment between residue index start and stop. 
    Better if the input is monochain
    If no pdbout, pdbin is replaced by the output
    """
    cmd = CMD_PDBCHAIN
    fd, tmpout = tempfile.mkstemp()
    os.close(fd)
    cmd += f" -i {pdbin} "
    cmd += f" -o {tmpout} "
    cmd += f" -c {chain} "
    cmd += f" -n {chain} "
    cmd += f" -res {start}:{stop} "
    if is_noter:
        cmd += " -noter "
    if is_noend:
        cmd += " -noend "
    if not log:
        cmd += " > /dev/null "
    cmd_status = os.system(cmd)
    if pdbout:
        os.system(f"mv {tmpout} {pdbout}")
    else:
        os.system(f"mv {tmpout} {pdbin}")
    return cmd_status

def alterchainPDB(pdbin, current_chain, new_chain, pdbout = None, is_noter=True, is_noend=True, log=None):
    """
    Change in pdbin the chain name. If not pdbout, replace pdbin 
    """
    cmd = CMD_PDBCHAIN
    fd, tmpout = tempfile.mkstemp()
    os.close(fd)
    cmd += f" -i {pdbin} "
    cmd += f" -o {tmpout} "
    cmd += f" -c {current_chain} "
    cmd += f" -n {new_chain} "
    if is_noter:
        cmd += " -noter "
    if is_noend:
        cmd += " -noend "
    if not log:
        cmd += " > /dev/null "
    
    cmd_status = os.system(cmd)
    #print(cmd)
    if pdbout:
        os.system(f"mv {tmpout} {pdbout}")
    else:
        os.system(f"mv {tmpout} {pdbin}")
    return cmd_status

def fusePDB(l_pdb, fused_filename=None, add_ter=False, add_end=False):
    """
    Generates a tempfile pdb concatenating the pdbin l_pdb
    Returns the name of the tempfile
    """
    fd, foutname = tempfile.mkstemp()
    os.close(fd)
    fout = open(foutname,"w")
    for pdbin in l_pdb:
        fin = open(pdbin, "r")
        fout.write(fin.read())
        if add_ter:
            fout.write("TER\n")
    if add_end:
        fout.write("END\n")
    fout.close()
    if not fused_filename:
        return foutname
    else:
        os.system(f"mv {foutname} {fused_filename}")
        return ""


def splitPDB(pdbin, l_chains, is_noter=True, is_noend=True, log=False):
    """
    Input : a multichain pdb and a list of chains to process
    Returns a list of pdb in which every chain was isolated
    """
    l_pdbout = []
    for chain in l_chains:
        cmd = CMD_PDBCHAIN
        fd, pdbout = tempfile.mkstemp()
        os.close(fd)
        cmd += f" -i {pdbin} "
        cmd += f" -o {pdbout} "
        cmd += f" -c {chain} "
        cmd += f" -n {chain} "
        if is_noter:
            cmd += " -noter "
        if is_noend:
            cmd += " -noend "
        if not log:
            cmd += " > /dev/null "
        #print(cmd)
        cmd_status = os.system(cmd)
        l_pdbout.append(pdbout)
    return l_pdbout,cmd_status
    
def getLastResidueIndex(pdbin):
    """
    Scan the 10 last lines of the pdb and try to find the last residue index
    """
    with open(pdbin, 'r') as f:
        split_file = f.readlines()
        for ii in range(1,10):
            last_line = split_file[-ii]
            if len(last_line) > 0 and last_line.split()[0] == 'ATOM':
                residue_index = last_line[22:26].strip()
                break
    return residue_index

