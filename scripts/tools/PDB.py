
import sys
import os

"""

AA3to1 [XXX] -> X        # Give one letter amino acids from 3 letters notation

class PDB(adress,analyse=True)
  -> analyse = False -> access to data_pdb
  -> analyse = True  -> sort of the data

Methods :
        # self.writePDB(chain,namefile)

Variables:
        # data_pdb ["order"]                           ->       List of chain
        #          [chain] -> ["order"]                ->       List of index residue in chain
        #                     [index]    ->  ["aa"]    ->       Amino acide name (3letters)
        #                                    ["coord"] ->       Coordinnate [x,y,z] (all are float)
        #                                    ["all"]   ->       The line corresponding of the residue number (PDB format)
        #                  -> ["all"]                  ->       List of all lines of the PDB file by chain
        #

        # nb_chain                      -> Number of chain
        # chain                         -> List of chain
        # nb_base        [chain]        -> Number of residue in the chain
        # index          [chain]        -> List index in the chain
        # sequence       [chain]        -> List of amino acids
        # index2sequence [chain][index] -> give the amino acid from the index number
        # index2coord    [chain][index] -> give the coordinate dictionnar from the index number
        #                               ->  [atom] -> [x,y,z] give the coordinate in a list

"""


AA_3to1 = {"ALA" :"A",\
           "CYS":"C",\
           "ASP":"D",\
           "GLU":"E",\
           "PHE":"F",\
           "GLY":"G",\
           "HIS":"H",\
           "ILE":"I",\
           "LYS":"K",\
           "LEU":"L",\
           "MET":"M",\
           "ASN":"N",\
           "PRO":"P",\
           "GLN":"Q",\
           "ARG":"R",\
           "SER":"S",\
           "THR":"T",\
           "VAL":"V",\
           "TRP":"W",\
           "TYR":"Y",\
           "HOH":"O"\
}
DNA_2to1 = {"DT" :"T",\
           "DA":"A",\
           "DG":"G",\
           "DC":"C",\
           }

class PDB:
    def __init__(self,adress,analyse=True,ignore_non_std=False,ignore_hydro=True):
        """ """
        # adress -> path of the pdb
        self.adress = adress                                # adress of pdb

        self.name   = os.path.basename(adress).split(".")[0]

        # Read the pdb
        # data_pdb ["order"]                           ->       List of chain
        #          [chain] -> ["order"]                ->       List of index residue in chain
        #                     [index]    ->  ["aa"]    ->       Amino acide name (3letters)
        #                                    ["coord"] ->       Coordinnate [x,y,z] (all are float)
        #                  -> ["all"]                  ->       List of each line of the PDB corresponding to the stipuled chain
        self.data_pdb = self.readPDB(adress,ignore_hydro)                # data in the pdb

        # Analyse pdb
        # nb_chain                      -> Number of chain
        # chain                         -> List of chain
        # nb_base        [chain]        -> Number of residue in the chain
        # index          [chain]        -> List index in the chain
        # sequence       [chain]        -> List of amino acids
        # index2sequence [chain][index] -> give the amino acid from the index number
        # index2coord    [chain][index] -> give the coordinate dictionnar from the index number
        # ignore_non_std : False by default, if True, all non-defined residues and nucleotides are ignored
        # ignore_hydro: True by default, if True, ignore all hydrogen atoms

        if analyse:
            self.nb_chain,self.chain,self.nb_base,self.index,self.sequence,self.index2sequence,self.index2coord = self.analysePDB(self.data_pdb,ignore_non_std)

        return


    def readPDB(self,adress,ignore_hydro):
        """Compute coords and sequence of protein"""

        ######## Ouverture du pdb
        pdb  = open(adress)

        data_pdb          = {}
        data_pdb["order"] = []

        for line in pdb:
            # If no 4 caracters we do not process the ATOM test
            #if len(line)<4:
            #    break

            # ATOM lines !
            if line[0:4] == "ATOM": # faster than re library

                # Scan information
                index        = line[22:26].strip()
                atom         = line[12:16].strip()
                res          = line[17:20].strip()
                coord        = [float(line[29:38]), float(line[38:46]), float(line[46:54])]
                chain        = line[21]

                if ignore_hydro:
                    if atom.strip().startswith("H"):
                        continue
                        
                # Stockage information
                if not data_pdb.has_key(chain):
                    data_pdb[chain] = {}
                    data_pdb["order"].append(chain)
                    data_pdb[chain]["order"] = []
                    data_pdb[chain]["all"]   = []
                if not data_pdb[chain].has_key(index):
                    data_pdb[chain][index]          = {}
                    data_pdb[chain][index]["aa"]    = res
                    data_pdb[chain][index]["coord"] = {}
                    data_pdb[chain][index]["all"]   = ""
                    data_pdb[chain]["order"].append(index)

                data_pdb[chain][index]["coord"][atom] = coord
                data_pdb[chain][index]["all"] += line
                data_pdb[chain]["all"].append(line)
        # Closing file
        pdb.close()

        return data_pdb

    def analysePDB(self,data_pdb,ignore_non_std=False):
        """ Analyse the feature of the pdb
        ignore_non_std: default = False, to ignore any non-standard residues or nucleotides
        """

        nb_chain       = len(data_pdb["order"])      # Number of chain
        chain          = data_pdb["order"]           # List of chains
        nb_base        = {}                          # Number of base in each chain
        index          = {}                          # Dictionnar of list of index
        new_index      = {}                          # same as index but ignores residues or nucleotides that aren't in AA_3to1 or DNA_2to1
        sequence       = {}                          # Dictionnar of list of residue
        index2sequence = {}                          # index -> AA
        index2coord    = {}                          # index -> coord

        for chain_current in data_pdb["order"]:
            nb_base[chain_current]         = len(data_pdb[chain_current]["order"])                    # Count the number of base in each chain
            index[chain_current]           = data_pdb[chain_current]["order"]                         # Get the list of index
            new_index[chain_current]       = []
            temp_aa                        = []                                                       # Temporar list for the sequence AA

            if not index2coord.has_key(chain_current):
                index2coord[chain_current] = {}

            for index_current in data_pdb[chain_current]["order"]:                                    # through the index residue
                if AA_3to1.has_key(data_pdb[chain_current][index_current]["aa"]):#This is not the case when the chain considered is RNA or DNA or when AA aren't identified (i.e. UKN)
                    temp_aa.append(AA_3to1[data_pdb[chain_current][index_current]["aa"]])
                    index2coord[chain_current][index_current] = data_pdb[chain_current][index_current]["coord"]          # Get coordinate
                    new_index[chain_current].append(index_current)
                if DNA_2to1.has_key(data_pdb[chain_current][index_current]["aa"]):#This is for DNA
                    temp_aa.append(DNA_2to1[data_pdb[chain_current][index_current]["aa"]])
                    index2coord[chain_current][index_current] = data_pdb[chain_current][index_current]["coord"]          # Get coordinate
                    new_index[chain_current].append(index_current)
            sequence[chain_current] = "".join(temp_aa)                                                # Get the sequence
            if ignore_non_std:
                index2sequence[chain_current]  = dict(zip(new_index[chain_current],sequence[chain_current]))
            else:
                index2sequence[chain_current]  = dict(zip(index[chain_current],sequence[chain_current]))

        if ignore_non_std:
            return nb_chain,chain,nb_base,new_index,sequence,index2sequence,index2coord
        else:
            return nb_chain,chain,nb_base,index,sequence,index2sequence,index2coord


    def writePDB(self,chain,file):
        """ Write a PDB corresponding to a chain"""
        if chain not in self.chain:
            print("No Chain %s in %"%(chain,self.adress))
            return

        f = open(file,"w")
        for line in self.data_pdb[chain]["all"]:
            f.write(line)

        f.close()
        return

