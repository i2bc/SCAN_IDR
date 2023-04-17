

import sys,re,tempfile,os

#
# Tools to extract or rename one or several chains
#
def usage():
    print(""" PDBChain.py -i <pdbfile or pdb_id> -o <outputfile> -c <chain1> <chain2> ...
                [-n <newchainname1> <newchainname2> ... : Rename the chains ]
                [-noter : Add no TER line when the chain changes in the newfile (default : a TER is added) ]
                [-noend : Add no END line as last line of the newfile (default : a END is added) ] 
                [-renumber : renumber the outputfile from 1 to n, forcing sequential numbering throughout every chains]
                [-force_over_chain : renumber from 1 to n over the whole pdb irrespective of the chain]
                [-haddock : add the chain in the good case
                [-res <resBegin:resEnd> : select the residue and create a structure. WARNING: select only the first index residue corresponding to the bound !!
                [-noH : Remove hydrogens]
                [-noHETATM : Remove HETATM lines]
                If no chain in PDBfile, use "_" after -c 
                if whatever chain is concerned, use "*" after -c
    
    Example :
        python PDBChain.py -i 1GXD.pdb -o select.pdb -c A B -n A C -renumber
    
    Example for residue selection :
        python PDBChain.py -i 1GXD.pdb -o select.pdb -c A -n A -res 2:50 -renumber
        
    Example for haddock :
        python PDBChain.py -i haddock.pdb -o select.pdb -c A B -n A B -res 2:50 -renumber
        
    !! NB: If the pdbfile is not present in the directory,
             the pdb_id code can be given (under the -i option)
             and the program will automatically fetch the pdb from the net
    
                
                
    """)    
    sys.exit()
    
class PDBChain:
    """ PDBChain.PDBChain(pdb, out, dicvar) -> dump an output pdbfile
     
            # INPUT DESCRIPTION
            dicvar have key       : "chains"    -> list of chains to extract
                   may have key   : "newchains" -> list of new chain names
                   may have key   : "noter"     -> if not None, will switch off writing of TER 
                   may have key   : "noend"     -> if not None, will switch off writing of END 
                   may have key   : "renumber"  -> if not None, will renumber the outputfile from 1 to n, 
                                                        forcing sequential numbering throughout the chains
                   may have key   : "haddock"   -> True if not use.
                   may have key   : "res"       -> resBegin:resEnd Extract these residue you have to stipule chains and newchains with the same for instance A and A
    """
    def __init__(self, pdb, out, dicvar):

        self.pdb     = pdb
        self.out     = out
        self.dicvar  = dicvar

        self.dicpdb  = {}        # key : every chain, value : list of corresponding pdblines
        self.chain_rename = None # if rename : dict mapping old and new chain names        
        try:
            f            = open(self.pdb,"r")
        except IOError:
            print("ERROR: Couldn't find pdb "+self.pdb)

            dtmp = tempfile.mkdtemp()
            cmd = "wget 'https://files.rcsb.org/download/{}.pdb' -O {}"

            os.system(cmd.format(self.pdb[:4].upper(),dtmp+"/"+self.pdb[:4].upper()+".pdb"))
            f = open(os.path.join(dtmp,self.pdb[:4].upper()+".pdb"),"r")

        self.lines   = f.readlines()
        f.close()
        
        
        if 'noter' in self.dicvar and self.dicvar['noter']!=None:
            self.writeTER = False
            print("No TER will be added at the end of every chain")
        else:
            self.writeTER = True
            print("A TER will be added at the end of every chain")

        if 'noend' in self.dicvar and self.dicvar['noend']!=None:
            self.writeEND = False
            print("No END will be added at the end of the file")
        else:
            self.writeEND = True
            print("A END will be added at the end of the file")
        
        
        self.fillOptionPDBChain()
        
        if self.dicvar['chains']==["*"]:
            self.dicvar['chains'] = self.getChains()
            self.dicvar['newchains'] = [self.dicvar['newchains'][0]] * len(self.dicvar['chains'])
        
        if 'newchains' in self.dicvar and self.dicvar['newchains']!=None:
            if self.dicvar["res_select"]:
                self.dicvar['newchains'] =self.dicvar["chains"]
            self.chain_rename = dict(zip(self.dicvar['chains'],self.dicvar['newchains']))
            
            # dict mapping old and new chain name
            
        self.run()
    
    
    def fillOptionPDBChain(self):
        """ Fill default options in the dicvar """ 
        
        if "haddock" not in self.dicvar:
            self.dicvar["haddock"] = False
            self.chain_col = 21
        else :
            self.chain_col = 72
            
        if 'res_select' not in self.dicvar:
            self.dicvar["res_select"] = None

        if "noH" not in self.dicvar:
            self.dicvar['noH'] = False

        if "noHETATM" not in self.dicvar:
            self.dicvar['noHETATM'] = False
            
        return

    def getChains(self):
        chlist = set()
        for line in self.lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chlist.add(line[21])
        return(list(chlist))

    def run(self):
        if self.dicvar["res_select"]:
            self.resFilterPDBChain()
        else:
            self.extract_chains()
        self.save_output()
        if self.dicvar["renumber"]:
            self.renumber()

    def resFilterPDBChain(self):
        """only keep residues that are within the limits defined by -res option and change chain name if need be"""
        for line in self.lines:
            if (line.strip()) and (line.startswith("ATOM") | line.startswith("HETATM")):
                if line.startswith("ATOM") and self.dicvar["noH"] and line[13] == "H":
                    continue
                if line.startswith("HETATM") and self.dicvar["noHETATM"]:
                    continue
                if ("OXT" not in line):
                    index        = int(line[22:26].strip())
                    if (index > self.dicvar["res_select"][1]): # we select just the first res respecting the conditions
                        break
                    if (index >= self.dicvar["res_select"][0]) and (index <= self.dicvar["res_select"][1]): # if index is within bounds, keep it and change chain name if need be
                            newline = line[:self.chain_col] + self.dicvar["newchains"][0] + line[self.chain_col+1:]
                            try:
                                self.dicpdb[self.dicvar["newchains"][0]].append(newline)
                            except KeyError:
                                self.dicpdb[self.dicvar["newchains"][0]] = []
                                self.dicpdb[self.dicvar["newchains"][0]].append(newline)   

            
    def extract_chains(self):
        """select chains of interest """
        for line in self.lines:
            if (line.strip()) and (line.startswith("ATOM") | line.startswith("HETATM")):
                if line.startswith("ATOM") and self.dicvar["noH"] and line[13] == "H":
                    continue
                if line.startswith("HETATM") and self.dicvar["noHETATM"]:
                    continue
                if (not "OXT" in line):
                    chain = line[self.chain_col]
                    if self.chain_rename:
                        if chain in self.chain_rename:
                            newline = line[:self.chain_col] + self.chain_rename[chain] + line[self.chain_col+1:]
                        else:
                            newline = line 
                    else:
                        newline = line 
                    try:
                        self.dicpdb[chain].append(newline)
                    except KeyError:
                        self.dicpdb[chain] = []
                        self.dicpdb[chain].append(newline)   
                                           
    def renumber(self):
        import RenumberPDB
        # firstres should be written as <chain1>:1 or <newchain1>:1
        for chain in self.dicvar["chains"]:
            if self.chain_rename:
                firstres = "%s:1" %(self.chain_rename[chain])
            else:
                firstres = "%s:1" %(chain)
            if not self.dicvar["force_over_chain"]:
                # Start every chain by one
                renumber_object = RenumberPDB.RenumberPDB(self.out,firstres,self.out,force_sequential = True, force_sequential_overchains = False)
            else:
                # Increment residue index over all the chains
                renumber_object = RenumberPDB.RenumberPDB(self.out,firstres,self.out,force_sequential = True, force_sequential_overchains = True)
            del(renumber_object)
        
    def save_output(self):
        fout = open(self.out,"w")
        for chain in self.dicvar['chains']:
            fout.writelines(self.dicpdb[chain])
            if self.writeTER:
                fout.write("TER\n")
        if self.writeEND:
            fout.write("END\n")
        fout.close()
            
if __name__ == '__main__':
    
    DEBUG = None
    if len(sys.argv)==1 and DEBUG == None:
        usage()
    elif DEBUG == None:
        pdb = None
        dicvar = {}

    print(sys.argv)
    # ARG1    
    try:
        pdb = sys.argv[sys.argv.index('-i')+1]
    except:
        usage()
    # ARG2
    try:
        out = sys.argv[sys.argv.index('-o')+1]
    except:
        usage()
    # ARG3      
    try:
        dicvar['chains'] = []
        for chain in sys.argv[sys.argv.index('-c')+1:]:
            if chain[0]=="-":
                break
            if chain == "_":
                dicvar['chains'].append(" ")
            else:
                dicvar['chains'].append(chain)
    except:
        usage()
    # ARG4 : rename ?      
    try:
        dicvar['newchains'] = []
        for chain in sys.argv[sys.argv.index('-n')+1:]:
            if chain[0]=="-":
                break
            dicvar['newchains'].append(chain)
    except:
        dicvar['newchains'] = None
    # ARG5
    try:
        dicvar["noter"] = sys.argv[sys.argv.index('-noter')]
    except:
        dicvar["noter"] = None
        print("Option noter was set to None")
    # ARG5b
    try:
        dicvar["noend"] = sys.argv[sys.argv.index('-noend')]
    except:
        dicvar["noend"] = None
        print("Option noend was set to None")
    # ARG6
    try:
        dicvar["renumber"] = sys.argv[sys.argv.index('-renumber')]
    except:
        dicvar["renumber"] = None

    try:
        dicvar["res_select"] = [int(i) for i in sys.argv[sys.argv.index('-res')+1].split(":")]
        #
        # TODO: several selection management of instance 1:3  50:70. And find a solution for multichain without chain annotation !
        #
    except:
        pass

    try:
        test = sys.argv[sys.argv.index('-noH')]
        dicvar["noH"] = True
    except:
        pass

    try:
        test = sys.argv[sys.argv.index('-noHETATM')]
        dicvar["noHETATM"] = True
    except:
        pass

    try:
        test = sys.argv[sys.argv.index('-haddock')]
        dicvar["haddock"] = True
    except:
        pass
        
    try:
        test = sys.argv[sys.argv.index('-force_over_chain')]
        dicvar["force_over_chain"] = True
    except:
        dicvar["force_over_chain"] = False


    P = PDBChain(pdb,out,dicvar)
