
import sys
import PDB
import math
import pmath

def usage():
    print("""Contact4CAPRI.py -pdb <pdb> -t <threshold>  
                            [-hc <file_out>         : write the hydrophobic contacts (See the header for the features)
                            [-cc <file_out>         : write the charged contacts (See the header for the features)
                            [-cp <file_out>         : write the cation pi  contacts (See the header for the features)
                            [-all                   : screen the full contact information
                            [-co <contact_type>     : Type of the contact AA-> All atoms, BB-> only BackBone, CA,CB, SC (default: AA)
                            [-c <chain1:chain2>     : compute the contact between chain1 and chain2 (default : the first and second in the pdb file) 
                            [-inter                 : Give residue on each interface
                           
                python Contact.py -pdb 1GXD_complex.pdb -t 5
                python Contact.py -pdb 1GXD_complex.pdb -t 4 -hc hydrophobic.contact
                
                
                threshold advices :
                    hydrophobic  4A
                    charged      6A or 3.5 for salt bridge
                    cationPi     3.5
                    
                
                
            @author: Guerois'Team LBSR CEA
            """)
    sys.exit(1)

class Contact:
    """
    Contact(dicvar)       have key    "pdb"           : adress of the pdb file
                          have key    "threshold"     : threshold for considering a contact
                          may key     "interface"     : True. Search the interface residue
    
    
    Avalaible variables :
    
        list_contact   ["contact"]   : list of contact A->B
                       ["interA"]    : list of residue of the interface in the partner A
                       ["interB"]    : list of residue of the interface in the partner B
                       
        data_pdb  [chain]->[index]->[x,y,z]   : pdb data
    
    """
    
    
    
    
    def __init__(self,dicvar):
        
        self.dicvar = dicvar
        
        self.fillOptionContact()

        self.f =lambda x,y:x-y
        
        
        
        
        if not "all" in self.dicvar.keys():
            self.dicvar["all"] = None


        # Interface and contact type
        self.readPdbContact()
        

        if self.dicvar["capri"]:
            self.computeCAPRIContact(self.forCAPRI)
            return
        
        
        
        # Contact with all informations
        if self.dicvar["all"]:
            self.pdb = PDB.PDB(self.dicvar["pdb"])
            # Give all informations about contact
            self.computeInformationContact(self.pdb)
        
            return

        # Search hydrophobic contact
        if self.dicvar["hc"]:
            self.pdb = PDB.PDB(self.dicvar["pdb"])
            # Give all informations about contact
            self.printer = self.seachHydrophobicContact(self.pdb,self.dicvar["hc"])
            return

        
        # Search Charged Interaction
        if self.dicvar["cc"]:
            self.pdb = PDB.PDB(self.dicvar["pdb"])
            # Give all informations about contact
            self.printer = self.seachChargedContact(self.pdb,self.dicvar["cc"])
            return
            
        # Search Cathion Pi Interaction
        if self.dicvar["cp"]:
            self.pdb = PDB.PDB(self.dicvar["pdb"])
            self.printer = self.searchCathionPiContact(self.pdb,self.dicvar["cp"])
            return
        
        # Interface and contact type
        #self.readPdbContact()
        

        #if self.dicvar["capri"]:
        #    self.computeCAPRIContact(self.forCAPRI)
        #    return
        
        
        
        self.computeContact(self.dicvar["type"])
        
        
        if self.dicvar["interface"]:
            self.computeInterfaceContact()
        if self.dicvar["dump"] and self.dicvar["interface"]:
            self.screenInterfaceContact()
        return

    
    def fillOptionContact(self):
        """ Fill Default Option  """

        if not "dump" in self.dicvar.keys():
            self.dicvar["dump"] = False
            self.contact_list = {"contact":[],"interA":[],"interB":[]}
        
        if not "interface" in self.dicvar.keys():
            self.dicvar["interface"] = False
        else:
            self.contact_list = {"contact":[],"interA":[],"interB":[]}



        if not "capri" in self.dicvar.keys():
            self.dicvar["capri"] = False
        else:
            return

        
        self.dicvar["threshold"] *= self.dicvar["threshold"] 
        
        if self.dicvar["threshold"]>400.0:
            self.dicvar["threshold_limit"] = 99999.0
        else:
            self.dicvar["threshold_limit"] = 400.0
        
        
        if not "type" in self.dicvar.keys():
            self.dicvar["type"] = self.isContactAA
        elif self.dicvar["type"] == "AA":
            self.dicvar["type"] = self.isContactAA
        elif self.dicvar["type"] == "CA":
            self.dicvar["type"] = self.isContactCA
        elif self.dicvar["type"] == "BB":
            self.dicvar["type"] = self.isContactBB
        elif self.dicvar["type"] == "SC":
            self.dicvar["type"] = self.isContactSC
       
        if not "all" in self.dicvar.keys():
            self.dicvar["all"] = None
        
        if not "hc" in self.dicvar.keys():
            self.dicvar["hc"] = None
        if not "cc" in self.dicvar.keys():
            self.dicvar["cc"] = None
        
        if not "cp" in self.dicvar.keys():
            self.dicvar["cp"] = None
        
        return
    
    
    def readPdbContact(self):
        """Compute coords and sequence of protein"""

        ######## Ouverture du pdb
        pdb  = open(self.dicvar["pdb"])
        
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
                atom         = line[13:16].strip()
                res          = line[17:20].strip()
                coord        = [float(line[29:38]), float(line[38:46]), float(line[46:54])]
                chain        = line[21]
                
                # Stockage information
                if not chain in data_pdb.keys():
                    data_pdb[chain] = {}
                    data_pdb["order"].append(chain)
                    data_pdb[chain]["order"] = []
                    
                if not index in data_pdb[chain].keys():
                    data_pdb[chain][index]          = {}
                    data_pdb[chain][index]["aa"]    = res
                    data_pdb[chain][index]["coord"] = {}
                    data_pdb[chain]["order"].append(index)
                    
                data_pdb[chain][index]["coord"][atom] = coord
        
        self.data_pdb = data_pdb
        
        # Closing file
        pdb.close()
        
        return data_pdb
        
        
        
        
    def computeInformationContact(self,pdb):
        """ """
        
        c = "AB"
        
        # List Inter
        list_couple_chain =  [[c[chain1],c[chain2]] for chain1 in range(len(c)) for chain2 in range(len(c))   if chain1 < chain2]
        
        # List intra
        list_couple_chain = [[c[chain1],c[chain2]] for chain1 in range(len(c)) for chain2 in range(len(c)) if chain1 <= chain2]
        
        
        
        for chains in list_couple_chain:
            chain1,chain2 = chains
            
            # Chain1
            for index1 in pdb.index[chain1]:
                coord1 = pdb.index2coord[chain1][index1]
                
                # Versus chain2
                for index2 in pdb.index[chain2]:
                    coord2 = pdb.index2coord[chain2][index2]
                    
                    if (chain1 == chain2) and (index1 >= index2) : 
                        continue
                    
                    info_contact = self.infoResContact(coord1,coord2)   # [True/False,[[atom1,atom2],[atom1,atom2],...]]
                    
                    if info_contact[0]:
                        for atoms in info_contact[1]:
                            print(chain1, index1,pdb.index2sequence[chain1][index1], atoms[0] , "|", chain2,index2,pdb.index2sequence[chain2][index2],atoms[1])
            
        return
        
    
    
    def seachHydrophobicContact(self,pdb,adress):
        """ Search the hydrophobic contacts"""
        
        # List intra
        c = pdb.chain
        list_couple_chain = [[c[chain1],c[chain2]] for chain1 in range(len(c)) for chain2 in range(len(c)) if chain1 <= chain2]
        
        
        # AA contains hydrophobic chain (at least 2 CH2 in the SC)
        hydrophobic = "GMVILACPWYFKR"
        bb = ["N","CA","C","O","NZ","NE","NH1","NH2"]
        
        
        
        printer = ""
        printer += "# Hydrophobic Contacts\n"
        printer += "# %s\n"%(self.dicvar["pdb"])
        printer += "# Hydrophobic : Side Chain of  MVILACPWYF   (G -> just the CA  (K -> without NZ (R -> without NE,CZ,NH1,NH2\n"
        printer += "# Threshold : %s\n"%(math.sqrt(self.dicvar["threshold"]))
        
        # Chains loop
        for chain in list_couple_chain:
            chain1,chain2 = chain
            
            # index of fisrt chain
            for index1 in pdb.index[chain1]:
                
                
                coord1 = pdb.index2coord[chain1][index1]
                aa1    = pdb.index2sequence[chain1][index1]
                
                
                # index of second chain
                for index2 in pdb.index[chain2]:
                    
                    
                    # keep after the redundancies
                    if chain1 == chain2 and int(index1) >= int(index2):
                        continue
                    
                    
                    coord2 = pdb.index2coord[chain2][index2]
                    aa2    = pdb.index2sequence[chain2][index2]
                    
                    # Just residue with hydrophobic part interesting
                    if (aa1 not in hydrophobic) or (aa2 not in hydrophobic):
                        continue
                    
                    for atom1 in coord1:
                        for atom2 in coord2:
                            
                            # For the glycine
                            if atom1 not in bb or (aa1 == "G" and atom1 == "CA"):
                                if atom2 not in bb or (aa2 == "G" and atom2 == "CA"):
                                    
                                    
                                    # Not the CZ for R !! (CZ is not in bb list cause the aromatic residues have got one !)
                                    if (aa1 == "R" and atom1 == "CZ") or (aa2 == "R" and atom2 == "C2"):
                                        continue                                   
                                    
                                    d = self.getDistanceContact(coord1[atom1],coord2[atom2])
                                    
                                    if d < self.dicvar["threshold"]:
                                        printer += "%s%s%s_%s %s%s%s_%s %.2f\n"%(aa1,index1,chain1,atom1,aa2,index2,chain2,atom2,math.sqrt(d))
                                        
        # For lib use                   
        if adress != "X":                             
            f = open(adress,"w")                        
            f.write(printer)
            f.close()
            
        return printer




    def seachChargedContact(self,pdb,adress):
        """ Search the hydrophobic contacts"""
        
        # List intra
        c = pdb.chain
        list_couple_chain = [[c[chain1],c[chain2]] for chain1 in range(len(c)) for chain2 in range(len(c)) if chain1 <= chain2]
        
        
        # AA contains hydrophobic chain (at least 2 CH2 in the SC)
        positively   = "RKH"
        negatively   = "DE" 
        bb = ["N","CA","C","O","NZ","NE","NH1","NH2"]
        
        charged_atoms = {"R":["NH1","NH2","NE"],\
                         "K":["NZ"],\
                         "E":["OE1","OE2"],\
                         "D":["OD1","OD2"],\
                         "H":["ND1","NE2"]\
                         }
        
        printer = ""
        printer += "# Charged Contacts\n"
        printer += "# %s\n"%(self.dicvar["pdb"])
        printer += "# Residues charged positively R (NH1,NH2,NE), K (NZ), H (ND1,NE2)\n"
        printer += "# Residues charged negatively D (OD1,OD2), E (OE1,OE2)\n"
        printer += "# The distance written is the lowest\n"
        printer += "# Threshold : %s\n"%(math.sqrt(self.dicvar["threshold"]))
        
        # Chains loop
        for chain in list_couple_chain:
            chain1,chain2 = chain
            
            # index of fisrt chain
            for index1 in pdb.index[chain1]:
                
                
                coord1 = pdb.index2coord[chain1][index1]
                aa1    = pdb.index2sequence[chain1][index1]
                
                
                # index of second chain
                for index2 in pdb.index[chain2]:
                    
                    
                    # keep after the redundancies
                    if chain1 == chain2 and int(index1) >= int(index2):
                        continue
                    
                    
                    coord2 = pdb.index2coord[chain2][index2]
                    aa2    = pdb.index2sequence[chain2][index2]
                    
                    # Charged complementarity ?
                    if (aa1 in positively and aa2 in negatively) or (aa1 in negatively and aa2 in positively):
                        list_distance = []
                        for atom1 in charged_atoms[aa1]:
                            for atom2 in charged_atoms[aa2]:
                                
                                # Test if the coordinate of these atoms exist
                                try:
                                    d = self.getDistanceContact(coord1[atom1],coord2[atom2])
                                    if d < self.dicvar["threshold"]:
                                        list_distance.append(d)
                                        
                                        # One couple of atoms match -> do not need to continue on these residue
                                except:
                                    pass
                        if len(list_distance) >0:
                            d = min(list_distance)
                            printer += "%s%s%s_%s %s%s%s_%s %.2f\n"%(aa1,index1,chain1,atom1,aa2,index2,chain2,atom2,math.sqrt(d))          
                            
        # For lib use
        if adress != "X":
            f = open(adress,"w")                        
            f.write(printer)
            f.close()
            
        return printer


    def searchCathionPiContact(self,pdb,adress):
        """ Search the Cation PI """
        
        
        # List intra
        c = pdb.chain
        list_couple_chain = [[c[chain1],c[chain2]] for chain1 in range(len(c)) for chain2 in range(len(c)) if chain1 <= chain2]
        
        # AA contains hydrophobic chain (at least 2 CH2 in the SC)
        positively   = "RKH"
        aromatic     = "FYW"
        #bb = ["N","CA","C","O","NZ","NE","NH1","NH2"]
        
        charged_atoms = {"R":["NH1","NH2","NE"],\
                         "K":["NZ"],\
                         "H":["ND1","NE2"]\
                         }
        last_atom = {"F":"CG_CD1_CD2_CE1_CE2_CZ",\
                     "W":"CG_CD1_CD2_CE1_CE2_CZ",\
                     "Y":"CG_CD1_CD2_CE1_CE2_CZ"\
                     }
        
        
        printer = ""
        printer += "# Cation Pi Contacts\n"
        printer += "# %s\n"%(self.dicvar["pdb"])
        printer += "# Residues charged positively R (NH1,NH2,NE), K (NZ), H (ND1,NE2)\n"
        printer += "# Residues Aromatic F,Y,W -> Barycentre between the cycle\n"
        printer += "# If the distance is below %s between one atoms charged positively and the barycentre of an aromatic\n"%(math.sqrt(self.dicvar["threshold"]))
        
        # Chains loop
        for chain in list_couple_chain:
            chain1,chain2 = chain
            
            # index of fisrt chain
            for index1 in pdb.index[chain1]:
                
                
                coord1 = pdb.index2coord[chain1][index1]
                aa1    = pdb.index2sequence[chain1][index1]
                
                
                # aa1 is a positif charged 
                if aa1 in positively:
                    complementar = aromatic
                    coord_aa1 = [pdb.index2coord[chain1][index1][i] for i in charged_atoms[aa1] if i in pdb.index2coord[chain1][index1].keys()] 
                    
                # aa1 is a aromatic
                elif aa1 in aromatic:
                    complementar = positively
                    coord_aa1   = [pmath.get_barycenter_from_list([pdb.index2coord[chain1][index1][i] for i in last_atom[aa1].split("_") if i in pdb.index2coord[chain1][index1].keys()])]
                    
                    
                # none we skip the loop
                else:
                    continue
                
                
                
                
                # index of second chain
                for index2 in pdb.index[chain2]:
                    
                    
                    # keep after the redundancies
                    if chain1 == chain2 and int(index1) >= int(index2):
                        continue
                    
                    
                    coord2 = pdb.index2coord[chain2][index2]
                    aa2    = pdb.index2sequence[chain2][index2]
                    
                    # aa2 is complementar for cation pi interaction
                    if aa2 not in complementar:
                        continue
    
                
                    if aa2 in positively:
                        coord_aa2 = [pdb.index2coord[chain2][index2][i] for i in charged_atoms[aa2] if i in pdb.index2coord[chain2][index2].keys()] 
                        
                    # aa1 is a aromatic
                    elif aa2 in aromatic:
                        coord_aa2 = [pmath.get_barycenter_from_list([pdb.index2coord[chain2][index2][i] for i in last_atom[aa2].split("_") if i in pdb.index2coord[chain2][index2].keys()])]
                    
                    
                    # Distance analyse
                    list_distance = []
                    for c1 in coord_aa1:
                        for c2 in coord_aa2:
                            d = self.getDistanceContact(c1,c2)
                            if d < self.dicvar["threshold"]:
                                list_distance.append(d)
                    
                    # If there is at least one distance in list         
                    if len(list_distance)!=0:
                        printer += "%s%s%s_%s %s%s%s_%s %.2f\n"%(aa1,index1,chain1,"X",aa2,index2,chain2,"X",math.sqrt(min(list_distance)))

        # For lib use
        if adress != "X":
            f = open(adress,"w")                        
            f.write(printer)
            f.close()
            
        return printer


    
    def infoResContact(self,coord1,coord2,type="AA"):
        """ """
        
        #backbone_atoms = ["N","O","C","CA"]
        
        try:
            if self.getDistanceContact(coord1["CA"], coord2["CA"]) > self.dicvar["threshold_limit"]:
                return False,None
        except:
            pass
        
        
        list_atom_contact = []
        for atom1 in coord1:
            
            for atom2 in coord2:
                
                d = self.getDistanceContact(coord1[atom1],coord2[atom2])
                
                if d < self.dicvar["threshold"]:
                    list_atom_contact.append([atom1,atom2])
        
        return True,list_atom_contact
        
    
        
    def computeContact(self,isContact):
        """ Search the contacts """
        
        c1,c2 = self.fillChainsContact()
        chain1 = self.data_pdb[c1]
        chain2 = self.data_pdb[c2]
        
        for index1 in chain1["order"]:
            for index2 in chain2["order"]:
                
                try: # if there is an issue like atoms missing
                    if isContact(chain1[index1],chain2[index2]):
                        
                        # IF Contact launch from the shell
                        if self.dicvar["dump"] and not self.dicvar["interface"]:
                            print("%s_%s %s_%s"%(c1,index1,c2,index2))
                        
                        
                        # IF not from the shell OR interface asked -> access to the contact in self.contact_list["contact"]
                        if self.dicvar["interface"] or not self.dicvar["dump"]:
                            self.contact_list["contact"].append("%s_%s"%(index1,index2))
                        
                        # IF interface asked
                        if self.dicvar["interface"]:
                            self.contact_list["interA"].append(index1)
                            self.contact_list["interB"].append(index2)
                except:
                    pass
                
        return

    
    def computeCAPRIContact(self,isContact):
        """ Search the contacts """

        c1,c2 = self.fillChainsContact()
        chain1 = self.data_pdb[c1]
        chain2 = self.data_pdb[c2]
        


        for index1 in chain1["order"]:
            for index2 in chain2["order"]:
                try: # if there is an issue like atoms missing
                    analyse = isContact(chain1[index1],chain2[index2])
                    if analyse[0]:
                        self.contact_list["interA"].append(index1)
                        self.contact_list["interB"].append(index2)
                        if analyse[1]:
                            self.contact_list["contact"].append("%s_%s"%(index1,index2))
                except:
                    print("exception in computeCAPRIContact")
                    pass
        
        self.contact_list["interA"] = list(set(self.contact_list["interA"]))
        self.contact_list["interB"] = list(set(self.contact_list["interB"]))
        
        
        return
         
    
    def fillChainsContact(self):
        """ """
        
        # If chains are not stipuled we take the first and the second in the PDB file
        if not "chains" in self.dicvar.keys():
            self.dicvar["chains"] = (self.data_pdb["order"][0],self.data_pdb["order"][1])
            return self.data_pdb["order"]
        
        # Else we extract
        #
        # TODO : multi chain management and same chain 
        #
        return self.dicvar["chains"].split(":")
            
        
    def forCAPRI(self,res1,res2):
        """ Compute if there is a contact AA"""
        
        CA = self.getDistanceContact(res1["coord"]["CA"],res2["coord"]["CA"])
        
        # Optimisation filter : we compute precisely the contact if the CA are below 20A
        if CA < 400.0:
            dinter = False
            for a in res1["coord"]:
                for b in res2["coord"]:
                    d = self.getDistanceContact(res1["coord"][a], res2["coord"][b])
                    if not dinter:
                        if d < 100.0:
                            dinter = True
                    if d <25.0:
                        return True,True
            if dinter:
                return True,False
        return False,False
        
    
    def isContactAA(self,res1,res2):
        """ Compute if there is a contact AA"""
        
        CA = self.getDistanceContact(res1["coord"]["CA"],res2["coord"]["CA"])
        
        # Optimisation filter : we compute precisely the contact if the CA are below 20A
        if CA < self.dicvar["threshold_limit"]:
            for a in res1["coord"]:
                for b in res2["coord"]:
                    if self.getDistanceContact(res1["coord"][a], res2["coord"][b]) < self.dicvar["threshold"]:
                        return True
        return False
    
    def isContactCA(self,res1,res2):
        """ Compute if there is a contact CA"""
        
        try: # Try if residue have CA atoms
            if self.getDistanceContact(res1["coord"]["CA"],res2["coord"]["CA"])< self.dicvar["threshold"]:
                return True
            else:
                return False
        except: # residu have no CA
            return False
    
    
    def isContactBB(self,res1,res2):
        """ Compute if there is a contact BB"""
        
        CA = self.getDistanceContact(res1["coord"]["CA"],res2["coord"]["CA"])
        
        # Optimisation filter : we compute precisely the contact if the CA are below 20A
        if CA < self.dicvar["threshold_limit"]:
            atoms_BB = ["CA","C","O","N"]
            
            for a in atoms_BB:
                for b in atoms_BB:
                    if self.getDistanceContact(res1["coord"][a],res2["coord"][b]) < self.dicvar["threshold"]:
                        return True
        return False

    
    def isContactSC(self,res1,res2):
        """ Compute if there is a contact SC"""
        
        CA = self.getDistanceContact(res1["coord"]["CA"],res2["coord"]["CA"])
        
        # Optimisation filter : we compute precisely the contact if the CA are below 20A
        
        if CA < self.dicvar["threshold_limit"]:
            
            # If there is a glycine we consider SC as a CA of the glycine
            gly1 = False
            gly2 = False
            if res1["aa"] == "GLY":
                atoms_SC_res1 = ["CA"]
                gly1 = True
            if res2["aa"] == "GLY":
                atoms_SC_res2 = ["CA"]
                gly2 = True
            
            # If no glycine we selecte the side chain residus
            if not gly1:
                atoms_SC_res1 = [at for at in res1["coord"] if at not in ["CA","C","O","N"]]
            if not gly2:
                atoms_SC_res2 = [at for at in res2["coord"] if at not in ["CA","C","O","N"]]
        
            for a in atoms_SC_res1:
                for b in atoms_SC_res2:
                    if self.getDistanceContact(res1["coord"][a],res2["coord"][b]) < self.dicvar["threshold"]:
                        return True
        return False
    
    
    
    def getDistanceContact(self,c1,c2):
        """ Get the square distance """
        b = list(map(self.f,c1,c2))
        return b[0]*b[0]+b[1]*b[1]+b[2]*b[2]
        

    def computeInterfaceContact(self):
        """ Format the interface residue """
        
        self.contact_list["interA"],self.contact_list["interB"] = map(lambda x:list(set(x)),[self.contact_list["interA"],self.contact_list["interB"]])
        
        return

    def screenInterfaceContact(self):
        """ Screen if launch from the shell """
        print("Chain1 ")
        for i in self.contact_list["interA"]:
            print("%s"%(i))
            
        print("Chain2 ")
        for i in self.contact_list["interB"]:
            print("%s"%(i))
        
        return
        
        
        
        
        
""" MAIN """

if __name__=="__main__":

    dicvar = {}
    
    dicvar["dump"] = True

    
    try:
        dicvar["pdb"] = sys.argv[sys.argv.index("-pdb")+1]
        dicvar["threshold"] = float(sys.argv[sys.argv.index("-t")+1])
    except:
        usage()
        sys.exit(1)
    
    try:
        dicvar["type"] = sys.argv[sys.argv.index("-co")+1]
    except:
        pass
    
    try:
        dicvar["chains"] = sys.argv[sys.argv.index("-c")+1]
    except:
        pass
    
    try:
        test = sys.argv[sys.argv.index("-inter")]
        dicvar["interface"] = True
    except:
        pass
    
    try:
        dicvar["all"] = sys.argv[sys.argv.index("-all")+1]
    except:
        pass
    
    try:
        dicvar["hc"] = sys.argv[sys.argv.index("-hc")+1]
    except:
        pass
    
    try:
        dicvar["cc"] = sys.argv[sys.argv.index("-cc")+1]
    except:
        pass
   
    try:
        dicvar["cp"] = sys.argv[sys.argv.index("-cp")+1]
    except:
        pass
   
   
   
    Contact(dicvar)
