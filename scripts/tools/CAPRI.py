import os,sys,glob,re
import tempfile
import configparser
#import subprocess

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(os.path.dirname(script_path))

config = configparser.ConfigParser()
config.read(f'{script_dir}/config.ini')

PROFIT = config['PROGRAMS']['PROFIT']

# import time
sys.path.append(os.path.join(script_dir, "tools"))
import ContactForCAPRI as Contact

try:
    import quick_fnat
    QUICK_FNAT_OK = True
except ImportError:
    QUICK_FNAT_OK = False
    print("Import Error: quick_fnat not found.")

def usage():
    print("""CAPRI.py -ref <ref_pdb> -mob <pdb_mobile> -chain <AB:AB>
                            [-mob can be a regex -> WARNING before each special character use \ i.e : complex\*
                      or    [-lm        : list of mobile structure
                            
                            if no arguments -> All steps are computed
                            [-lrms  : Compute the lrms and rms  WARNING : this step do not align sequence so, you want to align add -fr  
                            [-irms  : Compute the irms
                            [-fr    : Compute fractions Native Non Native and Interface Residues
                            
                            [-f    : write the results on the dump file
                                
                            [-quick_fnat : to use prefiltering with fnat (all decoys < 10% fnat are incorrect)]
                            [-quick_fnat_ref : optional alternative reference PDB to use when prefiltering with fnat (all decoys < 10% fnat are incorrect)]

                python CAPRI.py -ref 1GXD_complex.pdb -mob Complex_30g.pdb -f results.rmsd -chain AB:AB
                python CAPRI.py -ref 1GXD_complex.pdb -mob Complex_\*.pdb -f results.rmsd -chain AB:AB
                python CAPRI.py -ref 1GXD_complex.pdb -lm list_mobile.list -f results.rmsd -chain AB:AB
                
                
                STDOUT   FILE RMSAVSA RMSBSVB RMSABVSAB LRMSA LRMSB FRNAT FRNNAT FRIR IRMS CAPRIRANK DOCKQ DOCKQRANK
                
            DOCKQ implemented as in Basu & Wallner, PLOSONE, 2016  (https://doi.org/10.1371/journal.pone.0161879) with PROFIT alignments

            PS: -chain AB:XY means chain A in reference is considered as Receptor. Since CAPRI criteria see only L-RMS but not 
                R-RMS, results of -chain AB:XY differ from BA:YX
            @author: Guerois'Team LBSR CEA

            RMK: /!\ be careful, atom names have to start in column 14 (i.e. with a space) if name has less than 4 letters/numbers 
                     otherwise ProFit won't be able to select CA, C, O and N atoms and you'll end up with an error.
            """)
    sys.exit(1)


class CAPRI:
    """
    CAPRI(dicvar)  
                            have key  "ref"      : adress of reference pdb
                            have key  "chain"    : for instance AB:AB
                            have key  "mob"      : list of mobile or just adress of one mobile, or regex
                        or  have key   <filelist> : file containing lists of mobile structure
                            may key   <lrms>     ; compute lrms and rms
                            may key   <irms>     ; compute irms
                            may key   <fraction> : compute the fraction
                            may key   <f>        : adress of file for the results
                            
                
            @author: Guerois'Team LBSR CEA
    
    """

    def __init__(self,dicvar):

        self.dicvar = dicvar
        
        self.fillOptionCAPRI()

        self.runCAPRI()

        
        return
    
    
    def fillOptionCAPRI(self):
        """ Fill Options by default """
        
        if not "mobile" in self.dicvar.keys():
            self.dicvar["mobile"] = ""
        
        # REGEX list
        if "*" in self.dicvar["mobile"]:
            self.dicvar["mobile"] = glob.glob(self.dicvar["mobile"])
            
        # File
        elif "liste_mob" in self.dicvar.keys() :
            f = open(self.dicvar["liste_mob"],"r")
            l = f.readlines()
            self.dicvar["mobile"] = [i.strip() for i in l]
            f.close()
        
        # Unique 
        elif os.path.isfile(self.dicvar["mobile"]):
            self.dicvar["mobile"] = [self.dicvar["mobile"]]
            
        else:
            print("File not found")
            raise IOError  # Raise the error and stop the program
        
        if not "quick_fnat" in self.dicvar:
            self.dicvar["quick_fnat"] = False
        
        else:
            if not "quick_fnat_ref" in self.dicvar:
                self.dicvar["quick_fnat_ref"] = self.dicvar["reference"]
        
        # ELSE -> pdb file or list of pdb file
        
        #self.chain1 = "A"
        #self.chain2 = "B"
        self.chain1 = self.dicvar["chain"][0]
        self.chain2 = self.dicvar["chain"][1]
        self.mobchain1 = self.dicvar["chain"][3]            
        self.mobchain2 = self.dicvar["chain"][4]
        
        
        self.root_for_temp = os.path.basename(self.dicvar["reference"]).split(".")[0]
        
        
        # Options
        if not "lrms" in self.dicvar.keys() and not "fraction" in self.dicvar.keys() and not "irms" in self.dicvar.keys():
            self.dicvar["lrms"]     = True
            self.dicvar["fraction"] = True
            self.dicvar["irms"]     = True
        
        if not "lrms" in self.dicvar.keys():
            self.dicvar["lrms"] = False
        
        if not "fraction" in self.dicvar.keys():
            self.dicvar["fraction"] = False
        
        if not "irms" in self.dicvar.keys():
            self.dicvar["irms"] = False
        
        if not "dump" in self.dicvar.keys():
            self.dicvar["dump"] = False
        
        if not "file" in self.dicvar.keys():
            self.dicvar["file"] = None
        else:
            self.file = open(self.dicvar["file"],"w")
            self.file.write("#FILE RMSAVSA RMSBSVB RMSABVSAB LRMSA LRMSB FRNAT FRNNAT FRIR IRMS CAPRIRANK DOCKQ DOCKQRANK\n")
        
        if not "verbose" in self.dicvar.keys():
            self.dicvar["verbose"] = False
            
        return


    def filter_with_quick_fnat(self):

        if not QUICK_FNAT_OK:
            print("Warning: couldn't load quick_fnat module, proceeding without this option")
            return

        if self.dicvar["verbose"]:
            print("Analyzing %d decoys with quick_fnat"%len(self.dicvar["mobile"]))

        lst_decoys_fnat10p = []
        chainsR, chainsL = self.dicvar["chain"].split(":")
        chainR1, chainR2 = chainsR
        chainL1, chainL2 = chainsL
        contact_ref = quick_fnat.get_contact(self.dicvar["quick_fnat_ref"], chainR1, chainR2)
        #contact_ref = quick_fnat.get_contact(self.dicvar["reference"], chainR1, chainR2)
        print(contact_ref)
        for mobile in self.dicvar["mobile"]:
            # fnat 
            contact_pred = quick_fnat.get_contact(mobile, chainL1, chainL2)
            print(contact_pred)
            fnat = quick_fnat.compute_fnat(contact_ref, contact_pred)
            if fnat >= 10:
                lst_decoys_fnat10p.append(mobile)

        if self.dicvar["verbose"]:
            print("%d structures have fnat >= 10"%len(lst_decoys_fnat10p))

        self.dicvar["mobile"] = lst_decoys_fnat10p
        
        return

    
    def runCAPRI(self):
        """ Runing Step """
        
        # Extract the Chains
        self.getChainCAPRI()
        
        # Contact computing
        #self.ref_contact,self.ref_index, self.ref_data_pdb = self.getContactCAPRI(self.dicvar["reference"],self.chain1,self.chain2)
        
        self.ref_contact,self.ref_interface, self.ref_index, self.ref_data_pdb = self.getInterfaceAndContactCAPRI(self.dicvar["reference"],self.chain1,self.chain2)
        
        if self.dicvar["irms"]:
            # interace computing
            #self.ref_interface     = self.getInterfaceCAPRI(self.dicvar["reference"],self.chain1,self.chain2)
            self.ref_pdb_interface = self.writePdbInterfaceCAPRI(self.dicvar["reference"],self.ref_interface,self.ref_data_pdb)
            # A file which only contains the coordinates of the interface backbone atoms
            
        
        if "quick_fnat" in self.dicvar and self.dicvar["quick_fnat"]:
            self.filter_with_quick_fnat()


        for mobile in self.dicvar["mobile"]:
            
            # Contact and index of mobile structure
            #mob_contact,mob_index,mob_data_pdb = self.getContactCAPRI(mobile,self.chain1,self.chain2)
            mob_contact, mob_interface, mob_index, mob_data_pdb = self.getInterfaceAndContactCAPRI(mobile,self.mobchain1,self.mobchain2)
            #mob_contact, mob_interface, mob_index, mob_data_pdb = self.getInterfaceAndContactCAPRI(mobile,self.chain1,self.chain2)
            
            # Same index ? """ Test if the index number are the same to know if we need alignment or not """
            align = not self.sameIndexCAPRI(mob_index) 
            ### STEP 1
            # Compute the RMS and lRMS
            
            self.rmsAvsA,self.rmsBvsB,self.rmsABvsAB,self.LrmsA,self.LrmsB = self.getRmsAndLrmsCAPRI(mobile,align)

            # Get corresponding if structure are not the same
            if align:
                self.managementIndexAlignCAPRI(mob_index)
            
            
            
            if self.dicvar["fraction"]:
                # Compute the fractions : Native, Non Native, Interface Residue
                #mob_interface = self.getInterfaceCAPRI(mobile,self.chain1,self.chain2)
                self.fractionNat,self.fractionNonNat,self.fractionIR = self.fractionCAPRI(mob_contact,mob_interface,align)
                
            if self.dicvar["irms"]:
                # Compute the iRMS
                # Calculate the ref interface depending on the align parameter
                self.ref_pdb_interface = self.writePdbInterfaceCAPRI(self.dicvar["reference"],self.ref_interface,self.ref_data_pdb,refboolean=align)
                mob_pdb_interface = self.writePdbInterfaceCAPRI(mobile,self.ref_interface,mob_data_pdb,align,mobboolean=True)
                #mob_pdb_interface = self.writePdbInterfaceCAPRI(mobile,self.ref_interface,mob_data_pdb,align)
                self.Irms = self.getIrmsCAPRI(mob_pdb_interface)
                
            if self.dicvar["lrms"] and self.dicvar["irms"] and self.dicvar["fraction"]:
                # Capri ranking
                self.rank = self.getRankingCAPRI()
                # DockQ score (self.DockQ) and ranking
                self.rankQ = self.getDockQ()
            
            # Print on the file and on the screen the result
            self.screenerCAPRI(mobile,align)
        
        if self.dicvar["file"]:
            self.file.close()
        
        # We delete the pdb interface of the reference
        if self.dicvar["irms"]:
            os.system("rm %s"%(self.ref_pdb_interface))
        
        return
    
    def getChainCAPRI(self):
        """ """
        
        chain =  ["%s*:%s*"%(self.chain1,self.mobchain1),"%s*:%s*"%(self.chain2,self.mobchain2)]
        self.chain = chain
        return chain
    
            
    def getInterfaceAndContactCAPRI(self,pdb,chain1,chain2):
        """ Get contact between 2 chains, index and all data of pdb"""
        
        # Search contact 5A
        contact_object = Contact.Contact({"pdb":pdb,\
                                          "chains":chain1+":"+chain2,\
                                          "capri":True})
        # Stockage
        contact= contact_object.contact_list["contact"]
        
        # Index Stockage
        index = [contact_object.data_pdb[chain1]["order"],contact_object.data_pdb[chain2]["order"]]

        # Stockage
        inter = {chain1: contact_object.contact_list["interA"],chain2: contact_object.contact_list["interB"]}
        
        # All stokage
        data_pdb = contact_object.data_pdb
        
#         print contact
#         print inter
#         print index
#         print len(data_pdb)
        return contact,inter,index,data_pdb

    
    
    def writePdbInterfaceCAPRI(self,pdb,interface,data_pdb,align=False,refboolean=False,mobboolean=False):
        """ """
        def writeInterfacePdb(f_pdb_inter,inter,chain,mobchain,data_pdb):
            nb_atom = 1

            if align:
#                # New inter for mobile
#                interface_new  = [self.ref2mob[chain][i] for i in inter[chain]]
#                # Ref index
#                interface_ref  = inter[chain]
#                interface_cor = zip(interface_new,interface_ref)
#
                # New inter for mobile
                interface_new = []
                interface_newref = []
                for i in inter[chain]:
                    # We need to check that the position exists in both sequences
                    if i in self.ref2mob[chain].keys() and self.ref2mob[chain][i] in self.mob2ref[mobchain].keys():
                        interface_new.append(self.ref2mob[chain][i])
                        interface_newref.append(i)
                    else:
                        if i in self.ref2mob[chain].keys():
                            print(self.ref2mob[chain][i])
                # Ref index
                interface_ref  = interface_newref
                interface_cor = zip(interface_new,interface_ref)          
            else:
                if refboolean:
                    # This means we are calculating the reference for a pdb with a sequence different from the other pdb
                    interface_new = []
                    interface_newref = []
                    for i in inter[chain]:
                        if i in self.ref2mob[chain].keys() and self.ref2mob[chain][i] in self.mob2ref[mobchain].keys():
                            interface_new.append(i)
                            interface_newref.append(i)
                    interface_ref =  interface_newref
                    interface_cor = zip(interface_new,interface_ref)
                else:
                    # We keep the ref index
                    interface_new = inter[chain]
                    interface_ref = inter[chain]
                    interface_cor = zip(interface_new,interface_ref)
                
            for mob,ref in interface_cor:
                
                for at in  ["N","CA","C","O"]:
                    
                    
                    try:
                    
                        f_pdb_inter.write("ATOM%7d  %-3s %3s %s%4s    %8.3f%8.3f%8.3f  1.00  1.00\n"%(nb_atom,\
                                                                                   at,\
                                                                                   data_pdb[mob]["aa"],\
                                                                                   chain,\
                                                                                   ref,\
                                                                                   data_pdb[mob]["coord"][at][0],\
                                                                                   data_pdb[mob]["coord"][at][1],\
                                                                                   data_pdb[mob]["coord"][at][2]))
                        nb_atom += 1
                    except:
                        print("Couldn't write",chain,mob,ref,at)
                        pass
        
            return
        ###############################################################

        ad_pdb_inter = tempfile.mktemp(suffix = "IRMS_CAPRI",prefix= self.root_for_temp)
        f_pdb_inter  = open(ad_pdb_inter,"w")
        
        # write interface residue chainA
        if mobboolean:
            writeInterfacePdb(f_pdb_inter,interface,self.chain1,self.mobchain1,data_pdb[self.mobchain1])
        else:
            writeInterfacePdb(f_pdb_inter,interface,self.chain1,self.mobchain1,data_pdb[self.chain1])
        f_pdb_inter.write("TER\n\n")
        
        # write interface residue chainB
        if mobboolean:
            writeInterfacePdb(f_pdb_inter,interface,self.chain2,self.mobchain2,data_pdb[self.mobchain2])
        else:
            writeInterfacePdb(f_pdb_inter,interface,self.chain2,self.mobchain2,data_pdb[self.chain2])
        
        # Closing and removing the file
        f_pdb_inter.close()
        
        return ad_pdb_inter
        



    
    def getInterfaceContactCAPRI(self,pdb,chain1,chain2):
        """ Get information about interface residue and contact of reference"""
        
        # Search interface residue 10A
        contact_object = Contact.Contact({"pdb":pdb,\
                                          "threshold":10,\
                                          "interface":True})
        
        # Stockage
        inter = {chain1: contact_object.contact_list["interA"],chain2: contact_object.contact_list["interB"]}
        
        
        # Search contact 5A
        contact_object = Contact.Contact({"pdb":pdb,\
                                          "threshold":5\
                                      })
        # Stockage
        contact= contact_object.contact_list["contact"]
        
        # Index Stockage
        index = [contact_object.data_pdb[chain1]["order"],contact_object.data_pdb[chain2]["order"]]
        
        
        return inter,contact,index
    
    def sameIndexCAPRI(self,index_mobile):
        """ Test if the index number are the same to know if we need alignment or not """
        
        if index_mobile == self.ref_index:
            return True
        else:
            return False
        
        return
    
    
    def fractionCAPRI(self,mobile_contact,mob_interface,align=False):
        """ Compute any Fraction : Native, Non Native and Interface Residue"""
        def tradContact(contact,align=False):
            if align:
                res1,res2 = contact.split("_")
                # Old Contact.cpp
                #sp = contact.split("_")
                #chain1 = sp[0][0]
                #chain2 = sp[1][0]
                #res1   = sp[0][1:]
                #res2   = sp[1][1:]
                if res1 in self.ref2mob[self.chain1].keys() and res2 in self.ref2mob[self.chain2].keys(): 
                    return "%s_%s"%(self.ref2mob[self.chain1][res1],self.ref2mob[self.chain2][res2])
                    #return [self.ref2mob[self.chain1][res1],self.ref2mob[self.chain2][res2]]
                    #return "%s%s_%s%s"%(chain1,self.ref2mob[chain1][res1],chain2,self.ref2mob[chain2][res2])
                return "NONE"
            else:
                return contact
            return
        
        
        def tradRes(res,chain,align=False):
            if align:
                if res in self.ref2mob[chain].keys():
                    return self.ref2mob[chain][res]
                else:
                    return "NONE"
            else:
                return res
            return
        
        # For Native and Non Native fraction
        contact_ref = self.ref_contact
        contact_mob = mobile_contact
        
        """
        if not a :
            nbSameNat = 0
            for cur_contact in contact_ref:
                if tradContact(cur_contact,align) in contact_mob:
                    nbSameNat += 1
        else: 
        """
        # Optimised version
        nbSameNat = sum([1 for cur_contact in contact_ref if tradContact(cur_contact,align) in contact_mob])
        
        #nbSameNat = len(set(contact_mob).intersection(contact_ref))
        
        
        # We counted 2 times A->B and B->A before so for the amount of contact we have to divide by 2
        nbNat  = len(contact_ref)      #nbNat  = len(contact_pdb1)/2
        nbPred = len(contact_mob)      #nbPred = len(contact_pdb2)/2
            
        
        # Computing Native Fraction & Non Native Fraction
        try:
            fractionNat    = (float(nbSameNat)*100.0) / float(nbNat)
        except:
            fractionNat = 0.0
            
        try:
            fractionNonNat = (float(nbPred-nbSameNat)*100.0) / float(nbPred)
        except:
            fractionNonNat = 0.0
        
        ### For IR Fraction
        # Get Interface residue   "A"->[1,2,3...]  "B"->[1,2,3,...]
        interface_ref = self.ref_interface
        interface_mob = mob_interface
        
        # Nb Interface residue in reference and in Mobile structutre
        NbResInterRef = len(interface_ref[self.chain1])+len(interface_ref[self.chain2])
#         NbResInterMob = len(interface_mob[self.mobchain1])+len(interface_mob[self.mobchain2])
        
        
        # Nb Common Interface Residue  
        if align:
            NbCommonInter = 0
            for cur_int_res in interface_ref[self.chain1]:
                if tradRes(cur_int_res,self.chain1,align) in interface_mob[self.mobchain1]:
                    NbCommonInter += 1
            for cur_int_res in interface_ref[self.chain2]:
                if tradRes(cur_int_res,self.chain2,align) in interface_mob[self.mobchain2]:
                    NbCommonInter += 1
        else:
            # Optimised version !
            NbCommonInter = len(set(interface_ref[self.chain1]).intersection(interface_mob[self.mobchain1])) + len(set(interface_ref[self.chain2]).intersection(interface_mob[self.mobchain2]))
        
        # Nb Common Interface Residue  : Usage set class  set(list1).intersection(list2) give the common between list1 and list2
        #NbCommonInter = len(set(interface_ref[self.chain1]).intersection(interface_mob[self.chain1])) + len(set(interface_ref[self.chain2]).intersection(interface_mob[self.chain2]))
        
        try:
            fractionIR = (float(NbCommonInter)*100.0) / float(NbResInterRef)
        except:
            fractionIR = 0.0
        
        return fractionNat,fractionNonNat,fractionIR
        
        
        
    def getRmsAndLrmsCAPRI(self,mobile,align=False):
        """ Compute rms et Lrms """
        
        # TEMPORY FILE FOR FILLING THE PROFIT OPTIONS
        protocole_file = tempfile.mktemp(prefix=self.root_for_temp,suffix="Capri_ProFit")
        
        # TEMPORY FILE FOR THE NW ALIGNMENT GENERATE BY PROFIT
        self.alignment_file = tempfile.mktemp(prefix=self.root_for_temp,suffix="Capri_Alig")
        
        
        ### PROFIT OPTION FILE
        f_profit = open(protocole_file,"w")
        
        # OPTIONS PREPARATION FOR PROFIT
        protocole  = "REFERENCE %s\n"%(self.dicvar["reference"])        # ADRESS OF REFERENCE
        protocole += "MOBILE %s\n"%(mobile)           # ADRESS OF MOBILE STRUCTURE
        protocole += "ATOMS CA,C,O,N\n"                            # ATOMS SELECTED FOR FITTING
        protocole += "IGNOREMISSING\n"                  # because sometimes there is a missing 0 in the last res of Zdock decoy
        
        
        if align:
            protocole += "ALIGN %s\n"%(self.chain[0])                  # SELECTED RESIDUES OF THE FIRST CHAIN
            protocole += "IGNOREMISSING\n"                  # because sometimes there is a missing 0 in the last res of Zdock decoy
            protocole += "FIT\n"                                       # STRUCTURAL ALIGNMENT BETWEEN FIRST CHAIN
            protocole += "ALIGN %s\n"%(self.chain[1])                  # SELECTED RESIDUES OF THE SECOND CHAIN
            protocole += "NOFIT\n"                                     # NO FITTING ! WE KEEP THE LAST FITTING IN MEMORY
            protocole += "RMS\n"                                       # COMPUTING RMS ON THE CHAIN B WITH ALIGNEMENT ON CHAIN A --> LRMS  A is receptor
            protocole += "IGNOREMISSING\n"                  # because sometimes there is a missing 0 in the last res of Zdock decoy
            protocole += "FIT\n"                                       # STRUCTURAL ALIGNMENT BETWEEN SECOND CHAIN
            protocole += "ALIGN %s\n"%(self.chain[0])                  # SELECTED RESIDUES OF THE SECOND CHAIN
            protocole += "NOFIT\n"                                     # NO FITTING ! WE KEEP THE LAST FITTING IN MEMORY
            protocole += "RMS\n"                                       # COMPUTING RMS ON THE CHAIN B WITH ALIGNEMENT ON CHAIN A --> LRMS B is receptor
            protocole += "ALIGN %s\n"%(self.chain[0])       # corrected by Jinchao
            protocole += "ALIGN %s APPEND\n"%(self.chain[1])       # corrected by Jinchao
            protocole += "IGNOREMISSING\n"                  # because sometimes there is a missing 0 in the last res of Zdock decoy
            protocole += "FIT\n"                                       # STRUCTURAL ALIGNMENT BETWEEN ALL
            protocole += "PRINTALIGN FASTA %s\n"%(self.alignment_file) # PUT THE NW ALIGNMENT IN A FILE
      
        if not align:
            protocole += "ZONE %s\n"%(self.chain[0])                   # SELECTED RESIDUES OF THE FIRST CHAIN
            protocole += "IGNOREMISSING\n"                  # because sometimes there is a missing 0 in the last res of Zdock decoy
            protocole += "FIT\n"                                       # STRUCTURAL ALIGNMENT BETWEEN FIRST CHAIN
            protocole += "RZONE %s\n"%(self.chain[1])                  # LRMSA
            protocole += "ZONE CLEAR \n"                               # CLEAR ALL ZONE
            protocole += "ZONE %s\n"%(self.chain[1])                   # SELECTED RESIDUES OF THE SECOND CHAIN
            protocole += "IGNOREMISSING\n"                  # because sometimes there is a missing 0 in the last res of Zdock decoy
            protocole += "FIT\n"                                       # STRUCTURAL ALIGNMENT BETWEEN SECOND CHAIN
            protocole += "RZONE %s\n"%(self.chain[0])                  # LRMSB
            protocole += "ZONE CLEAR \n"                               # CLEAR
            protocole += "ZONE *\n"                                    # SELECTED ALL
            protocole += "IGNOREMISSING\n"                  # because sometimes there is a missing 0 in the last res of Zdock decoy
            protocole += "FIT\n"                                       # RMSDABVSAB
            
        
        
        f_profit.write(protocole)
        f_profit.close()
        ###
        
        ### PROFIT LAUNCH !
        out_profit = tempfile.mktemp(suffix="PROFIT",prefix=self.root_for_temp)
        os.system(PROFIT+" -f %s > %s"%(protocole_file,out_profit))
        f_out       = open(out_profit,"r")
        line_profit = f_out.read()
        f_out.close()
        #line_profit = os.popen(PROFIT+" -f %s"%(protocole_file)).read()
        
        ### RMS PARSING 
        # rms [AvsA,LRMS_recA,BvsB,LRMS_recB,ABvsAB]
        rms = [float(i) for i in re.compile("RMS: (.*?)\n").findall(line_profit)]
        
        if len(rms) < 5:
            print("Error: Problem with ProFit, couldn't extract rms from output files")
            print(protocole_file, out_profit)
            print(line_profit)
            sys.exit()

        # CLEAN THE TEMPORARY FILES
        os.system("rm %s"%(protocole_file))
        os.system("rm %s"%(out_profit))

        return rms[0],rms[2],rms[4],rms[1],rms[3]   # AvsA,BvsB,ABvsAB,LRMSrecA,LRMSrecB



    def getIrmsCAPRI(self,mobile):
        """ Compute the Irms """
        
        # TEMPORY FILE FOR FILLING THE PROFIT OPTIONS
        protocole_file = tempfile.mktemp(prefix=self.root_for_temp,suffix="Capri_ProFit")
        
        ### PROFIT OPTION FILE
        f_profit = open(protocole_file,"w")
        
        # OPTIONS PREPARATION FOR PROFIT
        protocole  = "REFERENCE %s\n"%(self.ref_pdb_interface)        # ADRESS OF REFERENCE
        protocole += "MOBILE %s\n"%(mobile)           # ADRESS OF MOBILE STRUCTURE
        protocole += "ZONE %s*:%s*\n"%(self.chain1,self.chain1)           # SELECT ZONES
        protocole += "ZONE %s*:%s*\n"%(self.chain2,self.chain2)           # SELECT ZONES
        protocole += "ATOMS CA,C,O,N\n"                            # ATOMS SELECTED FOR FITTING
        protocole += "IGNOREMISSING\n"                  # because sometimes there is a missing 0 in the last res of Zdock decoy
        protocole += "FIT\n"                                       # STRUCTURAL ALIGNMENT BETWEEN FIRST CHAIN
        
        f_profit.write(protocole)
        f_profit.close()
        
        ### PROFIT LAUNCH !
        out_profit = tempfile.mktemp(suffix="PROFIT",prefix=self.root_for_temp)
        os.system(PROFIT+" -f %s > %s"%(protocole_file,out_profit))
        f_out       = open(out_profit,"r")
        line_profit = f_out.read()
        f_out.close()
        
        # PDB Interface deleting
        os.system("rm %s"%(mobile))
        os.system("rm %s"%(out_profit))
        os.system("rm %s"%(protocole_file))
        
        try:
            return float(re.compile("RMS: (.*?)\n").findall(line_profit)[0])
        except:
            print("Error: No RMS found in PROFIT output")
            print(line_profit)
            print("IRMS set as 999.99 (error value)")
            return 999.99

    
    
    def managementIndexAlignCAPRI(self,index):
        """ Create correspondance between struct1 and struct2"""

        def getcorrespondance(indexref,indexmob,seqref,seqmob):
            
            iteseq = 0
            itemob = 0
            iteref = 0
            
            pdb2alignref     = {}  # index residue in the pdb file -> position in the alignement
            pdb2alignmob     = {}  # same for mobile structure
            ref2mob          = {}  # index residue in the reference pdb file -> index residue in the mobile pdb file
            mob2ref          = {}  # inverse of the previous dictionnary
            
            
            for sref,smob in zip(seqref,seqmob):
                if sref != "-":
                    pdb2alignref[indexref[iteref]] = iteseq
                    if smob != "-":
                        ref2mob[indexref[iteref]]  = indexmob[itemob]
                
                if smob != "-":
                    pdb2alignmob[indexmob[itemob]] = iteseq
                    if sref != "-":
                        mob2ref[indexmob[itemob]]  =  indexref[iteref]
                
                
                # iteration process
                iteseq += 1
                if sref != "-": iteref += 1
                if smob != "-": itemob += 1
                
            return pdb2alignref,pdb2alignmob,ref2mob,mob2ref
            

        
        # Alignment ["pdb1"/"pdb2"/"order"]["A"/"B"] -> sequences
        align_seq = self.alignmentTraitementCAPRI()  
        seq_ref    = align_seq[align_seq["order"][0]]
        seq_mob    = align_seq[align_seq["order"][1]]
        
        # index [ChainA,ChainB]
        index_ref = self.ref_index
        index_mob = index
        
        pdb2alignref = {}
        pdb2alignmob = {}
        ref2mob      = {}
        mob2ref      = {}
        
        pdb2alignref[self.chain1],pdb2alignmob[self.mobchain1],ref2mob[self.chain1],mob2ref[self.mobchain1] = getcorrespondance(index_ref[0],index_mob[0],seq_ref[self.chain1],seq_mob[self.mobchain1])
        pdb2alignref[self.chain2],pdb2alignmob[self.mobchain2],ref2mob[self.chain2],mob2ref[self.mobchain2] = getcorrespondance(index_ref[1],index_mob[1],seq_ref[self.chain2],seq_mob[self.mobchain2])
        #pdb2alignref["A"],pdb2alignmob["A"],ref2mob["A"],mob2ref["A"] = getcorrespondance(index_ref[0],index_mob[0],seq_ref["A"],seq_mob["A"])
        #pdb2alignref["B"],pdb2alignmob["B"],ref2mob["B"],mob2ref["B"] = getcorrespondance(index_ref[1],index_mob[1],seq_ref["B"],seq_mob["B"])
        
        self.pdb2alignref = pdb2alignref
        self.pdb2alignmob = pdb2alignmob
        self.ref2mob      = ref2mob
        self.mob2ref      = mob2ref

        # keep for visuAlignCAPRI
        self.align_seq = align_seq
        
        return 
    
    
    def alignmentTraitementCAPRI(self):
        """ Alignement process if the structures are different"""
        
        f_ali = open(self.alignment_file,"r")
        l_ali = f_ali.readlines()[2:]
        f_ali.close()
        
        
#         data = "".join([i.strip()+"\n" for i in l_ali[2:]]).replace("\n\n\n","\n").replace("\n\n","\n").split()
        data = [i.strip() for i in  l_ali if len(i.strip())>0]
        seq = None
        # Parsing the sequences
        sequences          = {}
        sequences["order"] = []
        nbSeq             = 0
        cur_id = None
        for line in data:
            if line[0] == ">":
                if cur_id != None:
                    sequences[cur_id] = seq
                    sequences["order"].append(cur_id)
                
                header = line.strip()
                cur_id = header
                
                nbSeq += 1
                seq = "" 
            else:
                seq += line.strip()
        sequences[cur_id] = seq
        sequences["order"].append(cur_id)
        
        
        # Sort the sequence
        seq = {}
        seq["order"] =  []
        for header in sequences["order"]:
            current_id = os.path.basename(header.split()[0][1:])
            if not current_id in seq.keys():
                seq[current_id] = {}
                seq["order"].append(current_id)
            
            chain = header.split("Chain")[1].replace("'","").strip()
            seq[current_id][chain] = sequences[header]
            
        
        return seq
        
    def getRankingCAPRI(self):
        """
        Ranking feature extract from : Mendez,Leplae,Lensing,Wodak PROTEINS 2005
        """
        rank = ""
        if self.dicvar["lrms"] and self.dicvar["irms"] and self.dicvar["fraction"]:
            if self.fractionNat >= 50.0 and (self.LrmsA<=1.0 or self.Irms<=1.0):
                rank = "High"
            if ((self.fractionNat >= 30.0 and self.fractionNat < 50.0) and (self.LrmsA<=5.0 or self.Irms<=2.0)) or (self.fractionNat >= 50.0 and self.LrmsA>1.0 and self.Irms>1.0):
                rank = "Medium"
            if ((self.fractionNat >= 10.0 and self.fractionNat < 30.0) and (self.LrmsA<=10.0 or self.Irms<=4.0)) or (self.fractionNat >= 30.0 and self.LrmsA>5.0 and self.Irms>2.0):
                rank = "Acceptable"
            if self.fractionNat<10.0 or (self.LrmsA>10.0 and self.Irms>4.0):
                rank = "Incorrect" 
        else:
            rank = "No Ranking"        
        
        return rank
    
    def getDockQ(self):
        """
        Calculate DockQ score based on formula in DockQ - Basu & Wallner 2016
        """
        if self.dicvar["lrms"] and self.dicvar["irms"] and self.dicvar["fraction"]:
            self.DockQ=(float(self.fractionNat)/100. + 1/(1+(self.Irms/1.5)*(self.Irms/1.5)) + 1/(1+(self.LrmsA/8.5)*(self.LrmsA/8.5)))/3
            rank = ""
            if self.DockQ >= 0.80:
                rank = "High"
            if 0.80 > self.DockQ >= 0.49:
                rank = "Medium"
            if 0.49 > self.DockQ >= 0.23:
                rank = "Acceptable"
            if 0.23 > self.DockQ:
                rank = "Incorrect" 
        else:
            rank = "No Ranking"        
        
        return rank        
        
    def screenerCAPRI(self,mobile,align=False):
        """ """

        if self.dicvar["verbose"]:
            print(os.path.basename(mobile))
            if self.dicvar["lrms"]:
                print("RMS AvsA               :     %5.2f  "%(self.rmsAvsA))
                print("RMS BvsB               :     %5.2f  "%(self.rmsBvsB))
                print("RMS ABvsAB             :     %5.2f  "%(self.rmsABvsAB))
                print("LRMS A constrains      :     %5.2f  "%(self.LrmsA))
                print("LRMS B constrains      :     %5.2f  "%(self.LrmsB))
            
            if self.dicvar["fraction"]:
                print("Fraction Native        :     %5.2f   "%(self.fractionNat))
                print("Fraction Non Native    :     %5.2f   "%(self.fractionNonNat))
                print("Fraction InterfaceRes  :     %5.2f   "%(self.fractionIR))
            
            if self.dicvar["irms"]:
                print("IRMS                   :     %5.3f   "%(self.Irms))
                
            if self.dicvar["lrms"] and self.dicvar["irms"] and self.dicvar["fraction"]:    
                print("CAPRI                  :     %s"%self.rank)
                print("DockQ                  :     %s"%self.DockQ)
                print("DockQ-rank             :     %s"%self.rankQ)

        
            print("-----------------------------------------------------------------------")
            return
        
        
        screener = ""
        screener += "%s "%(mobile)
        if self.dicvar["lrms"]:
            screener += "%5.2f %5.2f %5.2f %5.2f %5.2f "%(self.rmsAvsA,self.rmsBvsB,self.rmsABvsAB,self.LrmsA,self.LrmsB)
                
        if self.dicvar["fraction"]:
            screener += "%5.2f %5.2f %5.2f "%(self.fractionNat,self.fractionNonNat,self.fractionIR)
            
        if self.dicvar["irms"]:
            screener += "%5.2f "%(self.Irms)
            
        if self.dicvar["lrms"] and self.dicvar["irms"] and self.dicvar["fraction"]:    
            screener += "%s "%self.rank
            screener += "%.2f "%self.DockQ
            screener += "%s "%self.rankQ
        screener +="\n"
        
        print(screener)
        if self.dicvar["file"]:
            self.file.write(screener)
        
        return
        
        
if __name__=='__main__':
      
    dicvar = {}
    
    dicvar["dump"] = True
    try:
        dicvar["reference"] = sys.argv[sys.argv.index("-ref")+1]
        dicvar["chain"] = sys.argv[sys.argv.index("-chain")+1]
    except:
        usage()

    try:
        dicvar["liste_mob"] = sys.argv[sys.argv.index("-lm")+1]
    except:
        pass
    try:
        dicvar["mobile"] = sys.argv[sys.argv.index("-mob")+1]
    except:
        pass

    try:
        test = sys.argv[sys.argv.index("-quick_fnat")]
        dicvar["quick_fnat"] = True
    except:
        pass 

    try:
        dicvar["quick_fnat_ref"] = sys.argv[sys.argv.index("-quick_fnat_ref")+1]
    except:
        pass 

    try:
        test = sys.argv[sys.argv.index("-lrms")]
        dicvar["lrms"] = True
    except:
        pass
    
    try:
        test = sys.argv[sys.argv.index("-fr")]
        dicvar["fraction"] = True
    except:
        pass

    try:
        dicvar["file"] = sys.argv[sys.argv.index("-f")+1]
    except:
        pass
    
    try:
        test = sys.argv[sys.argv.index("-irms")]
        dicvar["irms"] = True
    except:
        pass
    
    try:
        test = sys.argv[sys.argv.index("-v")]
        dicvar["verbose"] = True
    except:
        pass
        
    
    
    CAPRI(dicvar)
        
