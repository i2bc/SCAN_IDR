import re, os, string, sys, glob

lst_AA = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG",
          "SER", "THR", "VAL", "TRP", "TYR", "5HP", "ABA", "PCA", "FGL", "BHD", "HTR", "MSE", "CEA", "ALS", "TRO",
          "TPQ", "MHO", "IAS", "HYP", "CGU", "CSE", "RON", "3GA","TYS","AYA", "FME", "CXM", "SAC", "CSO", "MME",
          "SEG", "HSE", "HSD","HSP", 'PAQ', 'AGM', 'PR3', 'DOH', 'CCS', 'GSC', 'GHG', 'OAS', 'MIS', 'SIN', 'TPL',
          'SAC', '4HT', 'FGP', 'HSO', 'LYZ', 'FGL', 'PRS', 'DCY', 'LYM', 'GPL', 'PYX', 'PCC', 'EHP', 'CHG', 'TPO',
          'DAS', 'AYA', 'TYN', 'SVA', 'SCY', 'BNN', '5HP', 'HAR', 'IAS', 'SNC', 'AHB', 'PTR', 'PHI', 'NPH', 'PHL',
          'SNN', 'A66', 'TYB', 'PHD', 'MAA', 'APN', 'TYY', 'TYT', 'TIH', 'TRG', 'CXM', 'DIV', 'TYS', 'DTH', 'MLE',
          'CME', 'SHR', 'OCY', 'DTY', '2AS', 'AEI', 'DTR', 'OCS', 'CMT', 'BET', 'NLP', 'LLY', 'SCH', 'CEA', 'LLP',
          'TRF', 'HMR', 'TYI', 'TRO', 'NLE', 'BMT', 'BUC', 'PEC', 'BUG', 'SCS', 'NLN', 'MHO', 'CSO', 'FTR', 'DLE',
          'TRN', 'CSE', 'CSD', 'OMT', 'CSA', 'DSP', 'CSB', 'DSN', 'SHC', 'CSX', 'YCM', 'CSZ', 'TRQ', 'CSW', 'EFC',
          'CSP', 'CSS', 'CSR', 'CZZ', 'MSO', 'BTR', 'HLU', 'MGN', 'HTI', 'TYQ', '4IN', 'M3L', 'C5C', 'HTR', 'MPQ',
          'KCX', 'GLH', 'DIL', 'ACA', 'NEM', '5CS', 'LYX', 'DVA', 'ACL', 'GLX', 'MLZ', 'GLZ', 'SME', 'SMC', 'DLY',
          'NEP', 'BCS', 'ASQ', 'SET', 'SEP', 'ASX', 'DGN', 'DGL', 'MHS', 'SEG', 'ASB', 'ASA', 'SEC', 'SEB', 'ASK',
          'GGL', 'ASI', 'SEL', 'CGU', 'C6C', 'ASL', 'LTR', 'CLD', 'CLE', 'GMA', '1LU', 'CLB', 'MVA', 'S1H', 'DNP',
          'SAR', 'FME', 'ALO', 'ALM', 'LEF', 'MEN', 'TPQ', 'NMC', 'SBD', 'ALY', 'MME', 'GL3', 'ALS', 'SBL', '2MR',
          'CAY', '3AH', 'DPR', 'CAS', 'NC1', 'HYP', 'FLA', 'LCX', 'MSE', 'IYR', 'DPN', 'BAL', 'CAF', 'MSA', 'AIB',
          'HIP', 'CYQ', 'PCA', 'DAL', 'BFD', 'DAH', 'HIC', 'CYG', 'DAR', 'CYD', 'IIL', 'CYM', 'CYL', 'CY3', 'CY1',
          'HAC', '143', 'DHI', 'CY4', 'YOF', 'HPQ', 'SOC', 'DHA', '2LU', 'MLY', 'TRW', 'STY', 'MCL', 'BHD', 'NRQ',
          'ARM', 'PRR', 'ARO', "5HP", "ABA", "PCA", "FGL", "BHD", "HTR", "MSE", "CEA", "ALS", "TRO", "TPQ", "MHO",
          "IAS", "HYP", "CGU", "CSE", "RON", "3GA", "TYS", "AYA", "FME", "CXM", "SAC", "CSO", "MME", "SEG", "HSE",
          "HSD", "HSP"]


def usage():
    print("USAGE : RenumberPDB.py -s <pdb> -r1 <chain:index first residue> (example : A:1)")
    print("                            [ -o <pdbout> (default = <pdb>_ren.pdb)]")
    print("                            [ -force : renumber forcing residue numbers sequentiality")
    print("                            [ -force_chain : renumber forcing residue numbers sequentiality through chains")
    print("")
    sys.exit()


class RenumberPDB:
    """
    RenumberPDB.RenumberPDB(pdbin, firstres, pdbout, force_sequential = False )
        -> dump a renumbered output pdb file

        ### INPUT DESCRIPTION
        - pdbin : <input pdbfile>
        - firstres <chain:index first residue> (example : A:1)
        - force_sequential : renumber forcing residue numbers sequentiality
        -force_chain : renumber forcing residue numbers sequentiality through chains
        ### COMMENTS
        - Does not renumber HETEROATOMS (TODO ?)
    """

    def __init__(self, pdbin, firstres, pdbout, force_sequential=False, force_sequential_overchains=False):
        self.pdbin = pdbin
        self.firstres = firstres
        self.pdbout = pdbout
        if force_sequential:
            self.forceseq = "YES"
        else:
            self.forceseq = "NO"
        if force_sequential_overchains:
            self.forceseq_overchains = "YES"
        else:
            self.forceseq_overchains = "NO"
        self.parse()
        self.renumb()

    def parse(self):
        (self.chain2mod, self.firstres2mod) = self.firstres.split(":")
        if self.chain2mod == "":
            self.chain2mod = "_"

    def renumb(self):

        f = open(self.pdbin)
        fin = f.readlines()
        f.close()
        fout = open(self.pdbout, "w")
        oldres = -10000
        oldrestype = ""
        oldextra = ""  # for AA insertion in PDB i.e. consecutive residues with same number
        diff = 10000
        if self.forceseq == "YES":
            residue_counter = int(self.firstres2mod) - 1
        GET_1stRES = 0
        for line in fin:
            split = line.split()
            if (len(split) > 4 and line.startswith('ATOM')) or (
                    len(split) > 4 and line.startswith('HETATM') and line[17:20].strip() in lst_AA):
                try:  # no chain in the pdb, residue index is found in pos 4
                    testchain = int(line[21])
                    chain = '_'
                    resn = line[22:26]
                    extra = line[26].strip()
                    restype = line[17:20].strip()
                except ValueError:
                    chain = line[21]
                    resn = line[22:26]
                    extra = line[26].strip()
                    restype = line[17:20].strip()

                if chain == self.chain2mod and GET_1stRES == 0:
                    GET_1stRES = 1  # Flag to capture the numbering
                    # shift only once when chain2mod is 1st encountered
                    diff = int(resn) - int(self.firstres2mod)
                if self.forceseq_overchains:
                    self.chain2mod = chain
                if chain == self.chain2mod:
                    if self.forceseq == "YES":
                        if (resn) != oldres or extra != oldextra or restype != oldrestype:
                            # if residue insertion (i.e. 15A then 15B), add 1 to counter;
                            # if residue number doesn't change but residue type changes, add 1 to counter
                            residue_counter += 1
                            oldres = (resn)
                            oldextra = extra
                            oldrestype = restype
                        newresn = "%4d" % (residue_counter)
                    else:
                        if (resn) != oldres:
                            oldres = resn
                            oldextra = extra
                            oldrestype = restype
                        elif extra != oldextra:
                            oldextra = extra
                            oldrestype = restype
                            diff -= 1
                        elif restype != oldrestype:
                            oldrestype = restype
                            diff -= 1
                        newresn = "%4d" % (int(resn) - diff)
                    newline = line[:22] + newresn + " " + line[27:]
                else:
                    newline = line
            else:
                newline = line
            fout.write(newline)
        fout.close()


## Create template pdb for loop calculation in rosetta
if __name__ == '__main__':
    ## Read input self.parameters
    if len(sys.argv) > 1:
        try:
            pdbin = sys.argv[sys.argv.index('-s') + 1]
        except:
            usage()
        try:
            firstres = sys.argv[sys.argv.index('-r1') + 1]
        except:
            usage()
        try:
            force = sys.argv[sys.argv.index('-force')]
        except:
            force = None
        try:
            force_chain = sys.argv[sys.argv.index('-force_chain')]
        except:
            force_chain = None
        try:
            pdbout = sys.argv[sys.argv.index('-o') + 1]
        except:
            pdbout = ".".join(pdbin.split(".")[:-1]) + "_ren.pdb"
    else:
        usage()

    RP = RenumberPDB(pdbin, firstres, pdbout, force_sequential=force, force_sequential_overchains=force_chain)
