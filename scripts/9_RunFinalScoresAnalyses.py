import os, re, copy, glob
import random
import configparser

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)

config = configparser.ConfigParser()
config.read(f'{script_dir}/config.ini')

WORKING_DIR = config['DEFAULT']['WORKING_DIR']
SCRIPTS_DIR = os.path.abspath(config['DEFAULT']['SCRIPTS_DIR'])
DATA_DIR = os.path.abspath(os.path.join(WORKING_DIR, 'data', ))


CUTMODELS_DIR = os.path.abspath(os.path.join(DATA_DIR, config['DEFAULT']['CUTMODELS_DIR']))
CAPRIEVAL_DIR = os.path.abspath(os.path.join(DATA_DIR, config['DEFAULT']['CAPRIEVAL_DIR']))
REFERENCE_STRUCT_DIR = os.path.abspath(os.path.join(DATA_DIR, config['DEFAULT']['REFERENCE_DIR']))

n_iterations = int(config['DEFAULT']['N_ITERATIONS'])

RESULTS_DIR = os.path.abspath(os.path.join(DATA_DIR, config['DEFAULT']['RESULTS_DIR']))
fglobaloutput_path = os.path.abspath(os.path.join(RESULTS_DIR, config['DEFAULT']['OUTPUT_GLOBAL']))

"""
This script is designed to analyze the output of the AF2 runs on the peptides
It prepares the data into tables in a result file directory that can be processed 
subsequently by the Jupyter notebook to perform the Figure display
"""

class AnalyseResults:

    def __init__(self,pathfile) -> None:
        self.Ncol = 23
        self.pathfile = pathfile
        self.dg     = dict()    # global dict storing all the data
        self.ditems = dict()    # d[item] = col_index
        self.dmodel = dict()    # d[f"{entry}_{protocol}_{version}_{model}"] = pdb_filepath
        self.pathout = DATA_DIR
        self.extension = ""
        self.list_exclude = []
        # some entries can be excluded from the analysis adding lines as self.list_exclude.append("<index>_<pdb>")

        #
        # Creation of a dictionary listing the protocols
        #
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
        d_protocol = dict()
        for proto in list_of_possible_protocols:
            d_protocol[proto] = dict()
        #
        # -0- We analyze each protocol individually
        #
        # Defines self.list_protocols with the series of protocols to consider

        for ii, protocol in enumerate(list_of_possible_protocols):
            print("****", protocol, "****")
            # Put all the results for the file into a dictionary which can be subsequently analyzed
            self.parse_results(list_protocols=[protocol])
            # Write the information of the best models for every protocol in a single file 
            if ii == 0: writing_option = "w"
            else:  writing_option = "a"
            self.write_best_models_info(list_protocols=[protocol], fileoutbest=f"Infos_bestModels{self.extension}.out", write_or_append=writing_option)
            # Remove the pdbs which are not to be considered for the analysis
            if len(self.list_exclude) > 0:
                self.exclude_pdb()
            # Statistics over the 5 iterations
            self.analyse_dic_bestCombScore_vs_bestDockQ(fileout_name=f"Stats_bestCombScore_vs_bestDockQ_{protocol}{self.extension}.out")
            # Writting an output file extracting the rank_1 results of every version in a limited set of protocols
            self.analyse_dic_iterations(fileout_rootname=f"Stats_5iterations{self.extension}", list_protocols=[protocol])

        #
        # -1- We analyze the results focusing on three protocols
        #
        print("****", "3 delim-delim protocols combined analysis", "****")
        # Defines self.list_protocols with the series of protocols to consider
        l_protocols = ['mixed_ali-delim-delim', 'unpaired_ali-delim-delim', 'unpaired_ali-delim-delim',]
        # Put all the results for the file into a dictionary which can be subsequently analyzed
        self.parse_results(list_protocols=l_protocols)
        # Remove the pdbs which are not correct anymore
        self.exclude_pdb()
        # Statistics over the 5 interations
        self.analyse_dic_bestCombScore_vs_bestDockQ(fileout_name=f"Stats_bestCombScore_vs_bestDockQ_3scores{self.extension}.out")
        # Writting an output file extracting the rank_1 results of every version in a limited set of protocols
        self.analyse_dic_iterations(fileout_rootname=f"Stats_5iterations{self.extension}", list_protocols=l_protocols)

        #
        # -2- We analyze the results focusing on four protocols 'mixed_ali-delim-delim', 'unpaired_ali-delim-delim', 'unpaired_ali-delim-delim', 'mixed_ali-fl-fl'
        #
        print("****", "3 delim-delim + 1 fl-fl protocols combined analysis", "****")
        # Defines self.list_protocols with the series of protocols to consider
        l_protocols = ['mixed_ali-delim-delim', 'unpaired_ali-delim-delim', 'unpaired_ali-delim-delim', 'mixed_ali-fl-fl']
        # Put all the results for the file into a dictionary which can be subsequently analyzed
        self.parse_results(list_protocols=l_protocols)
        # Remove the pdbs which are not correct anymore
        self.exclude_pdb()
        # Statistics over the 5 iterations
        self.analyse_dic_bestCombScore_vs_bestDockQ(fileout_name=f"Stats_bestCombScore_vs_bestDockQ_4scoresFL{self.extension}.out")
        # Writting an output file extracting the rank_1 results of every version in a limited set of protocols
        self.analyse_dic_iterations(fileout_rootname=f"Stats_5iterations{self.extension}", list_protocols=l_protocols)

        #
        # -3- We analyze the results focusing on two protocols 'mixed_ali-delim-delim' and 'mixed_ali-delim-200'
        #
        print("****", "1 delim-delim + 1 delim-200 protocols comparative analysis", "****")
        # Defines self.list_protocols with the series of protocols to consider
        l_protocols = ['mixed_ali-delim-delim', 'mixed_ali-delim-200',]
        # Put all the results for the file into a dictionary which can be subsequently analyzed
        self.parse_results(list_protocols=l_protocols)
        # Remove the pdbs which are not correct anymore
        self.exclude_pdb()
        # Statistics over the 5 iterations
        self.analyse_dic_bestCombScore_vs_bestDockQ(fileout_name=f"Stats_bestCombScore_vs_bestDockQ_2scores{self.extension}.out")
        # Writting an output file extracting the rank_1 results of every version in a limited set of protocols
        self.analyse_dic_iterations(fileout_rootname=f"Stats_5iterations{self.extension}", list_protocols=l_protocols)
        # Analysis of the cases successful at short length but not at length ~ 200
        self.compare_diff_cases(protocol1='mixed_ali-delim-delim', protocol2='mixed_ali-delim-200', fileout_rootname=f"CompareDiffSuccess{self.extension}")

        #
        # -4- We analyze the results focusing on four protocols 'mixed_ali-delim-delim', 'unpaired_ali-delim-delim', 'unpaired_ali-delim-delim', 'mixed_ali-delim-100',
        #
        print("****", "3 delim-delim + 1 delim-100 protocols sub-sampling analysis", "****")
        # Defines self.list_protocols with the series of protocols to consider
        l_protocols = ['mixed_ali-delim-delim', 'unpaired_ali-delim-delim', 'unpaired_ali-delim-delim', 'mixed_ali-delim-100',]
        # Put all the results for the file into a dictionary which can be subsequently analyzed
        self.parse_results(list_protocols=l_protocols)
        # Remove the pdbs which are not correct anymore
        self.exclude_pdb()
        # Statistics over the 5 iterations
        self.analyse_dic_bestCombScore_vs_bestDockQ(fileout_name=f"Stats_bestCombScore_vs_bestDockQ_4scores{self.extension}.out")
        # Writting an output file extracting the rank_1 results of every version in a limited set of protocols
        self.analyse_dic_iterations(fileout_rootname=f"Stats_5iterations{self.extension}", list_protocols=l_protocols)
        # We select 1 out of 5 version of every protocol below. Sample the results Ntimes and extract a mean success rate
        self.sample_protocols(list_protocols=['mixed_ali-delim-delim', 'mixed_ali-delim-delim', 'unpaired_ali-delim-delim', 'unpaired_ali-delim-delim', 'mixed_ali-delim-100'], repeat=100)


    def get_pdb_filepath(self, protocol, pdb, index, version, rank, model, dirpath='cut4capri'):
        """
        Recovers the path where the pdb can be found
        :return: pdb_filepath
        """
        if dirpath == 'cut4capri':
            root_path = CUTMODELS_DIR
            dir1 = index + '_' + pdb
            dir2 = f"af2_{dir1}_{protocol}"
            dir3 = f"MSA4af2_{dir1}_{protocol}_unrelaxed_rank_{rank}_model_{model}_v{version}.pdb"
            pdb_filepath = os.path.join(root_path, dir1, dir2, dir3)
        if dirpath == 'capriref':
            root_path = REFERENCE_STRUCT_DIR
            pdb = pdb.upper()
            l_pdb = glob.glob(f"{pdb}*pdb")
            pdb_filepath = []
            for f in l_pdb:
                pdb_filepath.append(os.path.join(root_path, f)) # result is a list
        return pdb_filepath

    def write_best_models_info(self, list_protocols=[], fileoutbest=None, write_or_append="w"):
        """
        After parsing the self.dg dictionnary contains the information for the best model in a set of protocols
        This information will be dumped for that model for every pdb entry
        """
        fbest_path = os.path.abspath(os.path.join(RESULTS_DIR, fileoutbest))
        fout = open(fbest_path, write_or_append)
        if write_or_append =="w":
            fout.write('# $1:index_pdb\n# $2:sorting_method\n# $3:protocol\n# $4:SizeModel\n# $5:CombinedSCORE\n# $6:DockQ\n# $7:CAPRIrankPeptide\n# $8:IRMS\n# $9:filepath\n')
        str_protocol = "_".join(list_protocols)
        list_best_scores = [
                            'sort_by_CombinedSCORE',
                            'sort_by_DockQ', 
                            'sort_by_IRMS',
                            ]
        for p in self.dg:
            for best_score in list_best_scores:
                sorting_method = '_'.join(best_score.split('_')[-2:])
                d_temp = self.dg[p][best_score][0][1]
                best_comb = d_temp['CombinedSCORE']
                if str_protocol in self.dg[p]:
                    best_size = self.dg[p][str_protocol]["SizeModel"]
                else:
                    best_size = str(-1)
                best_rcapri = d_temp['CAPRIrankPeptide']
                best_DockQ = d_temp['DockQ']
                best_irms = d_temp['IRMS']
                best_path = d_temp['path']
                fout.write(f"{p}\tbest_{sorting_method}\t{str_protocol}\t{best_size}\t{best_comb}\t{best_DockQ}\t{best_rcapri}\t{best_irms}\t{best_path}\n")
        fout.close()
            

    def parse_results(self, list_protocols=[]):
        """
        Creates a self.dg based on all protocols or on those listed in [list_protocols]
        self.dg[index_pdb][protocol]["SizeModel"] -> size of the models
        self.dg[index_pdb][protocol]["Reference"] -> path to the reference
        self.dg[index_pdb][protocol]["sort_by_IRMS"] -> sorted list [sorted_model_index(0 to 5)][dict[model_properties]]
        self.dg[index_pdb][protocol]["sort_by_DockQ"] -> sorted list [sorted_model_index(0 to 5)][dict[model_properties]]
        self.dg[index_pdb][protocol]["sort_by_CombinedSCORE"] -> sorted list  [sorted_model_index(0 to 5)][dict[model_properties]]
        """
        self.dg = dict()
        self.ditems = dict()
        self.dmodel = dict()
        # Analyze the commented lines at the beginning of the Result file to recover the column headers
        p1 = re.compile(r"\$(\d+):(\S+)")
        with open(self.pathfile) as f:
            for l in f.readlines()[:self.Ncol]:
                s = l.split()
                m = p1.match(s[1])
                index = int(m.group(1))
                label = m.group(2)
                self.ditems[label] = index-1
                #print(f"dtemp['{label}'] = s[{index-1}]")
        with open(self.pathfile) as f:
            for l in f.readlines()[self.Ncol:]:
                s = l.split()
                ipdb = s[0]
                pdb = ipdb.split("_")[1]
                index = ipdb.split("_")[0]
                protocol = s[1]
                size = s[2]
                version = s[3]
                model = s[4]
                AF2Rank = s[5]

                # If a list of protocols is provided the sorting of the solutions will be restricted to that set
                if len(list_protocols) > 0 and protocol not in list_protocols:
                    continue

                # We build the reference filepath name for every entry,
                try:
                    int(AF2Rank)
                    pdb_filepath = self.get_pdb_filepath(protocol, pdb, index, version, AF2Rank, model, dirpath='cut4capri')
                    self.dmodel[f"{pdb}_{protocol}_{version}_{model}"] = pdb_filepath
                    #print(pdb_filepath)
                except: # We are considering a consensus reported in the output file (deprecated since we do not compute this best anymore)
                    pdb_filepath = self.dmodel[f"{pdb}_{protocol}_{version}_{model}"]
                    protocol = AF2Rank
                dtemp = dict()
                dtemp['path'] = pdb_filepath
                dtemp['pLDDT'] = s[6]
                dtemp['pTMSCORE'] = s[7]
                dtemp['ipTMSCORE'] = s[8]
                dtemp['CombinedSCORE'] = s[9]
                dtemp['DockQ'] = s[10]
                dtemp['DockQrank'] = s[11]
                dtemp['CAPRIrank'] = s[12]
                dtemp['CAPRIrankPeptide'] = s[13]
                dtemp['IRMS'] = s[14]
                dtemp['RMSAvsA'] = s[15]
                dtemp['RMSBvsB'] = s[16]
                dtemp['RMSABvsAB'] = s[17]
                dtemp['LRMSA'] = s[18]
                dtemp['LRMSB'] = s[19]
                dtemp['FRNAT'] = s[20]
                dtemp['FRNNAT'] = s[21]
                dtemp['FRIR'] = s[22]
                if ipdb not in self.dg:
                    self.dg[ipdb] = dict()
                    self.dg[ipdb]["SizeModel"] = size
                    self.dg[ipdb]['sort_by_CombinedSCORE'] = []
                    self.dg[ipdb]['sort_by_DockQ'] = []
                    self.dg[ipdb]['sort_by_IRMS'] = []
                if protocol not in self.dg[ipdb]:
                    self.dg[ipdb][protocol] = dict()
                    self.dg[ipdb][protocol]["SizeModel"] = size
                    self.dg[ipdb][protocol]["Reference"] = self.get_pdb_filepath(protocol, pdb, index, version,
                                                                             AF2Rank, model, dirpath='capriref')
                    self.dg[ipdb][protocol]["sort_by_IRMS"] = []
                    self.dg[ipdb][protocol]["sort_by_DockQ"] = []
                    self.dg[ipdb][protocol]["sort_by_CombinedSCORE"] = []
                if version not in self.dg[ipdb][protocol]:
                    self.dg[ipdb][protocol][version] = dict()
                    self.dg[ipdb][protocol][version]["sort_by_IRMS"] = []
                    self.dg[ipdb][protocol][version]["sort_by_DockQ"] = []
                    self.dg[ipdb][protocol][version]["sort_by_CombinedSCORE"] = []
                # To sort and find best model for all protocols requested for a given pdb case, all protocols and iterations combined
                self.dg[ipdb]['sort_by_IRMS'].append([eval(dtemp['IRMS']), copy.deepcopy(dtemp)])
                self.dg[ipdb]['sort_by_DockQ'].append([eval(dtemp['DockQ']), copy.deepcopy(dtemp)])
                self.dg[ipdb]['sort_by_CombinedSCORE'].append([eval(dtemp['CombinedSCORE']), copy.deepcopy(dtemp)])
                # To sort and find best model for every protocol, all iterations combined
                self.dg[ipdb][protocol]["sort_by_IRMS"].append([eval(dtemp['IRMS']), copy.deepcopy(dtemp)])
                self.dg[ipdb][protocol]["sort_by_DockQ"].append([eval(dtemp['DockQ']), copy.deepcopy(dtemp)])
                self.dg[ipdb][protocol]["sort_by_CombinedSCORE"].append([eval(dtemp['CombinedSCORE']), copy.deepcopy(dtemp)])
                # To sort and find best model for every iteration 
                self.dg[ipdb][protocol][version]["sort_by_IRMS"].append([eval(dtemp['IRMS']), copy.deepcopy(dtemp)])
                self.dg[ipdb][protocol][version]["sort_by_DockQ"].append([eval(dtemp['DockQ']), copy.deepcopy(dtemp)])
                self.dg[ipdb][protocol][version]["sort_by_CombinedSCORE"].append([eval(dtemp['CombinedSCORE']), copy.deepcopy(dtemp)])

        for p in self.dg:
            # -1- We compute and store the best solution over several protocols+iterations
            self.dg[p]["sort_by_IRMS"].sort(key = lambda x: x[0])
            self.dg[p]["sort_by_DockQ"].sort(key = lambda x: x[0],reverse=True)
            self.dg[p]["sort_by_CombinedSCORE"].sort(key = lambda x: x[0],reverse=True)
            for c in list_protocols:
                # -2- We compute and store the best solution over several iterations
                self.dg[p][c]["sort_by_IRMS"].sort(key = lambda x: x[0])
                self.dg[p][c]["sort_by_DockQ"].sort(key = lambda x: x[0],reverse=True)
                self.dg[p][c]["sort_by_CombinedSCORE"].sort(key = lambda x: x[0],reverse=True)
                # -3- We store the best solution for every version
                for it in range(n_iterations):
                    v = str(it + 1)
                    if c in ["best_CombinedSCORE", "best_iSCORE","best_DockQ",
                             "best_iSCORE_of_all_protocols","best_combSCORE_of_all_protocols","best_DockQ_of_all_protocols"]:
                        continue
                    self.dg[p][c][v]["sort_by_IRMS"].sort(key=lambda x: x[0])
                    self.dg[p][c][v]["sort_by_DockQ"].sort(key=lambda x: x[0], reverse=True)
                    self.dg[p][c][v]["sort_by_CombinedSCORE"].sort(key=lambda x: x[0], reverse=True)


    def exclude_pdb(self):
        """
        Used to remove specific pdb from the analysis
        """
        for p in self.list_exclude:
            del self.dg[p]

    def analyse_dic_bestCombScore_vs_bestDockQ(self, scoringScheme1='sort_by_CombinedSCORE', scoringScheme2='sort_by_DockQ', fileout_name=None):

        counter_score = 0
        counter_dockQ = 0
        d_missed = dict()

        if not fileout_name:
            fileout = os.path.abspath(os.path.join(RESULTS_DIR, f"Stats_{scoringScheme1}_vs_{scoringScheme2}.out"))
        else:
            fileout = os.path.abspath(os.path.join(RESULTS_DIR, fileout_name))

        with open(os.path.join(fileout), "w") as fout:
            d_missed[scoringScheme1] = []
            d_missed[scoringScheme2] = []
            for p in self.dg:
                print(f"#### {p} ####")
                #print(f"Best Combined Score")
                best_comb = self.dg[p][scoringScheme1][0][1]['CombinedSCORE']
                best_comb_rcapri = self.dg[p][scoringScheme1][0][1]['CAPRIrankPeptide']
                best_comb_DockQ = self.dg[p][scoringScheme1][0][1]['DockQ']
                best_comb_path = self.dg[p][scoringScheme1][0][1]['path']
                if best_comb_rcapri in ['High', 'Medium', 'Acceptable']:
                    counter_score += 1
                else:
                    d_missed[scoringScheme1].append(p)
                #print(f"Best DockQ")
                best_dockq = self.dg[p][scoringScheme2][0][1]['CombinedSCORE']
                best_dockq_rcapri = self.dg[p][scoringScheme2][0][1]['CAPRIrankPeptide']
                best_dockq_DockQ = self.dg[p][scoringScheme2][0][1]['DockQ']
                best_dockq_path = self.dg[p][scoringScheme2][0][1]['path']
                if best_dockq_rcapri in ['High', 'Medium', 'Acceptable']:
                    counter_dockQ += 1
                else:
                    d_missed[scoringScheme2].append(p)
                fout.write(f"#CombScore: {p}\t{best_comb_rcapri}\t{best_comb}\t{best_comb_DockQ}\t{best_comb_path}\n")
                fout.write(f"#DockQScore: {p}\t{best_dockq_rcapri}\t{best_dockq}\t{best_dockq_DockQ}\t{best_dockq_path}\n")
            fout.write(f"#SuccessRateScore: {counter_score}\t#SuccessRateDockQ: {counter_dockQ}\t#TotalCases: {len(self.dg)}\n")
            for protocol in d_missed:
                for p in d_missed[protocol]:
                    fout.write(f"#Missed_{protocol}: {p}\n")
            #print(counter_score)
            #print(counter_dockQ)
            #print(len(self.dg))

    def compare_diff_cases(self, protocol1='mixed_ali-delim-delim',  protocol2='mixed_ali-delim-200', fileout_rootname=None):
        """
        Returns the cases with different performance between protocol1 and protocol2
        """

        ####
        str_protocol1 = '-'.join([protocol1])
        str_protocol2 = '-'.join([protocol2])

        if not fileout_rootname:
            fileout = os.path.abspath(os.path.join(RESULTS_DIR, f"Compare_Cond1_{str_protocol1}_Cond2_{str_protocol2}.out"))
        else:
            fileout = os.path.abspath(os.path.join(RESULTS_DIR, f"{fileout_rootname}_Cond1_{str_protocol1}_Cond2_{str_protocol2}.out"))

        with open(os.path.join(self.pathout, fileout), "w") as fout:
            fout.write(f"# $1:index_pdb\n# $2:protocol\n# $3:CAPRIrankPeptide\n# $4:CombinedSCORE\n# $5:DockQ\n# $6:IRMS\n# $7:Path2File\n")
            fout.write(f"# TotalCases: {len(self.dg)}\n")
            for p in self.dg:
                d_test_2protocols = {}
                d_best_comb_rcapri = {}
                d_best_comb = {}
                d_best_comb_DockQ = {}
                d_best_comb_irms = {}
                d_best_comb_path = {}
                for c in [protocol1, protocol2]:
                    counter_score = 0
                    d_test_2protocols[c] = False
                    DockQ_ref = 0.0
                    for it in range(n_iterations):
                        v = str(it + 1)
                        try:
                            d_best_model_AF2 = self.dg[p][c][v]['sort_by_CombinedSCORE'][0][1]

                            if eval(d_best_model_AF2['DockQ']) > DockQ_ref:
                                d_best_comb_rcapri[c] = d_best_model_AF2['CAPRIrankPeptide']
                                d_best_comb[c] = d_best_model_AF2['CombinedSCORE']
                                d_best_comb_DockQ[c] = d_best_model_AF2['DockQ']
                                d_best_comb_irms[c] = d_best_model_AF2['IRMS']
                                d_best_comb_path[c] = d_best_model_AF2['path']
                                DockQ_ref = eval(d_best_comb_DockQ[c])
                            if d_best_model_AF2['CAPRIrankPeptide'] in ['High', 'Medium', 'Acceptable']:
                                counter_score += 1
                        except KeyError:
                            #print("KeyError",c,v)
                            continue
                    if counter_score > 0:
                        d_test_2protocols[c] = True
                l_results = list(d_test_2protocols.items())
                #print(l_results)
                #print(d_best_comb_rcapri)
                if l_results[0][1] != l_results[1][1]:
                    c0 = l_results[0][0]
                    c1 = l_results[1][0]
                    for c in [c0, c1]:
                        #print(p)
                        fout.write(f"{p}\t{c}\t{d_best_comb_rcapri[c]}\t"
                               f"{d_best_comb[c]}\t{d_best_comb_DockQ[c]}\t{d_best_comb_irms[c]}\t{d_best_comb_path[c]}\n")


    def analyse_dic_iterations(self, list_protocols=['mixed_ali-delim-delim','unpaired_ali-delim-delim','unpaired_ali-delim-delim'], fileout_rootname=None):
        """
        Method to write a file dumping for every version, for the best combined score model out of 5, and for all the
        protocols listed in the input, the dockQ score, it's CAPRI rank and the AF2 scores
        """

        str_protocol = '-'.join(list_protocols)

        if not fileout_rootname:
            fileout = os.path.abspath(os.path.join(RESULTS_DIR, f"Stats_alliterations_{str_protocol}.out"))
        else:
            fileout = os.path.abspath(os.path.join(RESULTS_DIR, f"{fileout_rootname}_{str_protocol}.out"))


        with open(os.path.join(self.pathout, fileout), "w") as fout:
            fout.write(f"# $1:index_pdb\n# $2:protocol\n# $3:CAPRIrankPeptide\n# $4:CombinedSCORE\n# $5:DockQ\n# $6:IRMS\n# $7:Path2File\n")
            #fout.write(f"# TotalCases: {len(self.dg)}\n")
            for p in self.dg:
                for c in list_protocols:
                    for it in range(n_iterations):
                        v = str(it + 1)
                        for m in range(5): # sample around models
                            try:
                                d_best_model_AF2 = self.dg[p][c][v]['sort_by_CombinedSCORE'][m][1]
                            except KeyError:
                                continue
                            best_comb = d_best_model_AF2['CombinedSCORE']
                            best_comb_rcapri = d_best_model_AF2['CAPRIrankPeptide']
                            best_comb_DockQ = d_best_model_AF2['DockQ']
                            best_comb_irms = d_best_model_AF2['IRMS']
                            best_comb_path = d_best_model_AF2['path']
                            fout.write(f"{p}\t{c}\t{best_comb_rcapri}\t{best_comb}\t{best_comb_DockQ}\t{best_comb_irms}\t{best_comb_path}\n")

    def sample_protocols(self,list_protocols=['mixed_ali-delim-delim',
                                                'mixed_ali-delim-delim',
                                                'unpaired_ali-delim-delim',
                                                'unpaired_ali-delim-delim',
                                                'mixed_ali-delim-100',], repeat=100):

        l_counter_score = []
        for r in range(repeat):
            counter_score = 0
            for p in self.dg:
                l_rescombprotocol = list()
                list_selection = []
                v_num = random.randint(1, 5)
                selection = list_protocols[0] + '_v' + str(v_num)
                list_selection.append(selection)
                for c in list_protocols[1:]:
                    while selection in list_selection:
                        v_num = random.randint(1,5)
                        selection = c + '_v' + str(v_num)
                    list_selection.append(selection)
                    l_rescombprotocol.append(self.dg[p][c][str(v_num)]['sort_by_CombinedSCORE'][0])
                l_rescombprotocol.sort(key=lambda x: x[0], reverse=True)
                best_comb_rcapri = l_rescombprotocol[0][1]['CAPRIrankPeptide']
                if best_comb_rcapri in ['High', 'Medium', 'Acceptable']:
                    counter_score += 1
            #print(f"Retrieved {counter_score} out of {len(self.dg)}")
            l_counter_score.append(counter_score)
        sum_lst = sum(l_counter_score)
        lst_avg = sum_lst / len(l_counter_score)
        print(f"On average, retrieved {lst_avg} cases out of {len(self.dg)}")


if __name__ == '__main__':

    AR = AnalyseResults(fglobaloutput_path)

