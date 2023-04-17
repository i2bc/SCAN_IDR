#!/usr/bin/env python
import copy
import os,re
import argparse
import json


def extract_species(header):
    """Identify a species from a header
    """
    pat = re.compile("OS=(\S+\s\S*)")
    pat2 = re.compile("Tax=(\S+\s\S*)")
    
    try:
        species = pat.findall(header)[0]
    except IndexError:
        try:
            species = pat2.findall(header)[0]
        except IndexError:
            species = header.split()[0].split("_")[-1].split("/")[0]
    
    return species


def clean_strange_characters(msa):
    """Replace no conventional amino acids or insertion by -"""
    clean_msa=dict()
    set_authorized_letters = set(["-","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","\n"])
    for header,seqs in msa.items():
        clean_msa[header]=""
        for letter in seqs:
            if letter not in set_authorized_letters:
                clean_msa[header] += "-"
            else:
                clean_msa[header] += letter
    return clean_msa


def read_fasta(infile):
    """
    Parse msa file to dict {header:seq}
    Keep the first sequence data in a separate dict
    Delete the first sequence of the dictionary
    Args:
        infile: name of the input fasta or a3m file
    Returns: arg1 -> dictionary  all but the first sequence {header:seq}
             arg2 -> dictionary  first sequence  {header:seq}
    """
    results = {}
    first_sequence = 0
    ordered_list_headers = []
    for l in open(infile):
        if l.startswith(">"):
            name = l.strip()[1:]
            results[name] = ""
            if first_sequence == 0:
                input1 = name
            else:
                ordered_list_headers.append(name)
            first_sequence += 1
        else:
            results[name] += l.strip()
    
    results = clean_strange_characters(results)
    results['order'] = copy.deepcopy(ordered_list_headers)
    msa_input1 = {input1:results[input1]}
    del(results[input1])
    return results, msa_input1


def GetPositionMSA(seqtrunc):
    """Make a dic {position without gap:position in the gapped sequence}"""
    d_eq=dict()
    count=0#initialize count in seq without gap
    for ii, letter in enumerate(seqtrunc): #ii=count in seq with gap
        if letter!='-':
            d_eq[count]=ii
            count+=1
    return d_eq


def CutSequence(seq, list_nongappedindex):
    """write a new sequence that keep only index in list_nongappedindex of the input sequence """
    newseq=""
    for el in list_nongappedindex:
        if len(seq)>el:
            newseq+=seq[el]
        else:
            newseq+='-'
    return newseq


def NonGappedIndex(seq):
    """ return list of index of a sequence that aren't gap and the ungapped sequence"""
    list_nongappedindex=list()
    newseq=""
    for ii, letter in enumerate(seq):
        if letter !='-':
            newseq += letter
            list_nongappedindex.append(ii)
    return list_nongappedindex, newseq


def AnalyseGapsIn1stSeq(list_1st_sequence):
    """
    Keep all non gapped index of the 1st sequence. All following sequence in the MSA will keep only non gapped index
    Args:
        list_1st_sequence: a list of dict() for each alignment of the form {header_1st_seq:seq_1st_seq}

    Returns: d_nongappedindex -> dict {index_of_the_msa:list_of_non_gapped_index_positions}

    """
    d_nongappedindex = dict()
    for ii, msa in enumerate(list_1st_sequence):
        header, seq = list(msa.keys())[0], list(msa.values())[0]
        d_nongappedindex[ii], msa[header] = NonGappedIndex(seq)
    return d_nongappedindex


def Species2Header(msa):
    """
    Give the species for each header {species:header}. One species will be present several times. Save the order of initial msa
    """
    dsp2head = dict()
    dsp2head_order = dict() # for every species contains all the index of entries where the species was found
    dsp2head['order'] = []
    for ii, header in enumerate(msa['order']): # Even if dict are now ordered sicne py3.7 , we prefer using this list as a security
        sp = extract_species(header)
        if sp in dsp2head:
            dsp2head[sp].append(header)
            dsp2head_order[sp].append(ii)
        else:
            dsp2head[sp] = [header]
            dsp2head['order'].append(sp)
            dsp2head_order[sp] = [ii]
    return dsp2head, dsp2head_order


def GetCommonSpecies(dsp2head):
    """output a set of species found in all msa"""
    common_species=set(dsp2head[0].keys())
    for spMSA in dsp2head[1:]:
        common_species = common_species.intersection(spMSA.keys())
    common_species = sorted(common_species, key=lambda x: list(dsp2head[0]).index(x))
    return common_species


def WriteOutput_sep(listMSAinput, msa_common_species, msa_gapped_species,list_1st_sequence,jsonname):
    """
    Write msa of common species for each msa input
    Write msa of not common species with gapped sequence for each msa input
    Args:
        listMSAinput: list of all the msas for every partner
        msa_common_species: an index list of msa {index_of_msa: dict(header:gapped_sequence);
                    every msa_common_species[index_of_msa] contains a key 'order' to give the order of headers
        msa_gapped_species: an index list of msa {index_of_msa: dict(header:gapped_sequence);
                    every msa_gapped_species[index_of_msa] contains a key 'order' to give the order of headers
        list_1st_sequence: list of the form {header_1st_seq:seq_1st_seq}
        jsonname: root name of the output json

    Returns: Save a json file ouput containing all the data for next step

    """

    data = list()

    # Because dict were ordered only recently in py37 and that it is potentially dangerous to use them as ordered objects
    # if we run with py35 for instance, it is better to keep a track of the order of the keys()
    for ii, msa in enumerate(listMSAinput):
        d_msa = dict()
        d_msa['inputseq']={list(list_1st_sequence[ii].keys())[0]:list(list_1st_sequence[ii].values())[0]}
        d_msa['paired'] = dict()
        for species, seq in msa_common_species[ii].items():
            d_msa['paired'][species]=seq
        d_msa['unpaired'] = dict()
        for species, seq in msa_gapped_species[ii].items():
            #if species[0:8]=='NotFound':#header is corresponding species in other msa if species not in this msa 
            #    species = species.split('_')[-1]
            d_msa['unpaired'][species] = seq
        data.append(d_msa)
    with open(jsonname+'.json', 'w') as fj:
        json.dump(data, fj, indent=2, sort_keys=True) # Sure we need to use this sort_keys ?
    fj.close()


def RemoveNotFoundSpecies(d):
    """
    Create a unpaired msa with output entirely gapped sequences/"NotFound" in all msa.
    """
    nbSpecies2check = len(d[0])
    d_new_msa={x:{} for x in d}
    for i_sp in range(nbSpecies2check):
        Found=False
        for i_msa in d:
            if "NotFound" not in d[i_msa][i_sp][0]:
                Found=True
        if Found:
            for i_msa in d:
                d_new_msa[i_msa][i_sp] = d[i_msa][i_sp]
    return d_new_msa


def OrderSpecies(d, ordered_species_list):
    """To keep the sequence order of initial msas for unpaired species
    """
    l_sp = list()
    l_indexes = list()
    nr_l_sp = list() # non redundant version of the list of species
    # k = species
    # v = list of indexes in the ali where the sp was detected
    for k in ordered_species_list: #d.items():
        list_of_seq_index = d[k]
        for rank in list_of_seq_index:
            l_sp.append(k)
            l_indexes.append(rank)
    l_indexes, l_sp = (list(t) for t in zip(*sorted(zip(l_indexes, l_sp))))
    # l_indexes was used to sort species in the order in which they are appearing in the ali
    # Create the same sp list without redundancies
    set_sp = set()
    for sp in l_sp:
        if sp not in set_sp:
            nr_l_sp.append(sp)
            set_sp.add(sp)
    return l_sp, nr_l_sp


def FilterSpeciesAli(ncsname, listMSAinput, jsonname, number_common_species=False, paralogfree=False):
    """

    Args:
        ncsname: Filename of the file containing the number of species
        listMSAinput: list of the input MSA to analyze
        jsonname: name of the json file storing the information
        number_common_species: Flag to specify if we only compute the number of common species without creating the ali
        paralogfree: Do we only keep in unpaired ali seq of species not in the paired

    Returns: Writes the output files

    """
    # Input is such as: ['msa1test.fa', 'msa2test.fa']
    # parse msa files to a list of dict [({header_all_other:seq_all_other},{header_1st_seq:seq_1st_seq},)]
    list1seq_msas = [read_fasta(infile) for infile in listMSAinput]

    list_msas = [x[0] for x in list1seq_msas] # dict of all msas without 1st seq

    # Extracts list_dsp2head_0 as [{species:[headers]}, {species:[order of headers in msa]}]
    list_dsp2head_0 = [Species2Header(msa) for msa in list_msas]
    list_dsp2head = [msa[0] for msa in list_dsp2head_0] # {species:[headers]} for every input msa
    list_dsp2head_order = [msa[1] for msa in list_dsp2head_0] #{species:[order of headers in msa]} for every input msa

    common_species = GetCommonSpecies(list_dsp2head)
    
    # We only compute here the number of common species
    print("number of common species : ", len(common_species)-1) # we remove 1 because we've added the 'order' key in the list
    with open(f"{ncsname}.txt", "w") as f:
        f.write(f"{len(common_species) - 1}") # we remove 1 because we've added the order parameter
    
    if not number_common_species:
        # list1seq_msas is of the form {header_1st_seq:seq_1st_seq}
        list_1st_sequence = [x[1] for x in list1seq_msas] # only 1st seq msas

        # Get a dictionary for each ali with a list of positions without gaps {index_of_the_msa:list_of_non_gapped_index_positions}
        d_nongappedindex = AnalyseGapsIn1stSeq(list_1st_sequence)

        # Initialize {index_of_the_msa:dict()}
        msa_common_species = {x:{} for x in range(len(list_msas))}
        msa_gapped_species = {x:{} for x in range(len(list_msas))}
        #for x in range(len(list_msas)):
        #    msa_common_species[x]['order'] = []
        #    msa_gapped_species[x]['order'] = []
        
        NonRedundantName = 0
        # We enumerate over every input_MSA
        # Hyper dangerous to iterate over a dictionary whose keys are destroyed at the same time
        global_l_immutable = copy.deepcopy(list_dsp2head)
        global_l_immutable_order = copy.deepcopy(list_dsp2head_order)
        seqindex_common_species = 0
        seqindex_unpaired_species = 0
        for ii, dheader2sp in enumerate(global_l_immutable):
            # there is one dheader2sp for every chain
            # 'ii' is the index of the chain
            if len(list_dsp2head_order[ii]) > 0:
                # For a given ali, l_sp lists all the species in the observed order with possible redundancy
                # example [sp1, sp2, sp3, sp1, sp3, ...]
                # Second list is without redundancies
                l_sp, nr_l_sp = OrderSpecies(global_l_immutable_order[ii], dheader2sp['order'])
                if paralogfree:
                    l_sp = copy.deepcopy(nr_l_sp) # Order could be lost with the previous version using list(set()) !!
                for sp in l_sp:
                    if sp in common_species:
                        # for those species the whole set of sequences are present. Typically paired cases
                        for x in range(len(list_msas)):
                            #if x not in list_dsp2head or sp not in list_dsp2head[x]:
                            #    continue
                            # We recover the sequence trimmed from the gaps in the input sequence
                            next_header = list_dsp2head[x][sp][0]
                            seq_no_insert_in_inputseq = CutSequence(list_msas[x][next_header], d_nongappedindex[x])
                            msa_common_species[x][seqindex_common_species] = [next_header, seq_no_insert_in_inputseq]
                            #msa_common_species[x]['order'].append(list_dsp2head[x][sp][0])

                            if paralogfree:
                                del list_dsp2head[x][sp]
                                del list_dsp2head_order[x][sp]
                            else:
                                list_dsp2head[x][sp] = copy.deepcopy(list_dsp2head[x][sp][1:])
                                if list_dsp2head[x][sp] == []:
                                    del list_dsp2head[x][sp]
                                    del list_dsp2head_order[x][sp]
                        seqindex_common_species += 1
                        # always use the first occurrence of each species
                        common_species.remove(sp)
                    else:
                        for x in range(len(list_msas)):
                            if sp in list_dsp2head[x] and len(list_dsp2head[x][sp]) > 0:
                                #msa_gapped_species[x]['order'].append(list_dsp2head[x][sp][0])
                                # If a sequence is present for a given species, it is added instead of gaps
                                #if x not in list_dsp2head or sp not in list_dsp2head[x]: # Not sure this is required, but I prefer to add it by security.
                                #    continue
                                next_header = list_dsp2head[x][sp][0]
                                seq_no_insert_in_inputseq = CutSequence(list_msas[x][next_header], d_nongappedindex[x])
                                msa_gapped_species[x][seqindex_unpaired_species] = [next_header, seq_no_insert_in_inputseq]

                                if paralogfree:
                                    del list_dsp2head[x][sp]
                                    del list_dsp2head_order[x][sp]
                                else:
                                    list_dsp2head[x][sp] = list_dsp2head[x][sp][1:]
                                    if list_dsp2head[x][sp] == []:
                                        del list_dsp2head[x][sp]
                                        del list_dsp2head_order[x][sp]
                            else:
                                # Otherwise gaps over the length of the corresponding sequence are added
                                if sp in list_dsp2head[x] and list_dsp2head[x][sp] == []:
                                    del list_dsp2head[x][sp]
                                    del list_dsp2head_order[x][sp]
                                msa_gapped_species[x][seqindex_unpaired_species] = [f"NotFound{str(NonRedundantName)}_{sp}", "-"*len(d_nongappedindex[x])]
                                NonRedundantName += 1
                        seqindex_unpaired_species += 1
        msa_gapped_species = RemoveNotFoundSpecies(msa_gapped_species)
        WriteOutput_sep(listMSAinput, msa_common_species, msa_gapped_species, list_1st_sequence, jsonname)


def main(ListMSA, fname_json, fname_n_common_sp, only_n_common_sp=False, is_server=False,
         n_MSA='', is_fulllength=False, is_no_paralog=False,):
    if is_server:
        # With this option, filenames are predefined.
        # Number of MSA is required for indexing standard filenames
        if n_MSA=='':
            return "Please enter a value for --numberMSA"
        
        listMSAinput = list()
        if is_fulllength:
            for ii in range(int(n_MSA)):
                listMSAinput.append(f"sequence_{ii}_filtered_FL_mafali.fasta")
        else:
            for ii in range(int(n_MSA)):
                listMSAinput.append(f"sequence_{ii}_filtered.a3m")        
    else:
        listMSAinput = ListMSA.split(',')

    if len(listMSAinput)<2:
        print("Please, provide at least 2 input fasta files separated by a comma")
    
    FilterSpeciesAli(fname_n_common_sp, listMSAinput, fname_json, only_n_common_sp, is_no_paralog)


if __name__=="__main__":
    usage = " CommonAndGapSpecies.py -a MSA_1.fasta,MSA_2.fasta,MSA_3.fasta,etc \n"+\
            "Description: keep only sequences from common species between several MSA. The order of species in the first MSA controls the order in the others"
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument(
        '-msa',
        '--listMSA',
        type=str,
        help="path to fasta files separated by comma"
    )
    parser.add_argument(
        '-ncs',
        '--number_common_species',
        action='store_true',
        help="if -ncs : only print the number of common species, that will be the length of the final coalignment"
    )
    parser.add_argument(
        '-json',
        '--jsonname',
        type=str,
        default="msa",
        help="output name without the .json extension"
    )
    parser.add_argument(
        '-ncsname',
        '--ncsname',
        type=str,
        default="numbercommonspecies",
        help="path to write the number of common species"
    )
    parser.add_argument(
        '-s',
        '--server',
        action='store_true',
        help="if -s : impose the names of output files all being normalized (i.e. extension mafali if option FL). "
             "In server mode, -nb need to be defined. "
             "Also check if you want to modify full length option with -FL. "
    )
    parser.add_argument(
        '-nb',
        '--numberMSA',
        type=str,
        help="Number of msa to analyse. IRuns in server mode only, the script will automatically take the name of each msa path"
    )
    parser.add_argument(
        '-FL',
        '--fulllength',
        action='store_true',
        help="if -FL : Runs in server mode only (-s), the script will run on Full Length MSA"
    )
    parser.add_argument(
        '-pf',
        '--paralogfree',
        action='store_true',
        help="if -pf: the unpaired alignment will contain only once each species and will not contain species present in paired alignment"
    )
    FLAGS = parser.parse_args()

    main(FLAGS.listMSA,
         FLAGS.jsonname,
         FLAGS.ncsname,
         only_n_common_sp=FLAGS.number_common_species,
         is_server=FLAGS.server,
         n_MSA=FLAGS.numberMSA,
         is_fulllength=FLAGS.fulllength,
         is_no_paralog=FLAGS.paralogfree,
         )
