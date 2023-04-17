#!/usr/bin/env python


import argparse, os
import shutil
import json
import copy


def add_input_seq(l_inputseq, lengths, code_found_notfound):
    """

    Args:
        l_inputseq:
        lengths:
        code_found_notfound: string of the form 0_1_1 for instance

    Returns: the string to add in the fasta file with the header and the sequence of the input segments

    """
    str_input_seq = ""
    l_code = code_found_notfound.split('_')
    l_header = []
    l_seq = []
    count_1 = l_code.count('1')
    if count_1 == 1:
        single_sequence = True
    else:
        single_sequence = False
    for ii, code in enumerate(l_code):
        if code == '1':
            l_header.append(f"10{ii + 1}")
            header_input = list(l_inputseq[ii].keys())[0]
            l_seq.append(l_inputseq[ii][header_input])
        elif not single_sequence:
            # If it is above 2 sequences, means it is an interaction between two members so we write a header with several items
            l_header.append(f"NotFound-10{ii + 1}") # added to have the same number of fields in all headers
            l_seq.append("-" * lengths[ii])
        else:
            # It is a single sequence. The header should only contain one header otherwise colabfold will complain
            l_seq.append("-" * lengths[ii])
    str_input_seq += '>{}\n'.format("\t".join(l_header))
    str_input_seq += ''.join(l_seq)
    str_input_seq += '\n'

    return str_input_seq


def check_not_found(l_headers):
    """

    Args:
        l_headers: From the headers, identifies which partner is present or not present so as to add the required sequence headers

    Returns: a string of 0 or 1 of the form '0_1_1' to notify which one is present/absent
            None if only NotFound sequences
            'z' if all sequences were found => no need to add a specific input sequence.

    """
    l_code = []
    for header in l_headers:
        if header[:8] == "NotFound":
            l_code.append('0')
        else:
            l_code.append('1')
    if '0' not in l_code:  # only 1
        str_l_code = 'z' # used to sort it in the first place
    elif '1' not in l_code:
        str_l_code = None
    else:
        str_l_code = '_'.join(l_code)
    return str_l_code


def catenate(l_msa, fileout, l_inputseq=False, lengths=False, stoichiometry = False, add_top_general_header=True, add_top_individual_header=False, unpaired_mode=False):
    """join several MSA line by line. if lengths an stoichiometry: header are sep with tab. A first line is added: lengths of sequences sep by a ',' and stoichiometry"""
    #print("l_inputseq", l_inputseq)
    if lengths:
        separator = '\t'
    else:
        separator = '|'
    
    with open(fileout, 'w') as fout:
        #for ii, msa in enumerate(l_msa):
        if l_inputseq and add_top_general_header:
            # Write the input seqs at the top of the alignments for combined inputs
            if lengths:
                # lengths variable controls the format for colab input
                supplementary_line = "#"+",".join(str(x) for x in lengths) + '\t' + ",".join(str(x) for x in stoichiometry) + "\n"
                fout.write(supplementary_line)
                header = ">" + separator.join(f'10{x}' for x in range(1, len(stoichiometry)+1))+'\n'
            else:
                # Format not for colab std. Might be to be adjusted depending the format
                header = ">" + separator.join(list(x.keys())[0] for x in l_inputseq)+'\n'
            fout.write(header)
            seq = "".join(list(x.values())[0] for x in l_inputseq)+'\n'
            fout.write(seq)

        l_headers_to_sort = []
        for jj, head in enumerate(l_msa[0]):
            # We check the status of the not found sequences and add a specific header for the input sequence if required
            l_headers = [l_msa[x][head][0] for x in range(len(l_msa))]
            l_sequences = [l_msa[x][head][1] for x in range(len(l_msa))]
            code_found_notfound = check_not_found(l_headers)
            if code_found_notfound:
                # 'z'=1_1_1, then 1_1_0, 1_0_1 should come first and so on, respecting otherwise the original sequence order
                l_headers_to_sort.append([code_found_notfound, len(l_msa[0]) - jj, l_headers, l_sequences])

        # The alignment will be organized as a series of blocks of paired, partially paired, unpaired aligned sequences
        l_headers_to_sort.sort(reverse=True)
        
        previous_code_found_notfound = ""

        # In unpaired mode we don't want partially or totally paired sequence
        if not unpaired_mode:
            for jj, [ code_found_notfound, _, l_headers, l_sequences ] in enumerate(l_headers_to_sort):

                if code_found_notfound != 'z' and code_found_notfound != previous_code_found_notfound:
                    # We are entering a new block in the alignment
                    # Hence we write the two input lines such as  >101\tNotFound-102\t103, or >102\n... etc
                    # 1_1_1 case is added in the input at first. We don't need to repeat
                    str_input_seq = add_input_seq(l_inputseq, lengths, code_found_notfound)
                    fout.write(str_input_seq)
                    previous_code_found_notfound = code_found_notfound

                # for every header of the list. All ali have the same number of headers because the singleton are already prepared with gaps
                header = ">"+separator.join(l_headers)+'\n'
    
                try:
                    seq = "".join(l_sequences)+'\n'
                    if len(set(seq)) == 2: # If the sequence is only made of gaps (due to delimitations for instance, we don't consider it)
                        keep_record = False
                    else:
                        keep_record = True
                except:
                    print("Problem in recovering the sequences for ", header)
                if keep_record:
                    fout.write(header)
                    fout.write(seq)
        # In unpaired_mode we add all the sequences as single inputs shifted by gaps
        else:
            for jj, [ code_found_notfound, _, l_headers, l_sequences ] in enumerate(l_headers_to_sort):
                if code_found_notfound.count('1') == 1 and code_found_notfound != previous_code_found_notfound: #strict unpaired cases
                    # We are entering a new block in the alignment for single sequences
                    # Hence we write the two input lines such as  >101\n
                    str_input_seq = add_input_seq(l_inputseq, lengths, code_found_notfound)
                    fout.write(str_input_seq)
                    previous_code_found_notfound = code_found_notfound 
                # for every header of the list. All ali have the same number of headers because the singleton are already prepared with gaps
                header = ">"+separator.join(l_headers)+'\n'
                try:
                    seq = "".join(l_sequences)+'\n'
                    if len(set(seq)) == 2: # If the sequence is only made of gaps (due to delimitations for instance, we don't consider it)
                        keep_record = False
                    else:
                        keep_record = True
                except:
                    print("Problem in recovering the sequences for ", header)
                if keep_record:
                    fout.write(header)
                    fout.write(seq)
        fout.close()

def appendFile(file1,file2):
    """catenate 2 file. Equivalent to cat in bash"""

    with open(file1, "a") as f1:
        with open(file2, "r") as f2:
            for line in f2.readlines():
                f1.write(line)

def Writeallunpaired(l_msa, fileout, l_inputseq, lenghts, stoichiometry, delimitation):
    separator = '\t'
    with open(fileout, 'w') as fout:
        # for ii, msa in enumerate(l_msa):
        if l_inputseq:
            # Write the input seqs at the top of the alignments for combined inputs
            # format for colab input
            supplementary_line = "#" + ",".join(str(x) for x in lenghts) + '\t' + ",".join(
                str(x) for x in stoichiometry) + "\n"
            fout.write(supplementary_line)
            header = ">" + separator.join(f'10{x}' for x in range(1, len(stoichiometry) + 1)) + '\n'
            fout.write(header)
            seq = "".join(list(x.values())[0] for x in l_inputseq) + '\n'
            fout.write(seq)
        nb_msa=len(l_msa)
        for ii, msa in enumerate(l_msa):
            fout.write(f">10{ii}\n")
            nb_gap_before = sum(lenghts[0:ii])
            nb_gap_after = sum(lenghts[ii+1:])
            sequence = '-'*nb_gap_before + l_inputseq[ii][list(l_inputseq[ii].keys())[0]] + '-'*nb_gap_after # only one seq in l_inputseq[ii][list(l_inputseq[ii].keys())[0]]
            fout.write(f"{sequence} \n")
            print("msa", msa.keys())
            start, stop = delimitation[ii].split('_')
            if start == '':
                start = 1
            start = int(start) - 1
            if start < 0:
                start = 0
            if stop == '':
                stop = None
            else:
                stop = int(stop)
            if 'paired' in msa:
                for num_sequence in msa['paired']:
                    header = ">" + msa['paired'][num_sequence][0]+'\n'
                    fout.write(header)
                    nb_gap_before = sum(lenghts[0:ii])
                    nb_gap_after = sum(lenghts[ii + 1:])
                    sequence = '-' * nb_gap_before + msa['paired'][num_sequence][1][start:stop] + '-' * nb_gap_after
                    fout.write(f"{sequence} \n")
            if 'unpaired' in msa:
                for num_sequence in msa['unpaired']:
                    header = ">" + msa['unpaired'][num_sequence][0]+'\n'
                    if "NotFound" not in header:
                        fout.write(header)
                        nb_gap_before = sum(lenghts[0:ii])
                        nb_gap_after = sum(lenghts[ii + 1:])
                        sequence = '-' * nb_gap_before + msa['unpaired'][num_sequence][1][start:stop] + '-' * nb_gap_after
                        fout.write(f"{sequence} \n")


# def Writepairedgap(l_msa, l_msa_header, fileout, l_inputseq=False):
#     """Write an unpaired alignement with sequence present in the paired msa """
#     for ii,msa in enumerate(l_msa):
#         with open(fileout, 'w') as fout:
#             header=">"+"|".join(list(x.keys())[0] for x in l_inputseq)+'\n'
#             fout.write(header)
#             seq="".join(list(x.values())[0] for x in l_inputseq)+'\n'
#             fout.write(seq)
#
#             for x in range(len(l_msa_header)):
#                 for jj, head in enumerate(l_msa[0]):
#                     header=">"+"|".join(' ' for y in range(x))+l_msa_header[x][jj]+"|".join(' ' for x in range(len(l_msa_header)-x))+'\n'
#                     fout.write(header)
#                     seq=''.join("-"*len(list(l_inputseq[y].values())[0]) for y in range(x))+ l_msa[x][l_msa_header[x][jj]]+''.join("-"*len(list(l_inputseq[y+x].values())[0]) for y in range(1, len(l_msa_header)-x))+'\n'
#                     fout.write(seq)



def TruncMSA(d_Truncated, start, stop, is_input_seq=False):
    """Trunc each seq of a msa between start and stop index """
    l_headerTruncated = list()
    if is_input_seq:
        for k in d_Truncated:
            d_Truncated[k] = d_Truncated[k][start:stop]
            l_headerTruncated.append(k)
    else:
        for index, el in enumerate(d_Truncated):
            d_Truncated[el][1] = d_Truncated[el][1][start:stop]
            l_headerTruncated.append(d_Truncated[el][0])
    return d_Truncated, l_headerTruncated


def WriteHeaderAlphaFold(list_msa_trunc, fileout):
    """Write concatenation of first sequences joined by ':' """
    with open(fileout, "w") as f:
        f.write(":".join(list(msa.values())[0] for msa in list_msa_trunc))


def WriteAllInputsAlphaFoldColab(list_msa_trunc, lengths, fileout):
    """Write input seqs in an unpaired way """
    with open(fileout, "w") as f:
        for i, msa in enumerate(list_msa_trunc):
            f.write(f'>10{i+1}\n')
            before = "".join("-"*lengths[x] for x in range(i))
            after = "".join("-"*lengths[x] for x in range(i+1, len(lengths)))
            seq = list(msa.values())[0]
            f.write(before+seq+after+'\n')


def WriteSpecificInputsAlphaFoldColab(list_msa_trunc, index, lengths, fileout, writing_mode="a"):
    """Write input seqs in an unpaired way """
    with open(fileout, writing_mode) as f:
        msa = list_msa_trunc[index]
        #for i, msa in enumerate(list_msa_trunc):
        f.write(f'>10{index+1}\n')
        before = "".join("-"*lengths[x] for x in range(index))
        after = "".join("-"*lengths[x] for x in range(index+1, len(lengths)))
        seq = list(msa.values())[0]
        f.write(before+seq+after+'\n')


def delimiter(list_msa_input, delimitation, is_input_seq=False):
    """For a list of msa: Check validity of start and stop index 
    Call TruncMSA and return list of truncated msa"""
    list_msa_trunc=list()
    l_msa_headerTruncated=list()
    for ii, el in enumerate(delimitation):
        start, stop = el.split('_')
        if start == '':
            start = 1
        start = int(start)-1
        if start < 0:
            start = 0
        if stop == '':
            stop = None
        else:
            stop = int(stop)
        msa_trunc, l_headerTruncated = TruncMSA(list_msa_input[ii], start, stop, is_input_seq=is_input_seq)
        list_msa_trunc.append(msa_trunc)
        l_msa_headerTruncated.append(l_headerTruncated)
    return list_msa_trunc, l_msa_headerTruncated


def PrepareJson(output_json,l_msa, key):
    """complete data to be saved in json format"""
    
    for ii, msa in enumerate(l_msa):
        output_json[ii][key] = l_msa[ii]
    return output_json



def main(FLAGS):

    if os.path.isdir('output') and FLAGS.new_directory:
        shutil.rmtree('output')
    if not os.path.isdir('output'):
        os.mkdir('output')


    list_msa_input = (FLAGS.msa_order).split(',') #0,0,1,2
    delimitation = (FLAGS.delimitation).split(',')
    mode = FLAGS.output.split(',')
    print('mode', mode)
    specificname = FLAGS.name
    jsonnameoutput = FLAGS.jsonoutput
    
    if specificname != '':
        specificname = '_' + specificname
    
    filenameIN = FLAGS.jsonname+'.json'

    with open(filenameIN, 'r') as fin:
        data0 = json.load(fin) # [{'inputseq':{header:seq},'paired':{index:[headers,seqs]},'unpaired':{index: [headers, seqs]} }{...}{...}]
    
    data = list()
    for ii in list_msa_input:
        data.append(copy.deepcopy(data0[int(ii)]))
    
    # Treatment of the input sequence
    data_inputseq = [ copy.deepcopy(x['inputseq']) for x in data ] # list of [{header:seq},{header:seq}... ]
    list_msa_trunc_inputseq, _ = delimiter(data_inputseq, delimitation, is_input_seq=True) # list of truncated [{header:seq},{header:seq}... ]
    output_json = [dict() for x in range(len(list_msa_input))]
    output_json = PrepareJson(output_json, list_msa_trunc_inputseq, 'inputseq')
    
    if FLAGS.colab:
        # Format for the colab server
        lengths = [len(list(x.values())[0]) for x in list_msa_trunc_inputseq]
        if FLAGS.stoichiometry == '':
            stoichiometry = [1 for x in list_msa_trunc_inputseq]
        else:
            stoichiometry = (FLAGS.stoichiometry).split(',')
    else:
        stoichiometry = False
        lengths = False
    
    if 'paired' in mode:
        fileout = f"output/MSAforAlphafold_paired{specificname}.a3m"
        list_msa_input_paired = [ copy.deepcopy(x['paired']) for x in data ]

        # Apply the delimitations on the paired alignments:
        list_msa_trunc, l_msa_headerTruncated = delimiter(list_msa_input_paired, delimitation)
        catenate(list_msa_trunc, fileout, l_inputseq=list_msa_trunc_inputseq, lengths=lengths, stoichiometry=stoichiometry, unpaired_mode=False)
        output_json = PrepareJson(output_json, list_msa_trunc, 'paired')

        if FLAGS.colab:
            # Writing the individual input sequences in unpaired format.
            # As such ok for paired only alignment but not for unpaired ones
            #   since the input should a head of every individual alignment
            fileoutColab = f"output/Alphafold_input_Sequence_COLAB{specificname}"
            WriteAllInputsAlphaFoldColab(list_msa_trunc_inputseq, lengths, fileoutColab)
            # The input individual seq only is written at the end to the alignment
            appendFile(fileout, fileoutColab)
            if os.path.isfile(f"output/Alphafold_input_Sequence_COLAB{specificname}"):
                os.remove(f"output/Alphafold_input_Sequence_COLAB{specificname}")

    if 'mixed' in mode:
        fileout = f"output/MSAforAlphafold_mixed{specificname}.a3m"
        list_msa_input_paired = [copy.deepcopy(x['paired']) for x in data]

        # Apply the delimitations on the paired alignments:
        list_msa_trunc, l_msa_headerTruncated = delimiter(list_msa_input_paired, delimitation)
        catenate(list_msa_trunc, fileout, l_inputseq=list_msa_trunc_inputseq, lengths=lengths,
                 stoichiometry=stoichiometry,unpaired_mode=False)
        output_json = PrepareJson(output_json, list_msa_trunc, 'mixed')

        # Check if this is still functional => not still functionnal
        #if FLAGS.gappaired:
        #    fileout = f"output/MSAforAlphafold_pairedwithgap{specificname}.a3m"
        #    Writepairedgap(list_msa_trunc, l_msa_headerTruncated, fileout, l_inputseq=list_msa_trunc_inputseq)
        """
        if FLAGS.colab:
            # Writing the individual input sequences in unpaired format at the end of the alignment.
            # If we want partially paired alignment to be integrated with the forcepaired.sif singularity image we need this
            # Not sure it is compatible with the previous unchanged colab --> Needs to be checked
            fileoutColab = f"output/Alphafold_input_Sequence_COLAB{specificname}"
            WriteAllInputsAlphaFoldColab(list_msa_trunc_inputseq, lengths, fileoutColab)
            # The input individual seq only is written at the end to the alignment
            appendFile(fileout, fileoutColab)
            if os.path.isfile(f"output/Alphafold_input_Sequence_COLAB{specificname}"):
                os.remove(f"output/Alphafold_input_Sequence_COLAB{specificname}")
        """
        #file without the first sequence
        fileout = f"MSAforAlphafold_unpaired_noinputseq{specificname}.a3m"
        list_msa_input_unpaired = [ copy.deepcopy(x['unpaired']) for x in data]

        list_msa_trunc, l_msa_headerTruncated = delimiter(list_msa_input_unpaired, delimitation)
        catenate(list_msa_trunc, fileout, l_inputseq=list_msa_trunc_inputseq, lengths=lengths,
                 stoichiometry=stoichiometry, add_top_general_header=False, add_top_individual_header=True, unpaired_mode=False) #Nipponia & Eubucco ne sont plus semblables

        output_json = PrepareJson(output_json, list_msa_trunc, 'unpaired')
        #fileout = f"output/MSAforAlphafold_mixed{specificname}.a3m"
        
        with open(f"output/MSAforAlphafold_mixed{specificname}.a3m", 'a') as fout:
            #with open(f"output/MSAforAlphafold_paired{specificname}.a3m", 'r') as fin:
            #    fout.write("".join([line for line in fin.readlines()]))

            with open(f"MSAforAlphafold_unpaired_noinputseq{specificname}.a3m", 'r') as fin:
                fout.write("".join([line for line in fin.readlines()]))
        for index in range(len(data_inputseq)):
            WriteSpecificInputsAlphaFoldColab(list_msa_trunc_inputseq, index, lengths,
                                              f"output/MSAforAlphafold_mixed{specificname}.a3m", writing_mode="a")
    
    if 'unpaired' in mode:
        fileout = f"output/MSAforAlphafold_unpaired{specificname}.a3m"
        
        list_msa_input_unpaired = [ copy.deepcopy(x['unpaired']) for x in data]

        #print('list_msa_input_unpaired', list_msa_input_unpaired)
        list_msa_trunc, l_msa_headerTruncated = delimiter(list_msa_input_unpaired, delimitation)
        catenate(list_msa_trunc, fileout, l_inputseq=list_msa_trunc_inputseq, lengths=lengths,
                 stoichiometry=stoichiometry, add_top_general_header=True, add_top_individual_header=False, unpaired_mode=True)
        
        output_json = PrepareJson(output_json, list_msa_trunc, 'unpaired')

    if 'all_unpaired' in mode:
        list_msa_input_paired = [copy.deepcopy(x['paired']) for x in data]
        list_msa_input_unpaired = [ copy.deepcopy(x['unpaired']) for x in data]
        list_msa_input = list_msa_input_paired + list_msa_input_unpaired
        list_msa_trunc, l_msa_headerTruncated = delimiter(list_msa_input, delimitation)
        fileout = f"output/MSAforAlphafold_all_unpaired{specificname}.a3m"
        Writeallunpaired(l_msa = data, fileout = fileout, l_inputseq = list_msa_trunc_inputseq, lenghts = lengths, stoichiometry = stoichiometry, delimitation= delimitation)

    if 'single_pep' in mode:
        list_msa_input_paired = [copy.deepcopy(x['paired']) for x in data]
        list_msa_input_unpaired = [copy.deepcopy(x['unpaired']) for x in data]
        list_msa_input = list_msa_input_paired + list_msa_input_unpaired
        list_msa_trunc, l_msa_headerTruncated = delimiter(list_msa_input, delimitation)
        fileout = f"output/MSAforAlphafold_single_pep{specificname}.a3m"
        Writeallunpaired(l_msa = data, fileout = fileout, l_inputseq = list_msa_trunc_inputseq, lenghts = lengths, stoichiometry = stoichiometry, delimitation= delimitation)

    fileout_1stSeq = f"output/Alphafold_input_Sequence{specificname}"

    WriteHeaderAlphaFold(list_msa_trunc_inputseq, fileout_1stSeq)

    with open(f"output/{jsonnameoutput}.json", 'w') as fj:
        json.dump(output_json, fj, indent=4)
    
    if os.path.isfile(f"MSAforAlphafold_unpaired_noinputseq{specificname}.a3m"):
        os.remove(f"MSAforAlphafold_unpaired_noinputseq{specificname}.a3m")



if __name__=="__main__":
    usage = "GetAlphafoldInput.py -msa 0,0,0,1,2,etc -d 10_50,70_100,120_,_60,30_100 -o mixed -n name_suffix"
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument(
        '-msa',
        '--msa_order',
        type=str,
        required=True,
        help="indexes of msas, sep by a comma ex 0,0,1,2"
    )
    parser.add_argument(
        '-json',
        '--jsonname',
        type=str,
        required=True,
        help="json output name of CommonAndGapspecies script without the .json extension"
    )
    parser.add_argument(
        '-d',
        '--delimitation',
        type=str,
        required=True,
        help="Must have same len as msa. for each msa, start and stop index separated by _.  separation by a , between each msa : ex if 3msa: 1_150, _60, 30_. First index is 1"
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        choices=['paired', 'unpaired', 'mixed', 'all_unpaired', 'single_pep'],
        default='paired,unpaired,mixed',
        help="paired or unpaired or mixed"
    )
    parser.add_argument(
        '-n',
        '--name',
        type=str,
        default='',
        help="if -n <specicificname> specified, a suffix will be added to output file names"
    )
    parser.add_argument(
        '-jo',
        '--jsonoutput',
        type=str,
        required=True,
        help="output name of json file without the json extention"
    )
    parser.add_argument(
        '-colab',
        '--colab',
        action='store_true',
        help="if --colab, files are generated with specific header and format to run on colab scripts. if colab please complete -stoichiometry option"
    )
    parser.add_argument(
        '-st',
        '--stoichiometry',
        type=str,
        default='',
        help="stoichiometry of the complex sep by a ',': ex 1,2,1,1. Must have same lenght than msa_order and delimitation. By default: a stoichiometry of 1 for each chain is assigned"
    )
    parser.add_argument(
        '-newdir',
        '--new_directory',
        action='store_true',
        help="Delete the previous output directory and replace it by a new one."
    )
    FLAGS = parser.parse_args()
    main(FLAGS)
