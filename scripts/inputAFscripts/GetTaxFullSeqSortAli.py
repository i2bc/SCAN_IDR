#!/usr/bin/env python3

import os
import re
import time
import pandas as pd
import sys
import argparse




def sort_gid_from_qid_cov_tsv(input_tsv, keep_idcov = False):
    """
    Method which computes a score for every sequence of an alignment and sort them with respect to qid and cov
    This is required because the evalue of mmseqs-profile search are not reliable
    Args:
        input_tsv:

    Returns:

    """
    df = pd.read_csv(input_tsv, sep='\t')
    d_similarity = dict()
    # We assume that every non-aligned position in the query is equivalent to a nearly random match that we estimate at about 15%
    # which is about threshold limit for not significant percentage identity of none related sequences
    # This threshold value may be changed
    unmatched_pos_eq_id = 15.
    df['score2sort'] = (df['pident'] * df['qcov'] + unmatched_pos_eq_id * (1. - df['qcov']))
    dfsort = df.sort_values(["score2sort", "evalue"], ascending=(False, True))
    df['score'] = df['score2sort'].apply(lambda x: round(x, 2))
    list_sorted_gid = dfsort['target'].tolist()
    d_target_to_fullseq = df.set_index('target').to_dict('dict')['tseq']
    d_target_to_taxname = df.set_index('target').to_dict('dict')['taxname']
    if keep_idcov:
        # we keep only the first occurrence in case of duplicated or more entries
        #print(df[['target', 'pident']].drop_duplicates(subset='target').itertuples().head(3))
        d_similarity['id'] = {row.target: row.pident for row in df[['target', 'pident']].drop_duplicates(subset='target').itertuples()}

        d_similarity['cov'] = {row.target: row.qcov for row in df[['target', 'qcov']].drop_duplicates(subset='target').itertuples()}
        d_similarity['score'] = {row.target: row.score for row in df[['target', 'score']].drop_duplicates(subset='target').itertuples()}

    return list_sorted_gid, d_target_to_fullseq, d_target_to_taxname, d_similarity


def get_taxid_from_header(header):
    pat = re.compile("TaxID=(\d+)")
    hh = pat.search(header)
    if hh:
        taxid = hh.group(1)
    else:
        taxid = None
    return taxid


def change_definitions_in_ali(d_id2tax, input_ali, output_ali, sort_order=[], keep_first_sequence=True, d_similarity_info=dict()):
    """

    Args:


    """

    seq = 0
    d_input_fasta_file = {}
    d_input_fasta_file['order'] = []
    l_first_sequence = ["", ""] # header and sequence of the first line
    with open(input_ali, 'r') as fid_in:
        with open(output_ali, "w") as fid_out:

            previous_time = time.time()
            for ii, ll in enumerate(fid_in):
                if ll[0] == ">":
                    line_index_header = ii
                    if keep_first_sequence and seq == 0:
                        fid_out.write(ll)
                        l_first_sequence[0] = ll
                        seq += 1
                        continue
                    gid = ll.split()[0][1:]
                    if gid not in d_input_fasta_file:
                        # Several matches with same gids but different alignments can have been detected
                        # We only keep the one which came first
                        d_input_fasta_file[gid] = {}
                        d_input_fasta_file['order'].append(gid) # dict are now ordered in py3.8 but just in case
                        DO_keep_this_sequence = True
                    else:
                        DO_keep_this_sequence = False
                    seq += 1
                    if DO_keep_this_sequence:
                        taxname = d_id2tax[gid]

                        if len(d_similarity_info) == 0:
                            d_input_fasta_file[gid]['header'] = f"{gid} n=1 Tax={taxname}"
                        else:
                            id = d_similarity_info['id'][gid]
                            cov = d_similarity_info['cov'][gid]
                            score = d_similarity_info['score'][gid]
                            d_input_fasta_file[gid]['header'] = f"{gid} n=1 Tax={taxname} Id={id} Cov={cov} Score={score}"
                        d_input_fasta_file[gid]['sequence'] = ""
                    else:
                        pass
                    """
                    except KeyError:
                        #fid_out.write(">{}\n".format(gid))
                        d_input_fasta_file[gid]['header'] = ""
                        d_input_fasta_file[gid]['sequence'] = ""
                    """
                    #if seq % 10000000 == 0:
                    #    next_time = time.time()
                    #    timelapse = time.strftime("%H:%M:%S", time.gmtime(next_time - previous_time))
                    #    print("Processed {} entries from the original fasta file in {}".format(seq, timelapse))
                    #    previous_time = next_time
                elif len(ll) > 0 :
                    if keep_first_sequence and seq == 1:
                        # if first sequence is on multilines, we will pass several times here to collect all the lines
                        fid_out.write(ll)
                        l_first_sequence[1] += ll # we keep it as reference for full length file
                        continue
                    # add this constraint to prevent that trailing lines in the end of the file crash the last sequence
                    if DO_keep_this_sequence: # This is the first instance we see this entry
                        d_input_fasta_file[gid]['sequence'] += ll
            if len(sort_order) > 0:
                set_added_git = set()
                for gid in sort_order:
                    if gid in set_added_git:
                        continue
                    set_added_git.add(gid)
                    if not gid in d_input_fasta_file:
                        print("Problem could not recover one entry. Check bug in inputs format since it should not happen")
                        sys.exit()
                    header = d_input_fasta_file[gid]['header']
                    sequence = d_input_fasta_file[gid]['sequence']
                    fid_out.write(">{}\n{}".format(header, sequence))
            else:
                for gid in d_input_fasta_file:
                    header = d_input_fasta_file[gid]['header']
                    sequence = d_input_fasta_file[gid]['sequence']
                    fid_out.write(">{}\n{}".format(header, sequence))

    fid_out.close()

    return d_input_fasta_file, l_first_sequence


def save_full_length_file(d_gid2seq, d_gid2header, full_length_file, l_first_sequence, sort_order=None):
    """
    Args:
        d_gid2seq:
        d_gid2header:
        full_length_file: Name of the output file
        l_first_sequence:
        sort_order:
    Returns:

    """
    with open(full_length_file, "w") as fid_out:
        fid_out.write("{}{}".format(l_first_sequence[0], l_first_sequence[1]))
        if sort_order:
            set_added_git = set()
            for gid in sort_order:
                if gid in set_added_git:
                    continue
                set_added_git.add(gid)
                header = d_gid2header[gid]['header']
                seq = d_gid2seq[gid]
                fid_out.write(">{}\n{}\n".format(header, seq))
        else:
            for gid in d_gid2header['order']:
                header = d_gid2header[gid]['header']
                seq = d_gid2seq[gid]
                fid_out.write(">{}\n{}\n".format(header, seq))
    fid_out.close()


def main(input, auto_output=True, output='', full_length=False, file_qidcov=''):

    start_time = time.time()

    l_ifiles = input.split(',')
    print("Files which will be processed :")
    for ff in l_ifiles:
        print("\t",ff)

    if not auto_output:
        if len(output) == 0:
            print("Please provide the names of the outputfiles separated by commas")
            sys.exit()
        l_ofiles = output.split(',')

    if len(file_qidcov) > 0:
        l_qfiles = file_qidcov.split(',')
        d_in2qidcov = dict(zip(l_ifiles, l_qfiles))
        d_qfiles_to_sorted_gid_list = {}
        d_qfiles_to_id2fullseq = {}
        d_qfiles_to_id2taxname = {}
        d_qfiles_to_similarity = {}
        for qfile in l_qfiles:
            if qfile not in d_qfiles_to_sorted_gid_list:
                d_qfiles_to_sorted_gid_list[qfile], \
                d_qfiles_to_id2fullseq[qfile], \
                d_qfiles_to_id2taxname[qfile], \
                d_qfiles_to_similarity[qfile] = sort_gid_from_qid_cov_tsv(qfile, keep_idcov=True)

    # Modify every file in the input list
    for ii, ifile in enumerate(l_ifiles):
        print("Running", ifile)
        ipath, iext = os.path.splitext(ifile)
        if output:
            ofile = l_ofiles[ii]
        else:
            ofile = "".join([ipath, '_withdef', iext])
        if full_length:
            fl_file = "".join([ipath, '_FL.fasta'])

        # Important key corresponds to the name of the tsv file
        qfile = d_in2qidcov[ifile]

        print("Starting reading the fasta/a3m file")
        if len(file_qidcov) > 0:
            sorted_list_of_gid = d_qfiles_to_sorted_gid_list[qfile]
        else:
            sorted_list_of_gid = []
            # The lists recovered are useful for full length dumping

        d_gid2header, l_first_sequence = change_definitions_in_ali(d_qfiles_to_id2taxname[qfile], ifile, ofile,
                                                                sort_order=sorted_list_of_gid,
                                                                keep_first_sequence=True,
                                                                d_similarity_info=d_qfiles_to_similarity[qfile]
                                                                )

        if full_length:
            save_full_length_file(d_qfiles_to_id2fullseq[qfile], d_gid2header, fl_file, l_first_sequence, sort_order=sorted_list_of_gid)

    step2_time = time.time()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
    '-i',
    '--input',
    type=str,
    required=True,
    help='List of comma-separated fasta or a3m file.s Example ESA1_YEAST.fasta,EPL1_YEAST.a3m, ...'
    )
    parser.add_argument(
    '-o',
    '--output',
    type=str,
    required=False,
    default=False,
    help='List of comma-separated names of outputfile fasta or a3m files. '
         'Default is same as input with "_withdef" added in namefile'
    )
    parser.add_argument(
    '-l',
    '--full_length',
    required=False,
    default=False,
    action='store_true',
    help='If added, this option will trigger the recovery of full length sequences: '
    )
    parser.add_argument(
    '-q',
    '--qidcovfile',
    type=str,
    required=True,
    help='List of comma-separated tsv files containing qid and cov data to sort sequences for every sequence in the input. Example ESA1_YEAST.tsv,EPL1_YEAST.tsv, ...'
    )
    options = parser.parse_args()

    main(options.input, output=options.output, full_length=options.full_length, file_qidcov=options.qidcovfile)
