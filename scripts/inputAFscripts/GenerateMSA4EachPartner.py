import copy
import os
import sys
import re
import glob
import socket
import argparse
import configparser

import RunMmseqs
import GetTaxFullSeqSortAli

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(os.path.dirname(script_path))

config = configparser.ConfigParser()
config.read(f'{script_dir}/config.ini')

MMSEQS_BIN = config['PROGRAMS']['MMSEQS_BIN']
MMSEQS_DB = config['PROGRAMS']['MMSEQS_DB']

def main(input_fasta_directory, method, list_files='', name_filter='', run_mmseqs=True, skip_ali_generation=False):
    """
    Pipeline the generation of alignment for a number of N sequences

    """
    cur_dir = os.getcwd()
    print("Starting in directory: ", cur_dir)
    os.chdir(input_fasta_directory)
    print("Moving to directory: ", input_fasta_directory)

    print("Alignment method will be: ", method)

    if len(list_files) == 0:
        cmd = "*{}*.fas*".format(name_filter)
        list_fasta_files = glob.glob("*{}*.fas*".format(name_filter))
    else:
        list_fasta_files = list_files.split(",")

    print("Fasta files that will be used as inputs: {}".format(list_fasta_files))

    path_exe = MMSEQS_BIN

    list_files_to_map = []
    list_files_tsv_qidcov = []
    for ii, fasta_input in enumerate(list_fasta_files):
        print("Running mmseqs on :",fasta_input)
        generated_alis_formats = 'both'

        M = RunMmseqs.MMSEQS(fasta_input, output_file=None, ali_method=method,
                             executable_path=path_exe,
                             database=MMSEQS_DB, format=generated_alis_formats)
        if run_mmseqs:
            if not skip_ali_generation:
                M.main()
        if generated_alis_formats == 'both':
            list_ali = [M.out_fasta, M.out_a3m]
        elif generated_alis_formats == 'a3m':
            list_ali = [M.out_a3m]
        elif generated_alis_formats == 'fasta':
            list_ali = [M.out_fasta]
        else:
            print("Please define a format in the choices both, a3m or fasta")
            sys.exit()
        for elt in list_ali:
            if elt:
                list_files_to_map.append(elt)
                rootname, iext = os.path.splitext(elt)
                if iext == '.a3m':
                    linkname = "sequence_{}".format(ii) + iext
                else:
                    linkname = "alisequence_{}".format(ii) + iext
                    linkname_FLsequences = "sequence_{}_FL".format(ii) + iext
                    if os.path.islink(linkname_FLsequences):
                        os.remove(linkname_FLsequences)
                    os.symlink(rootname+"_FL.fasta", linkname_FLsequences) # Should be in phase with the extension added in MapTaxonomy...
                if os.path.islink(linkname):
                    os.remove(linkname)
                os.symlink(rootname+'_withdef'+iext, linkname)

            if M.method == 'mmseqs_profile':
                # The tsv files are used to sort the sequence with respect to a score calculated from a combination of qid and cov
                list_files_tsv_qidcov.append(M.out_tsv)
        linkname = "alisequence_{}.tsv".format(ii)
        if os.path.islink(linkname):
            os.remove(linkname)
        os.symlink(M.out_tsv, linkname)

    print("List of tsv files : {}".format(list_files_tsv_qidcov))


    # With this command we only load the huge UniRef100 db once for all files
    all_input_files = ",".join(list_files_to_map)
    print("List of ali files to be processed: {}".format(list_files_to_map))
    if M.method == 'mmseqs_profile':
        all_tsv_files = ",".join(list_files_tsv_qidcov)
    if run_mmseqs:
        # Next step will generate msa_xxx/xxx_withdef.a3m and msa_xxx/xxx_withdef.fasta
        # Then it will generate msa_xxx/xxx_FL.fasta
        GetTaxFullSeqSortAli.main(all_input_files, output=False, full_length=True, file_qidcov=all_tsv_files)
        print("Finished with GetTaxFullSeq")
    os.chdir(cur_dir)

if __name__ == "__main__":
    usage = " GenerateMSA4EachPartner.py -i <directory with all the fasta to process> \n" + \
            "Description: All the files containing a fas or fasta extension will be processed"+ \
            "Multiple sequence alignment will be generated with the chosen method"
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument(
        '-i',
        '--input_dir_fasta',
        type=str,
        help="Working directory where all the .fas or .fasta input files are stored. "
             "Filters or list of files can be added to limit the set of files to consider"
    )
    parser.add_argument(
        '-f',
        '--file_name_filter',
        type=str,
        default="",
        help="Filter string used with wild-cards to recover all the fasta files to consider as input"
    )
    parser.add_argument(
        '-l',
        '--list_of_fasta_files',
        type=str,
        default="",
        help="List of files written "
    )
    parser.add_argument(
        '-m',
        '--method',
        choices=['hhblits', 'mmseqs', 'mmseqs_profile'],
        default="mmseqs_profile",
        help="Select the method to use to generate the alignments"
    )
    parser.add_argument(
        '-n',
        '--no_run_mmseqs',
        action='store_true',
        help="Blank run in which mmseqs is not run but links and files are regenerated. Used when you need to test new associations of already computed MSA."
    )
    options = parser.parse_args()
    if options.no_run_mmseqs:
        run_mmseqs = False
    else:
        run_mmseqs = True
    main(options.input_dir_fasta, options.method, list_files=options.list_of_fasta_files, name_filter=options.file_name_filter, run_mmseqs=run_mmseqs)
