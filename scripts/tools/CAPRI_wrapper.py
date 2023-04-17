
'''
Created on 28 mar. 2014

@author: Jinchao Yu

This script is used for evaluating a large amount of complexes by calling CAPRI.py

'''
import sys
import os
import tempfile
import optparse
import glob
try:
    import commands
except ModuleNotFoundError:
    import subprocess as commands
import math
import re
import shutil
from collections import defaultdict

import configparser
#import subprocess

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(os.path.dirname(script_path))

config = configparser.ConfigParser()
config.read(f'{script_dir}/config.ini')
sys.path.append(script_path)

CAPRI_PY = os.path.join(script_dir, "tools", "CAPRI.py")
PYTHON_DEF = config['PROGRAMS']['PYTHON']

### tool box
def file_columns_2_defaultdict(list_file, key_column_number, value_column_number_list, delimiter="[\t ]+"):
    """Read a file containing a table and transfer one column to a key and multiple columns to a defaultdict value
    Any line beginning with # and any line containing only spaces and/or '\t' will be ignored
    """
    my_dict = defaultdict(list)
    print_warning = False
    try:
        FH = open(list_file)
    except:
        print("Input file not found.")
        return my_dict
    for line in FH:
        res_match = re.match("#",line)
        if res_match:
            continue
        res_match = re.match("^[\s\t]*$",line)
        if res_match:
            continue
        try:
            terms = re.split(delimiter,line.rstrip())
            my_key = terms[key_column_number-1]
            my_value_list = []
            for value_column_number in value_column_number_list:
                try:
                    my_value_list.append(terms[value_column_number-1])
                except:
                    print_warning = True
                    pass
            my_dict[my_key] = my_value_list
        except:
            print("Warning: this line is abnormal in the file %s"%list_file)
            print(line)
    if print_warning:
        print("File {} has less columns than expected".format(list_file))
    FH.close()
    return my_dict

def file_column_2_list(list_file, column_number, delimiter="[\t ]+"):
    """Read a file containing a table and transfer one column to a list
    Any line beginning with # and any line containing only spaces and/or \t will be ignored
    """
    my_list = []
    try:
        FH = open(list_file)
    except:
        print("Input file not found.")
        return my_list
    for line in FH:
        res_match = re.match("#",line)
        if res_match:
            continue
        res_match = re.match("^[\s\t]*$",line)
        if res_match:
            continue
        try:
            my_term = re.split(delimiter,line.rstrip())[column_number-1]
        except:
            print("Warning: this line is abnormal in the file %s"%list_file)
        my_list.append(my_term)
    FH.close()
    return my_list

def combine_CAPRI_results(CAPRI_file_lst, CAPRI_combined_outfile,reference_info):
    """# file RMSAvsA RMSBvsB RMSABvsAB LRMSA LRMSB FRNAT FRNNAT FRIR IRMS CAPRIrank
    """
    dict_rank_int = {"Incorrect":0,"Acceptable":1,"Medium":2,"High":3}
    lst_dict = []
    lst_decoy_names = set()
    for i in range(len(CAPRI_file_lst)):
        filepath = CAPRI_file_lst[i]
        print(filepath)
        print('')
        lst_dict.append(file_columns_2_defaultdict(filepath, 1, range(1,14)))
        print("Extraction of CAPRI results dict:")
        print(list(file_columns_2_defaultdict(filepath, 1, range(1,14)).items())[:10])
        FH = open(filepath)
        print(FH.readline())
        FH.close()
        # we have to include decoy names from all different references when quick_fnat is used
        lst_decoy_names.update(file_column_2_list(CAPRI_file_lst[i], 1))
    # merge and sort the list of decoy names
    lst_decoy_names = sorted(list(lst_decoy_names), key=get_Zdock_decoy_number)

    if len(CAPRI_file_lst)>1:
        content = "# multiple reference files: %s\n"%str(reference_info)
    else:
        content = ""      

    if len(lst_decoy_names) == 0: # no decoys in list (happens when pre-filtering with quick_fnat first and all decoys are bad)
        content += "# file RMSAvsA RMSBvsB RMSABvsAB LRMSA LRMSB FRNAT FRNNAT FRIR IRMS CAPRIrank DockQ DockQrank ref_index\n"
        FHo = open(CAPRI_combined_outfile,'w')
        FHo.write(content)
        FHo.close()
        return

    if len(lst_dict[0][lst_decoy_names[0]]) == 11:
        content += "# file RMSAvsA RMSBvsB RMSABvsAB LRMSA LRMSB FRNAT FRNNAT FRIR IRMS CAPRIrank ref_index\n"
    else: # len == 13 -> DockQ score and DockQ category
        content += "# file RMSAvsA RMSBvsB RMSABvsAB LRMSA LRMSB FRNAT FRNNAT FRIR IRMS CAPRIrank DockQ DockQrank ref_index\n"
    for decoy_name in lst_decoy_names:
        rank, irms = -1, 999999
        top = -1
        for i in range(len(CAPRI_file_lst)): # if several references
            # try:
            #     print "lst_dict[i][decoy_name][10]"
            #     print lst_dict[i][decoy_name][10]
            # except:
            #     print lst_dict[1]
            #     print lst_dict[0]
            #     print decoy_name,i
            #     print lst_dict[i][decoy_name]
            #     print len(lst_dict[i][decoy_name])
            try:
                my_rank, my_irms = dict_rank_int[lst_dict[i][decoy_name][10]],float(lst_dict[i][decoy_name][9])
                if my_rank > rank or (my_rank == rank and my_irms < irms):
                    top = i
                    rank, irms = my_rank, my_irms
            except:
                print("Warning: Problem with decoy %s"%decoy_name)
                continue
        if top != -1:
            content += " ".join(lst_dict[top][decoy_name]) + " %d \n"%top
    FHo = open(CAPRI_combined_outfile,'w')
    FHo.write(content)
    FHo.close()

def div_list(lst, divisor):
    """In a mpi program, when a complete list of tasks should be divided for each CPU,
    use div_list(lst, size)[rank]
    """
    size = int(math.ceil(float(len(lst))/divisor))
    my_list = [lst[i:i+size] for i in range(0, len(lst), size)]
    diff = divisor - len(my_list)# probably void sub-list should be added
    for i in range(diff):
        my_list.append([])
    return my_list

def get_pdb_file_list(dir_abspath):
    """Get the list of files ending at ".pdb" in the given directory.
    """
    pdb_file_list = []
    for path in os.listdir(dir_abspath):
        res_search = re.search("\.pdb$",path)
        if res_search:
            pdb_file_list.append(path)
    return pdb_file_list

def file_2_list(list_file):
    """Read a file containing a list, whose each term is located in a line
    Any line beginning with # and any line containing only spaces or \t will be ignored
    """
    my_list = []
    try:
        FH = open(list_file)
    except:
        print("Input file not found.")
        return my_list
    for line in FH:
        res_match = re.match("#",line)
        if res_match:
            continue
        res_match = re.match("^[\s\t]*$",line)
        if res_match:
            continue
        my_list.append(line.rstrip())
    FH.close()
    return my_list

def make_tempfile(suffix_str='',prefix_str='',parent_dir=None):
    """ create a tempfile in parent_dir (def = $HOME/tmp). Return the created tempfile absolute path.
    File handle problem has been solved.
    Arguments:
    1. suffix_str: string for tempfile name's suffix part (def = '')
    2. prefix_str: string for tempfile name's prefix part (def = '')
    3. parent_dir: parent directory path (def = $HOME/tmp)
    Return:
    the created tempdir absolute absolute path
    """
    if parent_dir == None:
        parent_dir = os.path.join(os.getenv("HOME"),"tmp")
        if not os.path.exists(parent_dir):
            try:
                os.mkdir(parent_dir)
            except:
                print("The program needs to create a tmp dir at $HOME/tmp/.")
                print("Creation failed.")
                raise IOError
        if not os.access(parent_dir, os.W_OK):
            print("No writing right in %s"% parent_dir)
            raise IOError
    fd, tmpfile = tempfile.mkstemp(suffix=suffix_str,prefix=prefix_str,dir=parent_dir)
    os.close(fd)
    return os.path.abspath(tmpfile)

def get_Zdock_decoy_number(decoy_name):
    """ used for sort the zdock decoys
    """
    res_search = re.search("[^\d](\d+)\.pdb",decoy_name)
    if res_search:
        return int(res_search.group(1))
    else:
        return -1

def extract_rank(name):
    res_rank = re.match("tmp_output_file_(\d+)",name)
    if res_rank:
        return int(res_rank.group(1))
    else:
        return 10000

def combine_multiple_output_files(l_output_files, global_ouput_filename):

    with open(global_ouput_filename,"w") as fout:
        for jj, file_out_tmp in enumerate(l_output_files):
            ftmp = open(file_out_tmp)
            for ll in ftmp:
                if jj == 0 and ll[0]=="#":
                    fout.write(ll)
                elif jj > 0 and ll[0]=="#":
                    continue
                else:
                    fout.write(ll)
            os.remove(file_out_tmp)
    fout.close()

def Main():
    """Main function is used on  server for launching a qsub job by calling itself with "mpi" argument
    """
    USAGE = "CAPRI_wrapper.py -r reference.pdb -t target.pdb" +\
            " -o outfile -c AB:XY \n"+\
            "    Option -t allows to evaluate only one structure locally.\n"+\
            "RMK: /!\ be careful, atom names have to start in column 14 (i.e. with a space)\n"+\
            "if name has less than 4 letters/numbers otherwise ProFit won't be able to select\n"+\
            "CA, C, O and N atoms and you'll end up with en error"
    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option('-r','--reference',action="store",dest='reference',type='string',
                      help="reference PDB file(s) in evaluation, multiple ref files "+\
                      "are given like ref1.pdb,ref2.pdb,ref3.pdb. In this case, the results with "+\
                      "the best rank and irms (according to diff ref) will be reported")
    parser.add_option('-t','--target_pdb',action="store",dest='target_pdb',type='string',
                      help="the target PDB file (only one) to be evaluated (def=None). ", default=None)
    parser.add_option('-d','--directory_pdb',action="store",dest='directory_pdb',type='string',
                      help="the directory containing target PDB files to be evaluated (def=None). ", default=None)
    parser.add_option('-c','--chain',action="store",dest='chain',type='string',
                      help="chain correspondence, for example, AB:XY means chain A and B in reference "+\
                      "correspond respectively chain X and Y of the models to be evaluated")
    parser.add_option('-o','--outfile',action="store",dest='outfile',type='string',
                      help="outfile, def = capri_evaluation.txt in current directory",
                      default = os.path.join(os.path.abspath(os.path.curdir), "capri_evaluation.txt"))
    parser.add_option('--quick_fnat',action="store_true",dest="quick_fnat",default=False,
                      help="use quick_fnat to prefilter decoys and accelerate calculation")
    parser.add_option('--quick_fnat_ref',action="store",dest='quick_fnat_ref',type='string',
                      help="optional reference PDB file(s) for quick_fnat evaluation (native unbound)."+\
                      "Multiple ref files can be given as unboundref1.pdb,unboundref2.pdb,unboundref3.pdb. "+\
                      "They MUST BE in the same order as the reference PDB files in the -r option.")
    parser.add_option('-w','--workdir',action="store",dest='workdir',type='string',
                      help="workdir where .o and .e files will be stored, def = current directory",
                      default = os.getcwd())
    parser.add_option('-a','--no_absolute_path',action="store_true",dest='no_absolute_path',default=False,
                      help="Do not use absolute path filenames in case they are too long and are deleterious for proper running",
                      )


    if len(sys.argv) == 1:
        print(USAGE)
        print("type -h or --help for more help information")
        sys.exit(1)
    (options, args) = parser.parse_args(sys.argv)
    if len(args) > 1: # sys.argv[0] is the name of this program
        print("Leftover arguments:" + str(args[1:]))
    if options.no_absolute_path:
        options.pdb_dir = os.path.curdir
        options.workdir = options.workdir
        
    else:
        options.pdb_dir = os.path.abspath(os.path.curdir)
        options.workdir = os.path.abspath(options.workdir)

    # serial job
    if options.target_pdb:
        l_target_pdb = [options.target_pdb]
    elif options.directory_pdb:
        l_target_pdb = glob.glob(os.path.join(options.directory_pdb, '*.pdb'))
    else:
        print(USAGE)
        print("type -h or --help for more help information")
        sys.exit(1)
    list_outputfile_indexed = []
    for ii, target_pdb in enumerate(l_target_pdb):
        ref_lst = []
        outputfile_indexed = options.outfile + f"__{ii}"
        list_outputfile_indexed.append(outputfile_indexed)
        if ',' in options.reference:
            lst = re.split(",",options.reference)
            for name in lst:
                if options.no_absolute_path:
                    ref_lst.append(name)
                else:
                    ref_lst.append(os.path.abspath(name))
        else:
            if options.no_absolute_path:
                ref_lst = [options.reference]
            else:
                ref_lst = [os.path.abspath(options.reference)]
        if options.quick_fnat_ref:
            if not options.quick_fnat:
                print("WARNING: quick_fnat_ref files will be ignored since quick_fnat option is not provided!")
            quick_fnat_ref_lst = []
            quick_fnat_lst = re.split(",",options.quick_fnat_ref)
            if len(quick_fnat_lst) != len(lst):
                sys.exit("The same number of quick_fnat_ref files as reference files should be given (and files should be in the same order).")
            for name in quick_fnat_lst:
                quick_fnat_ref_lst.append(os.path.abspath(name))
        else:
            quick_fnat_ref_lst = ref_lst
        tmp_final_file_lst = []
        for i in range(len(ref_lst)):
            reference_pdb = ref_lst[i]
            tmp_final_file = make_tempfile("_%d"%i)
            tmp_final_file_lst.append(tmp_final_file)
            cmd_CAPRI = PYTHON_DEF + " " + CAPRI_PY + " -ref " + reference_pdb + " -mob " + target_pdb +\
                    " -chain " + options.chain + " -f " + tmp_final_file
            if options.quick_fnat:
                quick_fnat_ref_pdb = quick_fnat_ref_lst[i]
                cmd_CAPRI += " -quick_fnat -quick_fnat_ref " + quick_fnat_ref_pdb
            if options.no_absolute_path:
                cmd_CAPRI += " --no_absolute_path "
            (stat, output) = commands.getstatusoutput(cmd_CAPRI)
            if stat != 0:
                print("status: " + str(stat))
                print(output)
        combine_CAPRI_results(tmp_final_file_lst, outputfile_indexed, options.reference)
        for filename in tmp_final_file_lst:
            os.remove(filename)
            print(filename)
        print("Output file is %s"%outputfile_indexed)
    combine_multiple_output_files(list_outputfile_indexed, options.outfile)


if __name__ == "__main__":
    # preparation
    if len(sys.argv) == 1:
        print("type -h or --help for usage.")
        sys.exit(1)
    print(CAPRI_PY)
    if not os.path.exists(CAPRI_PY):
        print("CAPRI.py not found")
        sys.exit(1)


    Main()
