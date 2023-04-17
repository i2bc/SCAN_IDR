import os
import shutil
import glob
import sys
import argparse
import socket
import configparser
#import subprocess

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(os.path.dirname(script_path))

config = configparser.ConfigParser()
config.read(f'{script_dir}/config.ini')

MMSEQS_BIN = config['PROGRAMS']['MMSEQS_BIN']
MMSEQS_DB = config['PROGRAMS']['MMSEQS_DB']

class MMSEQS:

    def __init__(self,input_fasta, output_file=None, ali_method='mmseqs_profile',
                 executable_path=MMSEQS_BIN,
                 database=MMSEQS_DB, format='both'):
        self.input_fasta = input_fasta
        self.method = ali_method
        self.path_exe = executable_path
        self.format = format

        self.filename = os.path.basename(self.input_fasta)
        self.pathname = os.path.dirname(self.input_fasta)
        self.rootname, self.extension = os.path.splitext(self.filename)

        if database:
            self.db = database

        if output_file:
            self.resultdir = os.path.dirname(self.output_msa)
            self.output_msa = output_file
        else:
            self.resultdir = os.path.join(self.pathname, "msa_{}".format(self.rootname))
            if not os.path.exists(self.resultdir):
                os.mkdir(self.resultdir)
            self.output_msa = None # no extension at that stage

        self.define_output_path()

    def main(self):
        if self.method == 'mmseqs_profile':
            self.compute_msas_with_mmseqs_indexed_profiles()

        else:
            print("Sorry the methods you requested were not implemented yet ...")
            sys.exit()

    def define_output_path(self):
        print("Files which will be generated :")
        if not self.output_msa:
            if self.format == 'both' or self.format == 'fasta':
                self.out_fasta = os.path.join(self.resultdir, self.rootname) + ".fasta"
                print("\t", self.out_fasta)
            if self.format == 'both' or self.format == 'a3m':
                self.out_a3m = os.path.join(self.resultdir, self.rootname) + ".a3m"
                print("\t", self.out_a3m)
        else:
            self.output_ali = self.output_msa
            print("\t", self.output_ali)

        if self.method == 'mmseqs_profile':
            self.out_tsv = os.path.join(self.resultdir, self.rootname) + ".tsv"
            print("\t", self.out_tsv)


    def compute_msas_with_mmseqs_indexed_profiles(self, create_outtsv=True):
        """
        Define the steps to run the mmseqs algo using expandable db of profiles
        Returns: Up to 3 files are generated -> self.out_fasta, self.out_a3m and self.out_tsv

        """

        qdb = self.create_querydb()
        res, tmp = self.search_homologs(qdb, self.db, num_iterations=3, db_load_mode=2,
                                   sensitivity=8, e_value=0.1, max_seqs=10000)

        dbidx = self.db + ".idx"

        res_exp = self.expand_ali(qdb, res, dbidx, db_load_mode=2, expansion_mode=0, expand_eval=float("inf"),
                                  expand_filter_clusters=1, max_seq_id=1.0)

        last_profile_origin = os.path.join(tmp, "latest", "profile_1")
        last_profile_target = os.path.join(self.resultdir, "prof_res")
        self.mv_db(last_profile_origin, last_profile_target)

        header_file_origin = os.path.join(self.resultdir, "qdb_h")
        header_file_target = os.path.join(self.resultdir, "prof_res_h")
        self.ln_db(header_file_origin, header_file_target)


        res_exp_realign = self.realign_homologs(last_profile_target, dbidx, res_exp, db_load_mode=2,
                              max_accept=100000, align_eval=10, alt_ali=10)

        res_exp_realign_filter = self.filter_result(qdb, dbidx, res_exp_realign, db_load_mode=2, qid=0, qsc=0.8, diff=0,
                                                    max_seq_id=1.0, filter_min_enable=100)

        if self.format == 'both' or self.format == 'fasta':
            self.result2msa(qdb, dbidx, res_exp_realign_filter, self.out_fasta, msa_format_mode='fasta', db_load_mode=2,
                            qid=None, qsc=None, max_seq_id=1.0, filter_min_enable=None, diff=None)
            print("{} was generated".format(self.out_fasta))

        if self.format == 'both' or self.format == 'a3m':
            self.result2msa(qdb, dbidx, res_exp_realign_filter, self.out_a3m, msa_format_mode='a3m', db_load_mode=2,
                            qid=None, qsc=None, max_seq_id=1.0, filter_min_enable=None, diff=None)
            print("{} was generated".format(self.out_a3m))
        if create_outtsv:
            selected_tabs = "query,target,pident,qcov,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,cigar,taxid,taxname,taxlineage,theader,tseq"
            self.convertalis(last_profile_target, dbidx, res_exp_realign_filter, self.out_tsv, db_load_mode=2,
                                        format_mode=4, tabs=selected_tabs)
            print("{} was generated".format(self.out_tsv))
        else:
            self.out_tsv = None

        print("Running :\nCleaning up of all temporary and intermediate files")
        for elt in [res, res_exp, res_exp_realign, res_exp_realign_filter]:
            self.remove_files(elt)
        for elt in [qdb, header_file_origin]:
            self.remove_files(elt)
        shutil.rmtree(tmp, ignore_errors=True)
        l_last_profile_target = glob.glob(last_profile_target+"*")
        print("prof_res files to remove:",l_last_profile_target)
        for elt in l_last_profile_target:
            os.system("rm {}".format(elt))
            #shutil.rmtree(elt, ignore_errors=True)


    def create_querydb(self):
        """
        Format the query as a mmseqs db
        Typical cmd:
        "${MMSEQS}" createdb "${QUERY}" "${BASE}/qdb"
        Returns: qdb name
        """
        qdb = os.path.join(self.resultdir,"qdb")
        exe = ["{} createdb".format(self.path_exe)]
        args = [
            self.input_fasta,
            qdb
        ]
        cmd = " ".join(exe+args)
        print("Running :\n{}".format(cmd))
        os.system(cmd)
        return qdb

    def search_homologs(self, qdb, db, num_iterations=3, db_load_mode=2, sensitivity=8, e_value=0.1, max_seqs=10000):
        """
        Search qdb against the mmseqs db
        typical cmd:
        "${MMSEQS}" search "${RES}/qdb" "${DBBASE}/${DB1}" "${RES}/res" "${BASE}/tmp" $SEARCH_PARAM
        with : SEARCH_PARAM="--num-iterations 3 --db-load-mode 2 -a -s 8 -e 0.1 --max-seqs 10000"
        Returns: res (stands for 'results not expanded')
                 tmp (path to directory in which some files can be exploited to expand clusters)
        """

        res = os.path.join(self.resultdir,"res")
        tmp = os.path.join(self.resultdir,"tmp")
        exe = ["{} search".format(self.path_exe)]
        args = [
            qdb,
            db,
            res,
            tmp,
            "--num-iterations {}".format(num_iterations),
            "--db-load-mode {}".format(db_load_mode),
            "-a",
            "-s {}".format(sensitivity),
            "-e {}".format(e_value),
            "--max-seqs {}".format(max_seqs),
        ]
        cmd = " ".join(exe + args)
        print("Running :\n{}".format(cmd))
        os.system(cmd)
        return res, tmp

    def expand_ali(self, qdb, res, db, db_load_mode=2, expansion_mode=0, expand_eval=float("inf"),
                   expand_filter_clusters=1, max_seq_id=1.0):
        """
        Expand homolog detection from head of clusters to all sequences
        mmseqs expandaln <i:queryDB> <i:targetDB> <i:resultDB> <i:resultDB|ca3mDB> <o:alignmentDB> [options]

        Typical cmd :
        "${BASE}/qdb" "${DBBASE}/${DB1}.idx" "${BASE}/res" "${DBBASE}/${DB1}.idx" "${BASE}/res_exp" --db-load-mode
        2 ${EXPAND_PARAM}
        with : EXPAND_PARAM="--expansion-mode 0 -e ${EXPAND_EVAL} --expand-filter-clusters ${FILTER} --max-seq-id 1.0"
         --expansion-mode INT : Update score, E-value, and sequence identity by 0: input alignment 1: rescoring the inferred backtrace [0]
         --expand-filter-clusters INT : Filter each target cluster during expansion 0: no filter 1: filter [0]
        Returns: res_exp (stands for 'results expanded')
        """
        res_exp = os.path.join(self.resultdir,"res_exp")

        exe = ["{} expandaln".format(self.path_exe)]
        args = [
            qdb,
            db,
            res,
            db,
            res_exp,
            "--db-load-mode {}".format(db_load_mode),
            "--expansion-mode {}".format(expansion_mode),
            "-e {}".format(expand_eval),
            "--expand-filter-clusters {}".format(expand_filter_clusters),
            "--max-seq-id {}".format(max_seq_id),
        ]
        cmd = " ".join(exe + args)
        print("Running :\n{}".format(cmd))
        os.system(cmd)
        return res_exp

        pass

    def mv_db(self, origin, target):
        """
        Move dbs from some directories
        """
        exe = ["{} mvdb".format(self.path_exe)]
        args = [
            origin,
            target,
        ]
        cmd = " ".join(exe+args)
        print("Running :\n{}".format(cmd))
        os.system(cmd)

    def ln_db(self, origin, target):
        """
        Create links for instance to bind header files
        """
        exe = ["{} lndb".format(self.path_exe)]
        args = [
            origin,
            target,
        ]
        cmd = " ".join(exe+args)
        print("Running :\n{}".format(cmd))
        os.system(cmd)

    def realign_homologs(self, profile, db, cluster_heads, db_load_mode=2, max_accept=1000000, align_eval=10, alt_ali=10):
        """
        Realignment of the sequences detected after the profile search
        Typical command: mmseqs align <i:queryDB> <i:targetDB> <i:resultDB> <o:alignmentDB> [options]
       "${MMSEQS}" align "${BASE}/prof_res" "${DBBASE}/${DB1}.idx" "${BASE}/res_exp" "${BASE}/res_exp_realign"
         --db-load-mode 2 -e ${ALIGN_EVAL} --max-accept ${MAX_ACCEPT} --alt-ali 10 -a
        with :
                ALIGN_EVAL=10  List matches below this E-value (range 0.0-inf) [1.000E-03]
                DIFF=3000
                --alt-ali Show up to this many alternative alignments
                MAX_ACCEPT=1000000
                -e : e-value ()
                -a is for  Add backtrace string (convert to alignments with mmseqs convertalis module)
        Args:
            profile: latest profile calculated which was recovered from the tmp directory
            db: path to the database
            cluster_heads: correspond to the output of the expanded alignment
            db_load_mode: 2 -> keep in memory

        Returns: res_exp_realign (alignment realigned with a Smith-Waterman backtrack)

        """
        res_exp_realign = os.path.join(self.resultdir, "res_exp_realign")
        exe = ["{} align".format(self.path_exe)]
        args = [
            profile,
            db,
            cluster_heads,
            res_exp_realign,
            "--db-load-mode {}".format(db_load_mode),
            "--max-accept {}".format(max_accept),
            "-e {}".format(align_eval),
            "--alt-ali {}".format(alt_ali),
            "-a",
        ]
        cmd = " ".join(exe + args)
        print("Running :\n{}".format(cmd))
        os.system(cmd)
        return res_exp_realign

    def filter_result(self, qdb, db, res_exp_align, db_load_mode=2, qid=0, qsc=-20.0, diff=None,
                      max_seq_id=1.0, filter_min_enable=1000000):
        """
        Method to filter the results of an alignment using a variety of possible filters
        mmseqs filterresult <i:queryDB> <i:targetDB> <i:resultDB> <o:resultDB> [options]

        Typical cmd : "${MMSEQS}" filterresult "${BASE}/qdb" "${DBBASE}/${DB1}.idx" "${BASE}/res_exp_realign"
        "${BASE}/res_exp_realign_filter"
        --filter-min-enable INT   Only filter MSAs with more than N sequences, 0 always filters [0]
        Returns: path of the output file -> res_exp_realign_filter
        """
        #print("DEBUG")
        #print("{}-{}-{}-{}-{}-{}".format(db_load_mode, qid, qsc, diff, max_seq_id, filter_min_enable))

        res_exp_realign_filter = os.path.join(self.resultdir, "res_exp_realign_filter")
        exe = ["{} filterresult".format(self.path_exe)]
        args = [
            qdb,
            db,
            res_exp_align,
            res_exp_realign_filter,
            "--db-load-mode {}".format(db_load_mode),
            "--qid {}".format(qid),
            "--qsc {}".format(qsc),
            "--max-seq-id {}".format(max_seq_id),
        ]
        if diff is not None:
            print("filter_result uses diff")
            args.append("--diff {}".format(diff))
        if filter_min_enable is not None:
            args.append("--filter-min-enable {}".format(filter_min_enable))
        cmd = " ".join(exe + args)
        print("Running :\n{}".format(cmd))
        os.system(cmd)
        return res_exp_realign_filter

    def result2msa(self, qdb, db, res, filename_output, msa_format_mode='fasta', db_load_mode=2, qid=None, qsc=None, diff=None,
                   max_seq_id=1.0, filter_min_enable=None):
        """
        Format the output alignment into either fasta or a3m format
        mmseqs result2msa <i:queryDB> <i:targetDB> <i:resultDB> <o:msaDB> [options]

        Typical cmd: ${MMSEQS} result2msa "${BASE}/qdb" "${DBBASE}/${DB1}.idx" "${BASE}/res_exp_realign_filter"
           "${BASE}/${ROOT_NAME}_uniref.fasta" --msa-format-mode 2 --db-load-mode 2 ${FILTER_PARAM}
        with : FILTER_PARAM="--filter-msa ${FILTER} --filter-min-enable 1000 --diff ${DIFF}
            --qid 0.0,0.2,0.4,0.6,0.8,1.0 --qsc 0 --max-seq-id 0.95"
        Returns: path of the output file -> filename_output
        """
        d = dict([('fasta',2),('a3m',5)])

        exe = ["{} result2msa".format(self.path_exe)]
        args = [
            qdb,
            db,
            res,
            filename_output,
            "--db-load-mode {}".format(db_load_mode),
            "--max-seq-id {}".format(max_seq_id),
            "--msa-format-mode {}".format(d[msa_format_mode])
        ]
        if qid is not None:
            args.append("--qid {}".format(qid))
        if qsc is not None:
            args.append("--qsc {}".format(qsc))
        if diff is not None:
            args.append("--diff {}".format(diff))
        if filter_min_enable is not None:
            args.append("--filter-min-enable {}".format(filter_min_enable))
        cmd = " ".join(exe + args)
        print("Running :\n{}".format(cmd))
        os.system(cmd)

    def convertalis(self, prof, db, res, tab_filename_output, db_load_mode=2, format_mode=4,
                    tabs="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"):
        """
         Generates an output file summarizing the properties of the alignment into a table
         mmseqs convertalis <i:queryDb> <i:targetDb> <i:alignmentDB> <o:alignmentFile> [options]

        Typical cmd : "${MMSEQS}" convertalis "${BASE}/prof_res" "${DBBASE}/${DB1}.idx" "${BASE}/res_exp_realign_filter"
         "${BASE}/res_exp_realign_filter.m8"
         --format-output query,target,pident,qcov,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,cigar
         --format-mode 4 --db-load-mode 2
        with :
        --format-output STR : Choose comma separated list of output columns from: query,target,evalue,gapopen,pident,
        fident,nident,qstart,qend,qlen, tstart,tend,tlen,alnlen,raw,bits,cigar,qseq,tseq,qheader,theader,qaln,taln,
        qframe,tframe,mismatch,qcov,tcov, qset,qsetid,tset,tsetid,taxid,taxname,taxlineage,qorfstart,qorfend,
        torfstart,torfend [query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits]
         --format-mode INT Output format:
                           0: BLAST-TAB
                           1: SAM
                           2: BLAST-TAB + query/db length
                           3: Pretty HTML
                           4: BLAST-TAB + column headers
                           BLAST-TAB (0) and BLAST-TAB + column headers (4) support custom output formats (--format-output) [0]
        Returns: path of the output file -> blast-tab filename output
        """

        exe = ["{} convertalis".format(self.path_exe)]
        args = [
            prof,
            db,
            res,
            tab_filename_output,
            "--db-load-mode {}".format(db_load_mode),
            "--format-output {}".format(tabs),
            "--format-mode {}".format(format_mode)
        ]

        cmd = " ".join(exe + args)
        print("Running :\n{}".format(cmd))
        os.system(cmd)


    def remove_files(self, db_file):
        """
        "${MMSEQS}" rmdb "${BASE}/res_exp_realign"

        Returns: None
        """
        exe = ["{} rmdb".format(self.path_exe)]
        args = [
            db_file,
        ]
        cmd = " ".join(exe+args)
        print("Running :\n{}".format(cmd))
        os.system(cmd)


if __name__ == "__main__":
    usage = " RunMmseqs.py -i <fasta_input> \n" + \
            "Description: All the files containing a fas or fasta extension will be processed"+ \
            "Multiple sequence alignment will be generated with the chosen method"
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument(
        '-i',
        '--input_file',
        type=str,
        help="Path to the fasta input file"
    )
    parser.add_argument(
        '-o',
        '--output_file',
        type=str,
        default=None,
        help="Path to the alignment output file. Not required. By default a directory msa_<root_name> "
             "will be created containing the msa in format as defined in format parameter"
    )
    parser.add_argument(
        '-f',
        '--format',
        choices=['a3m', 'fasta', 'both'],
        help="Select either fasta, a3m or both formats as outputs"
    )
    parser.add_argument(
        '-m',
        '--method',
        choices=['hhblits', 'mmseqs', 'mmseqs_profile'],
        default='mmseqs_profile',
        help="Select which method to use to generate the alignments. "
             "Mmseqs should be installed and memory capacities checked."
             "Current default is mmseqs_prof."
    )
    parser.add_argument(
        '-p',
        '--pathexe',
        default=None,
        help="Path to the executable, example /alpha/programs/gits/MMseqs2/ on node47"
    )
    parser.add_argument(
        '-d',
        '--db',
        default=None,
        help="Path to the database, example /alpha/database/uniref30_2202/uniref30_2202_db"
    )

    options = parser.parse_args()

    hostname = socket.gethostname()
    # Get/Define the executable for mmseqs
    if not options.pathexe and hostname in ['node47']:
        path_exe = "/alpha/programs/gits/MMseqs2/bin/mmseqs"
    else:
        print("You need to define a valid program on the node where you wish to compute.")
        path_exe = None

    # Get/Define the database on which to search
    if not options.db and hostname in ['node47']:
        path_db = "/alpha/database/uniref30_2202"
        db = os.path.join(path_db,"uniref30_2202_db")
    elif options.db:
        db = options.db
    else:
        print("You need to define a valid db on the node where you wish to compute.")
        raise IOError

    path_exe = options.pathexe

    if path_exe:
        M = MMSEQS(options.input_fasta, output_file=options.output_file, ali_method=options.method,
                   executable_path=path_exe, database=db, format=options.format)
        M.main()


