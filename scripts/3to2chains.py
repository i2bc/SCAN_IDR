from argparse import ArgumentParser
from glob import glob
from biopandas.pdb import PandasPdb


def resnumber(d, res):
    return d[res]


def built_dict(prec,list_res):
    d={list_res[0] : prec}
    res_prec=list_res[0]
    for res in list_res[1:]:
        d[res]=prec+(res-res_prec)
        res_prec=res
        prec=prec+(res-res_prec)
    return prec+200, d


def changeChain(fileout, ppdb, c1, c2):
    pdb = ppdb.df['ATOM']
    print(pdb.head(10))
    print(list(pdb["residue_number"][pdb["chain_id"] == c1]))
    resmax = max(list(pdb["residue_number"][pdb["chain_id"] == c1]))
    print(resmax)
    pdb["residue_number"][pdb["chain_id"] == c2] = pdb["residue_number"][pdb["chain_id"] == c2].apply(lambda x:x+resmax+200)
    pdb["chain_id"][pdb["chain_id"] == c2] = pdb["chain_id"][pdb["chain_id"] == c2].apply(lambda x:c1)
    ppdb.df['ATOM'] = pdb
    ppdb.to_pdb(path=f'{fileout}.pdb', 
            records=['ATOM'], 
            append_newline=True)

def main():
    parser = ArgumentParser()
    parser.add_argument(
        "-c",
        default="A,B",
        help="name of the two chains sep by ','"
    )
    parser.add_argument(
        "-pdb",
        help="Path of the pdb"
    )
    args = parser.parse_args()
    c1,c2=args.c.split(',')
    ppdb = PandasPdb().read_pdb(args.pdb)
    fileout = args.pdb.replace('.pdb', '')+'_ch3to2'
    changeChain(fileout, ppdb, c1, c2)

if __name__ == '__main__':
    main()


"""

python scripts/3to2chains.py -c B,C -pdb /store/EQUIPES/AMIG/PROJECTS/COALI4AF/Recent_peptides/af2_runs/af2_6SAT_gapPeptide/af2_predictions_nopep_v1/MSAforAlphafold_mixed_noPep_unrelaxed_rank_1_model_1_ch3to2.pdb
python scripts/3to2chains.py -c B,C -pdb af2_runs/af2_6SAT_DELIM/af2_predictions_delim_v1/MSAforAlphafold_mixed_DELIM_unrelaxed_rank_2_model_2.pdb
python scripts/3to2chains.py -c B,C -pdb af2_runs/af2_6SAT_DELIM/af2_predictions_delim_v1/MSAforAlphafold_mixed_DELIM_unrelaxed_rank_3_model_3.pdb
python scripts/3to2chains.py -c B,C -pdb af2_runs/af2_6SAT_DELIM/af2_predictions_delim_v1/MSAforAlphafold_mixed_DELIM_unrelaxed_rank_4_model_4.pdb
python scripts/3to2chains.py -c B,C -pdb af2_runs/af2_6SAT_DELIM/af2_predictions_delim_v1/MSAforAlphafold_mixed_DELIM_unrelaxed_rank_5_model_5.pdb

python scripts/3to2chains.py -c A,B -pdb /store/EQUIPES/AMIG/PROJECTS/COALI4AF/Recent_peptides/data/ref_capri/6SAT_A-C.pdb

"""
