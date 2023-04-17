#!/bin/bash

ROOT_NAME1=XERC_COLI
ROOT_NAME2=XERD_COLI

PROTEIN1=${ROOT_NAME1}.fasta
PROTEIN2=${ROOT_NAME2}.fasta

STOICHIOMETRY="1,1"

DO_BUILD_MSA=true
DO_RUN_WITH_MMSEQS=true
DO_RUN_WITH_HHBLITS=false
DO_FILTER_ALI=true
DO_FETCH_FL_SEQS=false
DO_FILTER_MMSEQS_WITH_FULL_LENGTH_SEQS=true
DO_RUN_MAFFT=true
DO_MATCH_SPECIES=true
DO_PAIR_SEQ_DELIM=true

FULL_LENGTH_OPTION='-FL'
PARALOG_FREE_OPTION='-pf'
RUN_MIXED_ALI_OPTION=true
RUN_PAIRED_ALI_OPTION=true
RUN_ALLUNPAIRED_ALI_OPTION=true

QID1=45
COV1=50
QID2=45
COV2=50
DELIM1="_"
DELIM2="_"

TERM1=">"
TERM2="XerC"
TERM3="XerD"

IT=1
NCPUS=8
MAXRES=20001
ID=100
MAX_NSEQ_IN_MAFFT=10000
TOPLINES=$((${MAX_NSEQ_IN_MAFFT} * 2))
echo Only the top ${MAX_NSEQ_IN_MAFFT} sequences will be processed by mafft
MAX_SEQ_HHBLITS=100000
# Not sure we need to increase this that much:
MAX_ALI=30000

PYTHON=python
#PROG=/store/EQUIPES/AMIG/PROGRAMMES/
#HHBLITS_BIN=${PROG}/git/hh-suite/build/bin/
#HHBLITS_DB=/store/EQUIPES/AMIG/DATABASES/UniRef30_2020_06/UniRef30_2020_06

SCRIPTS=.
ALPHASCRIPTS=${SCRIPTS}/inputAFscripts/
WORKING_DIRECTORY=.

index_partners=(0 1)
input_list=(${PROTEIN1} ${PROTEIN2})
rootname_list=(${ROOT_NAME1} ${ROOT_NAME2})
qid_list=(${QID1} ${QID2})
cov_list=(${COV1} ${COV2})
delim_list=(${DELIM1} ${DELIM2})


if ${DO_RUN_WITH_MMSEQS}
then
  if ${DO_BUILD_MSA}
  then
      ${PYTHON} ${ALPHASCRIPTS}/GenerateMSA4EachPartner.py -i ${WORKING_DIRECTORY} -l ${PROTEIN1},${PROTEIN2} -m mmseqs_profile
  else
      # We only change the links to the proper files. Keep the previous alignments
      ${PYTHON} ${ALPHASCRIPTS}/GenerateMSA4EachPartner.py -i ${WORKING_DIRECTORY} -l ${PROTEIN1},${PROTEIN2} -m mmseqs_profile --no_run_mmseqs
fi

for i in ${index_partners[@]}; do

    PROTEIN=${input_list[$i]}
    ROOTNAME=${rootname_list[$i]}
    QID=${qid_list[$i]}
    COV=${cov_list[$i]}

    echo RUNNING PROTEIN:${PROTEIN}

    cd ${WORKING_DIRECTORY}

    #######################################################################################################################################################################
    #FOR EACH PROTEIN SEPARATELY
    #######################################################################################################################################################################

    rm -f sequence_${i}.fasta
    ln -s ${PROTEIN} sequence_${i}.fasta

    if ${DO_BUILD_MSA} && ${DO_RUN_WITH_HHBLITS}
    then

    # 1 - Run HHblits
    echo 'hhblits'
    #${HHBLITS_BIN}/hhblits -i sequence_${i}.fasta -o sequence_${i}.hhr -oa3m sequence_${i}.a3m -n ${IT} -d ${HHBLITS_DB} -cpu ${NCPUS} -all -B ${MAX_ALI} -Z ${MAX_ALI} -maxseq ${MAX_SEQ_HHBLITS} -maxfilt ${MAX_SEQ_HHBLITS} -realign_max ${MAX_SEQ_HHBLITS} -maxres ${MAXRES}

    fi
    
    if ${DO_FILTER_ALI}
    then
    
    # 2 - Run HHfilter
        #######################################################################################
                #USER INPUTS : PARAMETERS OF HHFILTER 
        #######################################################################################

    # ${HHBLITS_BIN}/hhfilter -i sequence_${i}.a3m -o sequence_${i}_filtered.a3m -id ${ID} -cov ${COV} -qid ${QID}


    # This script is not running on full-length sequences but on the output alignment
    ${PYTHON} ${ALPHASCRIPTS}/toolsFilterMmseqsOut.py -i alisequence_${i}.tsv -a alisequence_${i}.fasta -f sequence_${i}_filtered.fasta -c ${COV} -q ${QID}

    if ${DO_RUN_WITH_MMSEQS}
    then
       cp sequence_${i}_filtered.fasta msa_${ROOTNAME}/${ROOTNAME}_filtered.fasta
    fi

    # Test the terms required for decision of the best filter
    echo grep -c ${TERM1} sequence_${i}_filtered.fasta
    grep -c  ${TERM1} sequence_${i}_filtered.fasta
    echo grep -c ${TERM2} sequence_${i}_filtered.fasta
    grep -c  ${TERM2} sequence_${i}_filtered.fasta
    echo grep -c ${TERM3} sequence_${i}_filtered.fasta
    grep -c  ${TERM3} sequence_${i}_filtered.fasta

    fi
    
    #head -n 1000 sequence_${i}_filtered.a3m > sequence_${i}_filtered_head1000.a3m
    #if head => change in python /store/EQUIPES/AMIG/SCRIPTS/workspace/sequence/RetrieveIDMappingUniprot.py -i sequence_${i}_filtered.a3m -o sequence_${i}_filtered_FL.fasta

    if ${DO_FETCH_FL_SEQS}
    then
    # 3 - Retrieve sequence FL
    ${PYTHON} ${SCRIPTS}/RetrieveIDMappingUniprot.py -i sequence_${i}_filtered.a3m -o sequence_${i}_filtered_FL.fasta
    fi

    if ${DO_FILTER_MMSEQS_WITH_FULL_LENGTH_SEQS} && ${DO_RUN_WITH_MMSEQS}
    then
        # Ideally we should run this set on the common species once everything has been filtered out.
        # We save the intermediary files in the msa directory
        rm sequence_${i}_filtered_FL.fasta
        ${PYTHON} ${ALPHASCRIPTS}/toolsFilterMmseqsOut.py -i alisequence_${i}.tsv -a sequence_${i}_FL.fasta -f sequence_${i}_filtered_FL.fasta -c ${COV} -q ${QID}
        cp sequence_${i}_filtered_FL.fasta msa_${ROOTNAME}/${ROOTNAME}_filtered_FL.fasta
    else
        # Since we ask no MSA building, it means it already exists
        rm sequence_${i}_filtered_FL.fasta
        ln -s msa_${ROOTNAME}/${ROOTNAME}_filtered_FL.fasta sequence_${i}_filtered_FL.fasta
    fi

    if ${DO_RUN_MAFFT}
       then
         # 4a - Keep only the first best hits not to overload mafft
         #head -${TOPLINES} sequence_${i}_filtered_FL.fasta > sequence_${i}_filtered_FL_top${MAX_NSEQ_IN_MAFFT}.fasta
         rm sequence_${i}_filtered_FL_top${MAX_NSEQ_IN_MAFFT}.fasta
         rm sequence_${i}_filtered_FL_mafali.fasta
         head -${TOPLINES} sequence_${i}_filtered_FL.fasta > sequence_${i}_filtered_FL_top${MAX_NSEQ_IN_MAFFT}.fasta

         # 4 - Run mafft
         #mafft --randomseed 1 --thread -1 sequence_${i}_filtered_FL_top${MAX_NSEQ_IN_MAFFT}.fasta > sequence_${i}_filtered_FL_mafali.fasta

         if ${DO_RUN_WITH_MMSEQS}
            then
            cp sequence_${i}_filtered_FL_top${MAX_NSEQ_IN_MAFFT}.fasta msa_${ROOTNAME}/${ROOTNAME}_filtered_FL_top${MAX_NSEQ_IN_MAFFT}.fasta
            cp sequence_${i}_filtered_FL_mafali.fasta msa_${ROOTNAME}/${ROOTNAME}_filtered_FL_mafali.fasta
         fi

    elif ${DO_RUN_WITH_MMSEQS}
       then
         # We assume mafft has already been run before. We have to properly link the correct files with respect to the root_name inputs
         rm sequence_${i}_filtered_FL_top${MAX_NSEQ_IN_MAFFT}.fasta
         rm sequence_${i}_filtered_FL_mafali.fasta
         ln -s msa_${ROOTNAME}/${ROOTNAME}_filtered_FL_top${MAX_NSEQ_IN_MAFFT}.fasta sequence_${i}_filtered_FL_top${MAX_NSEQ_IN_MAFFT}.fasta
         ln -s msa_${ROOTNAME}/${ROOTNAME}_filtered_FL_mafali.fasta sequence_${i}_filtered_FL_mafali.fasta

    fi

done


if ${DO_MATCH_SPECIES}
then

#######################################################################################################################################################################
#MERGE ALL PROTEIN
#######################################################################################################################################################################

# 5 - Match species
${PYTHON} ${ALPHASCRIPTS}/CommonAndGapSpecies.py -nb ${#index_partners[@]} --server -json json_step1 -ncsname numbercommonspecies ${FULL_LENGTH_OPTION} ${PARALOG_FREE_OPTION}
# ${#index_partners[@]} expands to 2
# input: sequence_0_filtered_FL_mafali.fasta, sequence_1_filtered_FL_mafali.fasta
# input without -FL: sequence_0_filtered.fasta, sequence_1_filtered.fasta
# optional parameter: -pf (paralogfree). Only one instance per species of every partner
# output stdout: a count of the common species
# output json_step1.json: [for msa0: {"inputseq":msa, "paired": msa, "unpaired":msa}, for msa1: {"inputseq":msa, "paired": msa, "unpaired":msa}]
fi


if ${DO_PAIR_SEQ_DELIM}
then

	#######################################################################################
			#USER INPUTS : PARAMETERS OF GetAlphaFoldInput 
	#######################################################################################
  # 6 - write Alphafold inputs
  if ${RUN_MIXED_ALI_OPTION}
  then
  ${PYTHON} ${ALPHASCRIPTS}/GetAlphaFoldInputs.py -msa 0,1 --delimitation ${delim_list[0]},${delim_list[1]} --jsonname json_step1 --jsonoutput json_step2 --colab -st ${STOICHIOMETRY} -o mixed
  fi

  if ${RUN_PAIRED_ALI_OPTION}
  then
  ${PYTHON} ${ALPHASCRIPTS}/GetAlphaFoldInputs.py -msa 0,1 --delimitation ${delim_list[0]},${delim_list[1]} --jsonname json_step1 --jsonoutput json_step2 --colab -st ${STOICHIOMETRY} -o paired
  fi

  if ${RUN_ALLUNPAIRED_ALI_OPTION}
  then
  ${PYTHON} ${ALPHASCRIPTS}/GetAlphaFoldInputs.py -msa 0,1 --delimitation ${delim_list[0]},${delim_list[1]} --jsonname json_step1 --jsonoutput json_step2 --colab -st ${STOICHIOMETRY} -o all_unpaired
  fi

  #output Alphafold_input_Sequence, MSAforAlphafold_mixed.fasta, MSAforAlphafold_paired.fasta, MSAforAlphafold_all_unpaired.fasta
  #output output/json_step2.json: [for msa0: {"inputseq":msa cropped, "paired": msa cropped, "unpaired":msa cropped}, for msa1: {"inputseq":msa cropped, "paired": msa cropped, "unpaired":msa cropped}]

fi
