#!/bin/sh

# Inspired from https://github.com/deepmind/alphafold/issues/10
#
# singularity build adapted from https://github.com/tubiana/colabfold_singularity/blob/main/colabfold_receipe.def
#
########## CHANGE THE OPTIONS BELOW THIS LINE
COLAB_VERSION=/opt/alphafold/colabfold_container_v13.sif

FASTA_FILE=$1 # change the name of the a3m file paired/mixed or other
export SINGULARITYENV_CUDA_VISIBLE_DEVICES=0 # Set to 0 or 1 to select a specific GPU card or remove the line
PARAMETERS='AlphaFold2-multimer-v2' # use PARAMETERS='AlphaFold2-ptm' for running the monomer-based parameters
NUMBER_OF_MODELS=5
NUMBER_OF_RECYCLES=3
NUMBER_OF_REPEATS=5 # changing the random seed

######### IN GENERAL, BELOW THIS LINE, DO NOT MODIFY UNLESS EXPERT USAGE

FASTA_DIR=$(pwd) # If the current directory is not the working one, define here the path to the working directory

OUTPUT_DIR=${FASTA_DIR}/af2_output
MSA_DIR=${FASTA_DIR}/af2_msas
PRED_DIR=${FASTA_DIR}/af2_predictions

mkdir -p ${FASTA_DIR} &> /dev/null
mkdir -p ${OUTPUT_DIR} &> /dev/null
mkdir -p ${MSA_DIR} &> /dev/null
mkdir -p ${PRED_DIR} &> /dev/null

cp ${FASTA_FILE} ${MSA_DIR} # comment if mmseqs search is run below


for (( c=1; c<=${NUMBER_OF_REPEATS}; c++ ))
do
RANDOMSEED=$RANDOM
echo "Running colabfold_batch with random seed : ${RANDOMSEED}" >> ${PRED_DIR}/log.txt
singularity exec --env TF_FORCE_UNIFIED_MEMORY=1,XLA_PYTHON_CLIENT_MEM_FRACTION=4.0,OPENMM_CPU_THREADS=8 \
            -B ${OUTPUT_DIR}:/inout/output -B ${FASTA_DIR}:/inout/fasta -B ${MSA_DIR}:/inout/msas -B ${PRED_DIR}:/inout/predictions \
            --nv ${COLAB_VERSION} \
            colabfold_batch \
            --model-type ${PARAMETERS} \
            --num-recycle ${NUMBER_OF_RECYCLES} \
            --num-models ${NUMBER_OF_MODELS} \
            --random-seed ${RANDOMSEED} \
            /inout/msas/${FASTA_FILE} /inout/predictions
	mv af2_predictions af2_predictions_v${c}
	cd af2_predictions_v${c}
	for i in *.pdb
	do
		mv -v "${i}" "${i%.*}_v${c}.${i##*.}"
	done
	cd ..
	mkdir -p ${PRED_DIR} &> /dev/null
done



