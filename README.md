# Running the scripts in SCAN_IDR

The pipeline goes through 10 steps including:<br>
1. the retrieval of the pdb for the cases to be tested.
2. the generation of AlphaFold2 runs for all these cases using [**10 different protocols**](#protocols) which sample different alignment and delimitation conditions.
3. the evaluation of the generated models.
4. the visualisation of the results using Jupyter Notebooks provided in `scripts/post_process/` directory.

At each step, different folders and files are generated and the user can tune how the commands should be run in its own environment (parallel jobs, GPU usage).<br>
<br>
### Content of the archive 
- `scripts/` contains the scripts used to generate the dataset, run the predictions and analyse the results <br />
- `data/` contains inputs processed by the python scripts to run the dataset<br>
- `scan_idr.yml` defines the environment to run the scan_idr pipeline<br> 

Packages required to run the `scan_idr` pipeline are listed in the conda environment file `scan_idr.yml`. Conda environment can be activated with the instruction: <br>

    conda env create -n scan_idr --file scan_idr.yml

`<WORKING_DIR>` can be defined in `scripts/config.ini`. Default `<WORKING_DIR> = ../` <br>

The list of pdb codes used as input is : `<WORKING_DIR>/data/list_cif.txt`, comma separated list of pdb codes to retrieve the reference files from the  rcsb<br>
<br>

## Step 1 : Retrieve mmcif files from rscb
    
    cd SCAN_IDR/scripts
    ./1_batch_download.sh -f <WORKING_DIR>/data/list_cif.txt -o <WORKING_DIR>/data/cif -c

- Output: PDB files will be stored in the directory `<WORKING_DIR>/data/cif/` with the mmcif format  

## Step 2 : Parse mmcif files and get uniprot ID and delimitations

    python 2_GetUnicodeDelim_mmcif.py


- Output: File listing the uniprot ID and delimitations of the chains : `File_Listing_Uniprot_Inputs.txt`

The uniprot and delimitations are automatically extracted from the `cif`file. <br>
In case the information is not found, the value is replaced by `???`.  
File can be edited and corrected manually following the format :

    # PDBCODE: 5NCL
    CHAIN:A	UNIPROT:P53894	START:251	STOP:756
    CHAIN:B	UNIPROT:P43563	START:46	STOP:287
    CHAIN:D	UNIPROT:P24276	START:205	STOP:214

 Edit and correct the file `File_Listing_Uniprot_Inputs.txt` if required.
 
## Step 3 : Generate the MSAs for every unicode entry


Command to run: <br>

    python 3_CmdRetrieveFastaAliFile.py

Use the information in `File_Listing_Uniprot_Inputs.txt`<br>
- Retrieve the fasta file for every entry from uniprot db and store it in `<WORKING_DIR>/data/fasta_msa`
- Create a command file `<WORKING_DIR>/data/fasta_msa/cmd_create_alignments.sh` that will : <br>
    - Run mmseqs on every fasta input
    - Retrieve __full-length sequence__ of every homolog from the resulting MSA
    - Realign __full-length sequences__ using mafft

Move to `<WORKING_DIR>/data/fasta_msa/` where the `cmd_create_alignments.sh` script was created. <br>
This script can be edited if required to modify the generic values of parameters such as QID, COV for specific entries <br>
Then, run:

    sh cmd_create_alignments.sh

For every uniprot entry, a directory `msa_<uniprot ID>` will be created in `<WORKING_DIR>/data/fasta_msa/`<br>
These directories contain the multiple sequence alignments that will be used to generate the concatenated alignments in subsequent steps.<br>
A copy of the MSA files generated in these directories or the 42 tested cases is provided as an `.a3m` in the final repository distributed as *__Available Data__*  

## Step 4 : Create the concatenated MSA for every protocol

Command to run in `<WORKING_DIR>/scripts`: 

    python 4_PrepareDirectoryforAF2.py -i all 

With the option `-i all`, the input files for ten protocols wil be generated. 

The list of <a id="protocols">10 protocols</a> to be run is provided below: <br>
- **mixed_ali-delim-delim**: paired+unpaired alignment delimited as in the pdb for both receptor and ligand.
- **mixed_ali-fl-fl**: paired+unpaired alignment using the full-length sequences of both receptor and ligand.
- **mixed_ali-delim-fl**: paired+unpaired alignment delimited as in the pdb for the receptor but using the full-length sequence of the ligand.
- **mixed_ali-delim-100**: paired+unpaired alignment delimited as in the pdb for the receptor and extending the size of the PDB ligand sequence by 100 residues.
- **mixed_ali-delim-200**: paired+unpaired alignment delimited as in the pdb for the receptor and extending the size of the PDB ligand sequence by 200 residues.
- **unpaired_ali-delim-delim**: unpaired alignment using the full-length sequences of both receptor and ligand.
- **unpaired_ali-delim-fl**: unpaired alignment delimited as in the pdb for the receptor but using the full-length sequence of the ligand.
- **single_pep-delim-delim**: alignment for the receptor and only single sequence for the ligand using the delimitations as in the pdb for both receptor and ligand.
- **single_pep-delim-100**: alignment for the receptor and only single sequence for the ligand using the delimitations as in the pdb for the receptor and extending the size of the PDB ligand sequence by 100 residues.
- **single_pep-delim-200**: alignment for the receptor and only single sequence for the ligand using the delimitations as in the pdb for the receptor and extending the size of the PDB ligand sequence by 200 residues.
  

- __Ouput:__
    - Creation of a directory `<WORKING_DIR>/data/af2_runs/`
    - Creation of a directory `<WORKING_DIR>/data/af2_runs/<CASE_INDEX>_<PDBCODE>/`
    - For every protocol generation of an af2 directory containing the concatenated MSA and the `run_colab.sh` script to run colabfold on this MSA
    - Generation of the output file `<WORKING_DIR>/data/af2_runs/cmd_global_af2runs.sh` which lists all the commands to run the protocols for every entry<br>
      
Move to `<WORKING_DIR>/data/af2_runs/` to run :

    sh cmd_global_af2runs.sh

*NB: If needed, each of these protocols can be generated individually using : `python 4_PrepareDirectoryforAF2.py -i <protocol_name>`*<br>

- Explanations of the format of the protocol names: `<MSA_mode>-<Receptor_Delimitations>-<Ligand_Delimitations>`<br>
`<MSA_mode>`: can be either `mixed` (paired+unpaired), `unpaired` or `single_pep` (no MSA, only single sequence)<br>
`<Delimitations>`: can be `delim` (same delimitations as in the reference PDB), `fl` (full-length sequence), `100` or `200` (delimitations of the reference PDB extended by this number)

## Step 5 : After AF2 has finished, cut of the models to the same delimitations as in the reference PDB to enable CAPRI-like evaluation

Command to run in `<WORKING_DIR>/scripts`: 

    python 5_CutModelsForCAPRI.py 

- __Ouput:__
    - Creation of a directory `<WORKING_DIR>/data/cutmodels_for_caprieval/`
    - Creation of a directory `<WORKING_DIR>/data/cutmodels_for_caprieval/<CASE_INDEX>_<PDBCODE>/`
    - For every protocol generation of an af2 directory containing the 25 af2 models in which receptors and ligands are cut following the delimitation of the reference PDB<br>

## Step 6 : Run the evaluation of the models using CAPRI criteria

Command to run in `<WORKING_DIR>/scripts`: 

    python 6_RunCapri.py 

- __Ouput:__
  - Creation of a directory `<WORKING_DIR>/data/caprieval_on_cutmodels/`
  - Creation of a directory `<WORKING_DIR>/data/caprieval_on_cutmodels/<CASE_INDEX>_<PDBCODE>/`
  - Creation of the script `cmd_run_capri_all_models.sh` in `<WORKING_DIR>/data`<br>

Move to `<WORKING_DIR>/data/` to run :

    sh cmd_run_capri_all_models.sh

## Step 7 : Post-processing - Retrieval of the AF2 scores

Command to run in `<WORKING_DIR>/scripts`: 

    python 7_GetAF2Scores.py 

- __Ouput:__
  - Creation of a directory `<WORKING_DIR>/data/results_files/`
  - Generation of a file `AF2_AllProtocols_Scores.out` listing all the AF2 scores of all the models in `<WORKING_DIR>/data/results_files/`<br>


## Step 8 : Post-processing - Retrieval of the CAPRI evaluation scores

Command to run in `<WORKING_DIR>/scripts`: 

    python 8_GenerateGlobalOutputTable.py 

- __Ouput:__
  - Generation of a file listing all the capri evaluations and the AF2 scores for all the models `Global_results_AF2_CAPRI.out` in `<WORKING_DIR>/data/results_files/`<br>

## Step 9 : Post-processing - Retrieval of the CAPRI evaluation scores

Command to run in `<WORKING_DIR>/scripts`: 

    python 9_RunFinalScoresAnalyses.py

- __Ouput:__
    - For every protocols, generate several files comparing the number of successful predictions over all the `Global_results_AF2_CAPRI.out` in `<WORKING_DIR>/data/results_files/`<br>
*These files will be used as inputs for the graphical analyses performed in the Jupyter Notebooks scripts located in `<WORKING_DIR>/scripts/post_process` directory.*

## Step 10 : Plotting the results

Several Jupyter notebooks are provided to plot the results of the analysis and are provided in: `<WORKING_DIR>/scripts/post_process` <br>

Move to `<WORKING_DIR>/scripts/post_process` to run each of these notebooks: <br>

1. File `Report_SuccessRates_InvidualMethods.ipynb`<br>
   * *Plots as stacked bars the success rates of different protocols*<br> 
2. File `Report_DockQVsPDB.ipynb`<br>
   * *Plots as histogram for every pdb, the best DockQ scores of the model*<br>
3. File `Report_ScoresVsCAPRI.ipynb`<br>
   * *Plots as box plots the distribution of the AF2 scores vs DockQ scores*<br<>>
4. File `Report_DatasetProperties_BoxPlots_lengths.ipynb`<br>
   * *Plots as box plots the length distributions of the systems studied*<br> 
