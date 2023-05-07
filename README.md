# Running the scripts in SCAN_IDR

SCAN_IDR is a pipeline developed to evaluate AlphaFold2 success rates in predicting the structures of protein assemblies interacting through an intrinsically disordered region. The pipeline evaluates how the delimitations of the inputs and the content of the multiple sequence alignments affect the reliability of the generated models. The study was evaluated on 42 complexes that have no sequence or structural similarity to the complexes used in the AlphaFold2 training.  

The pipeline goes through 10 steps including:<br>
1. the retrieval of the pdb for the cases to be tested (*[Steps 1](#step1) & [2](#step2)*).
2. the generation of AlphaFold2 runs for all these cases using [**10 different protocols**](#protocols) which sample different alignment and delimitation conditions (*[Steps 3](#step3) & [4](#step4)*).
3. the evaluation of the generated models (*[Steps 5](#step5) to [9](#step9)*).
4. the visualisation of the results using Jupyter Notebooks provided in `scripts/post_process/` directory (*[Step 10](#step10)*).

At each step, different folders and command files are generated which can be run by the user [as described below](#steps). A user can also set the way commands are to be executed in their own environment if necessary (parallel work, GPU usage) by adapting these command files. The use of servers equipped with GPU cards is strongly recommended to obtain best performance from AlphaFold2. <br>
### Content of the archive 
- `scripts/`: contains the scripts used to generate the dataset, run the predictions and analyse the results <br />
- `data/`: contains inputs processed by the python scripts to run the full dataset<br>
- `data_demo/`: contains inputs and outputs obtained for two examples from the dataset<br>
- `scan_idr.yml`: defines the environment to run the scan_idr pipeline<br>
- `colabfold_v1.3.0.def`: a singularity definition file to install the required version of [ColabFold](https://github.com/sokrypton/ColabFold) and AlphaFold2 <br>
- `install_uniref30_2202_db.sh`: a script to install the [ColabFold sequence database](https://colabfold.mmseqs.com/)  <br>

Packages required to run the `scan_idr` pipeline are listed in the conda environment file `scan_idr.yml` (install time ~ 5 minutes). Conda environment can be activated with the instruction: <br>

    conda env create -n scan_idr --file scan_idr.yml

`<WORKING_DIR>` can be defined in `scripts/config.ini`. Default `<WORKING_DIR> = ../` <br>

The list of pdb codes used as input is : `<WORKING_DIR>/data/list_cif.txt`, comma separated list of pdb codes to retrieve the reference files from the  rcsb<br>
<br>

- The pipeline can be run from [**Steps 1 to 10**](#steps) on all the 42 entries using the input files provided in the `data/` folder.<br>

- The full dataset of sequences, alignments, reference structures and structural models and generated for the 42 cases using the scripts of the `scan_idr` pipeline can be 
retrieved from the `scanidr_data_repository.tar` archive available at (https://zenodo.org/record/7838024). <br>
- A `data_demo/` folder containing inputs for 2 test cases with their expected outputs is [presented below](#demo_section)<br>

### Prerequisites

Most of the dependencies are installed upon the creation of the `scan_idr` conda environment explained above.<br>

Independent installation of [MMseqs2](https://github.com/soedinglab/MMseqs2), [ColabFold](https://github.com/sokrypton/ColabFold),  and [ProFit](http://www.bioinf.org.uk/software/profit/) softwares are required to run the full pipeline (install time ~ 20 minutes).<br>

- ***[MMseqs2](https://github.com/soedinglab/MMseqs2)***: This software is used in [Step 3](#step3) to generate the initial multiple sequence alignments. A suitably advanced version of MMseqs2 must be obtained from the github repository.<br>
  - Retrieve the program code with the commands below. Instructions taken from [compile-from-source-under-linux](https://github.com/soedinglab/MMseqs2/wiki#compile-from-source-under-linux) <br>
  - Modify the line `<your_path>/MMseqs2/bin/mmseqs` with the correct path in the `SCAN_IDR/scripts/config.ini` file<br>


    git clone https://github.com/soedinglab/MMseqs2.git
    cd MMseqs2/
    git checkout tags/14-7e284
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
    make
    make install



- ***[Colabfold](https://github.com/sokrypton/ColabFold)***: The `scan_idr` pipeline uses a colabfold singularity image to run AlphaFold2 in [Step 4](#step4). <br>
  - The singularity software v3.8.6 is installed at the creation of the `scan_idr` conda environment. <br>
  - A singularity definition file is provided in `SCAN_IDR/colabfold_v1_3_0.def` to build the `colabfold_container_v13.sif` image file (adapted from [Tubiana](https://github.com/tubiana/colabfold_singularity)). The generated singularity file `colabfold_container_v13.sif` corresponds to an image of the colabfold version 1.3 compatible with the `scan_idr` pipeline.<br>
  - To create the colabfold singularity image, root privilege is needed and the command below may more easily be performed on a local linux machine even without GPUs.<br>
  - Once built, the file `colabfold_container_v13.sif` can be copied to the machine running the `scan_idr` pipeline such as a server or a cluster with GPUs. <br>
  - Modify the line `<your_path>/colabfold_container_v13.sif` with the correct path in the `SCAN_IDR/scripts/config.ini` file<br>
  - Build the `.sif` image with :


    sudo singularity build colabfold_container_v13.sif colabfold_v1_3_0.def

*NB: Other implementations of colabfold may be used, provided the command line used in the script Step 4 is adapted.* <br>

- ***[ProFit](http://www.bioinf.org.uk/software/profit/)***: This software is used to calculate different RMSD metrics between the models and the reference structures in [Step 6](#step6). The `scan_idr` was initially developed to run with ProFit software version 3.1 and was tested with the most recent version 3.3.<br>
The ProFit software is freely available at http://www.bioinf.org.uk/software/profit/: <br>

    - Download the **ProFit V3.3 Linux Binary Distribution** and unpack the tar file.<br>
    - The ProFit executable will be accessible in `<your_path>/ProFitV3.3/bin/profit`
    - Modify the line `<your_path>/ProFitV3.3/bin/profit` with the correct path in the `SCAN_IDR/scripts/config.ini` file<br>
    - Create the environment variables HELPDIR and DATADIR which should both point to the top ProFit directory where the
       files ProFit.help and mdm78.mat are stored: <br>
                               

    export HELPDIR=<your_path>/ProFitV3.3
    export DATADIR=<your_path>/ProFitV3.3

- ***[uniref30_2202_db database](https://colabfold.mmseqs.com/)***: Installation of this ColabFoldDB is  required to build the alignments in [Step 3](#step3).<br>
  - The script `install_uniref30_2202_db.sh` is provided in `SCAN_IDR/`. <br>
  - This script is directly adapted from the ColabFold script [`setup_databases.sh`](https://github.com/sokrypton/ColabFold/blob/main/setup_databases.sh) <br>
  - Modify the line `<your_path>/uniref30_2202/uniref30_2202_db` with the correct path in the `SCAN_IDR/scripts/config.ini` file. <br>
  - Searches against this ColabFoldDB requires a machine with ~ 128GB RAM.<br>
  - Run the following command to download and install the database: <br>

          sh install_uniref30_2202_db.sh





## How to run the pipeline:

## Step 1: Retrieve mmcif files from rscb <a id="steps"></a> <a id="step1"></a>
    
    cd SCAN_IDR/scripts
    ./1_batch_download.sh -f <WORKING_DIR>/data/list_cif.txt -o <WORKING_DIR>/data/cif -c

- Output: PDB files will be stored in the directory `<WORKING_DIR>/data/cif/` with the mmcif format  

## Step 2: Parse mmcif files and get uniprot ID and delimitations <a id="step2"></a>

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
 
## Step 3: Generate the MSAs for every unicode entry <a id="step3"></a>


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

## Step 4: Create the concatenated MSA for every protocol and generate the models <a id="step4"></a>

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

## Step 5: After AF2 has finished, cut of the models to the same delimitations as in the reference PDB to enable CAPRI-like evaluation <a id="step5"></a>

Command to run in `<WORKING_DIR>/scripts`: 

    python 5_CutModelsForCAPRI.py 

- __Ouput:__
    - Creation of a directory `<WORKING_DIR>/data/cutmodels_for_caprieval/`
    - Creation of a directory `<WORKING_DIR>/data/cutmodels_for_caprieval/<CASE_INDEX>_<PDBCODE>/`
    - For every protocol generation of an af2 directory containing the 25 af2 models in which receptors and ligands are cut following the delimitation of the reference PDB<br>

## Step 6: Run the evaluation of the models using CAPRI criteria <a id="step6"></a>

Command to run in `<WORKING_DIR>/scripts`: 

    python 6_RunCapri.py 

- __Ouput:__
  - Creation of a directory `<WORKING_DIR>/data/caprieval_on_cutmodels/`
  - Creation of a directory `<WORKING_DIR>/data/caprieval_on_cutmodels/<CASE_INDEX>_<PDBCODE>/`
  - Creation of the script `cmd_run_capri_all_models.sh` in `<WORKING_DIR>/data`<br>

Move to `<WORKING_DIR>/data/` to run :

    sh cmd_run_capri_all_models.sh

## Step 7: Post-processing - Retrieval of the AF2 scores <a id="step7"></a>

Command to run in `<WORKING_DIR>/scripts`: 

    python 7_GetAF2Scores.py 

- __Ouput:__
  - Creation of a directory `<WORKING_DIR>/data/results_files/`
  - Generation of a file `AF2_AllProtocols_Scores.out` listing all the AF2 scores of all the models in `<WORKING_DIR>/data/results_files/`<br>


## Step 8: Post-processing - Retrieval of the CAPRI evaluation scores <a id="step8"></a>

Command to run in `<WORKING_DIR>/scripts`: 

    python 8_GenerateGlobalOutputTable.py 

- __Ouput:__
  - Generation of a file listing all the capri evaluations and the AF2 scores for all the models `Global_results_AF2_CAPRI.out` in `<WORKING_DIR>/data/results_files/`<br>

## Step 9: Post-processing - Retrieval of the CAPRI evaluation scores <a id="step9"></a>

Command to run in `<WORKING_DIR>/scripts`: 

    python 9_RunFinalScoresAnalyses.py

- __Ouput:__
    - For every protocols, generate several files comparing the number of successful predictions over all the `Global_results_AF2_CAPRI.out` in `<WORKING_DIR>/data/results_files/`<br>
*These files will be used as inputs for the graphical analyses performed in the Jupyter Notebooks scripts located in `<WORKING_DIR>/scripts/post_process` directory.*

## Step 10: Plotting the results <a id="step10"></a>

Several Jupyter notebooks are provided to plot the results of the analysis and are provided in: `<WORKING_DIR>/scripts/post_process` <br>

Move to `<WORKING_DIR>/scripts/post_process` to run each of these notebooks: <br>

1. File `Report_SuccessRates_InvidualMethods.ipynb`<br>
   * *Plots as stacked bars the success rates of different protocols*<br> 
2. File `Report_DockQVsPDB.ipynb`<br>
   * *Plots as histogram for every pdb, the best DockQ scores of the model*<br>
3. File `Report_ScoresVsCAPRI.ipynb`<br>
   * *Plots as box plots the distribution of the AF2 scores vs DockQ scores*<br>
4. File `Report_DatasetProperties_BoxPlots_lengths.ipynb`<br>
   * *Plots as box plots the length distributions of the systems studied*<br> 

*NB: Depending if you changed the definition of your `WORKING_DIR` in the `config.ini` file, the path defining WORKING_DIR variable in the Jupyter notebooks may need to be adjusted* 

## Running a <a id="demo_section">demo example</a> 

To test the pipeline, it is possible to work on one or two pdb codes rather than running all 42 entries of the dataset. <br><br> 
To do so, edit the `SCAN_IDR/data/list_cif.txt` and change the comma separated list of PDB codes. <br>

We provide a `SCAN_IDR/data_demo/` folder as an example of a simple input with two PDB entries with all the outputs that are expected to be generated from Step 1 to 10: <br>
- First, uncompress the large directories in `SCAN_IDR/data_demo/`: `fasta_msa.tar.gz`, `af2_runs.tar.gz` and `cutmodels_for_caprieval.tar.gz`  
- The file `SCAN_IDR/data_demo/list_cif.txt` can be copied to your `<WORKING_DIR>/data/` and used as input of [**Steps 1 to 10**](#steps) to run only 2 cases. <br> 
- If you want to run the pipeline starting at a specific Step, the entire content of `SCAN_IDR/data_demo/` can be copied to `<WORKING_DIR>/data/`. Thanks to that, all the files and folders required to run every Step independently will be accessible. For instance, it is then possible to start running the pipeline at Step 5 and all next steps subsequently. <br>
  *NB: Only to run the last step Step 10, you need to have at least run Step 9 on your system.*<br> 

The files and folders from the `SCAN_IDR/data_demo/` can be compared to the outputs of the pipeline. <br>
- Indicative times for the execution of the demo: <br>
  - Step 1: < 1 min <br>
  - Step 2: < 1 min <br>
  - Step 3: 5 min (including python script and MSA generation on 80 CPU) <br>
  - Step 4: 1-2 min for the python script + 5 hours to run AlphaFold on a single A100 GPU (can be bypassed by using data_demo af2_runs files)  <br>
  - Step 5: 10 min <br>
  - Step 6: <1 min for the python script + 5 min for running CAPRI evaluation <br>
  - Step 7: <1 min <br>
  - Step 8: <1 min <br>
  - Step 9: <1 min <br>
  - Step 10: a few minutes to run each Jupyter notebook <br>


## Acknowledgements

- We would like to thank the [ColabFold](https://github.com/sokrypton/ColabFold) and [AlphaFold](https://github.com/deepmind/alphafold) teams for providing open access to their softwares.<br>
- We are grateful to [Martin, A.C.R.](http://www.bioinf.org.uk/software/profit/) for the development of the ProFit software based on the McLachlan algorithm *(McLachlan, A.D., 1982, Acta Cryst A38, 871-873)*.<br>
- Thanks to [Thibault Tubiana](https://github.com/tubiana/colabfold_singularity) and Chloé Quignot for providing their singularity definition file for ColabFold.<br>
