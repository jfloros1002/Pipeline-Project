# Pipeline-Project
Repository for COMP483 Pipeline Project

# Required Dependences:
os

suprocess

NCBI Datasets (NCBI Command Line Tools)

kallisto

sleuth

dplyr

bowtie2

SPAdes

BLAST+

# Aquisition of Data
Data was aquired from individuals two and six days post infection (dpi) with Human Cytomegalovirus (HCMV). Raw files were found using SRA acession numbers and downloaded using wget with the following links:

Donor 1 2dpi: https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030

Donor 2 2dip: https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044

Donor 1 6dpi: https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033

Donor 2 6dip: https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045

Files were placed in repository under acession numbers:

SRR5660030

SRR5660044

SRR5660033

SRR5660045

Each file is a binary file that can then be turned into paired end FASTQ reads using the fasterq-dump command on the command line:

fasterq-dump SRR5660030

fasterq-dump SRR5660033

fasterq-dump SRR5660044

fasterq-dump SRR5660045

These paired end flies are then used for the analysis using the wrapper_script.py command

NOTE: It is important that the files names contain the SRR56600XX.fastq file name, as this is used in the pipeline. Samples files will contain the same title.

# Running the Pipeline

To run the pipeline: Eight sample sequences must be present in the same directory as the wrapper script. These can be either the sample files present in the github repository, or the original files downloaded as described above using wget and fasterq_dump. Additionally the sleuth_command.R file must also be downloaded, as this is used in the wrapper script.

Once the eight files, the sleuth_command.R, and wrapper_script.py are downloaded, analysis can be completed by writing python wrapper_script.py to the command line. A directory named PipelineProject_John_Floros will be created, and all output will be placed in this folder. If the pipeline has already been ran and additional runs are producing error, it is reccomended to delete the entire output directory using rm -r PipelineProject_John_Floros and run the wrapper script again with a fresh environment.
