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
