#Import Necessary Dependencies
import os
import subprocess

#Create a directory for output of all files
os.system("mkdir PipelineProject_John_Floros")
os.chdir("PipelineProject_John_Floros")

#Create an file for the output of each command in the pipeline
os.system("echo > PipelineProject.log")

#Create  a genome for HCMV using NCBI command datasets
HCMV_cds = os.system("datasets download virus genome accession NC_006273.2 --include cds")
#Unzip the file and go inside to analyze
os.system("unzip ncbi_dataset.zip")
os.chdir("ncbi_dataset/data")
#Find total number of CDS for first output to .log file using subprocess module, output if captured in genomic_CDS_num
genomic_CDS_num_command = subprocess.run("grep -c '>' cds.fna", shell=True, capture_output=True, text=True)
genomic_CDS_num = genomic_CDS_num_command.stdout.strip()
#Write first output of .log file
os.system("mv cds.fna ~/Comp_Bio/Pipeline-Project/PipelineProject_John_Floros")
os.chdir("/home/2025/jfloros/Comp_Bio/Pipeline-Project/PipelineProject_John_Floros")
print(genomic_CDS_num)
with open("PipelineProject.log","a+") as f:
	f.write("The HCMV genome (NC_006273.2) has " +str( genomic_CDS_num) + " CDS.")
reference_transcriptome = "kallisto index -i index.idx cds.fna"
os.system(reference_transcriptome)
