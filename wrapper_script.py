#Import Necessary Dependencies
import os
import subprocess
#Create a directory for output of all files
os.system("mkdir PipelineProject_John_Floros")
os.chdir("PipelineProject_John_Floros")
os.system("echo > PipelineProject.log")

HCMV_cds = os.system("datasets download virus  genome accession NC_006273.2 --include cds")
os.system("unzip ncbi_dataset.zip")
os.chdir("ncbi_dataset/data")
genomic_CDS_num_command = subprocess.run("grep -c '>' cds.fna", shell=True, capture_output=True, text=True)
genomic_CDS_num = genomic_CDS_num_command.stdout.strip()
print(genomic_CDS_num)
#reference_transcriptome = "kallisto index -1 index.idx " +
