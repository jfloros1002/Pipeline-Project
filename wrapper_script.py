#Import Necessary Dependencies
import os
#Create a directory for output of all files
os.system("mkdir PipelineProject_John_Floros")
os.chdir("PipelineProject_John_Floros")
HCMV_cds = os.system("datasets download virus  genome accession NC_006273.2 --include cds")

#reference_transcriptome = "kallisto index -1 index.idx " +
