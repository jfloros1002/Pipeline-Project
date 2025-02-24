#Import Necessary Dependencies
import os
import subprocess
#Create a directory for output of all files
os.system("mkdir PipelineProject_John_Floros")
os.system("cp *.fastq PipelineProject_John_Floros") 
os.chdir("PipelineProject_John_Floros")

#Create an file for the output of each command in the pipeline
os.system("echo > PipelineProject.log")

#Create  a genome for HCMV using NCBI command datasets
HCMV_cds = os.system("datasets download virus genome accession NC_006273.2 --include cds")
#Unzip the file and go inside to analyze
os.system("unzip ncbi_dataset.zip")
os.system("cp ncbi_dataset/data/cds.fna .")
#Find total number of CDS for first output to .log file using subprocess module, output if captured in genomic_CDS_num
genomic_CDS_num_command = subprocess.run("grep -c '>' cds.fna", shell=True, capture_output=True, text=True)
genomic_CDS_num = genomic_CDS_num_command.stdout.strip()
#Write first output of .log file
with open("PipelineProject.log","a+") as f:
	f.write("The HCMV genome (NC_006273.2) has " +str( genomic_CDS_num) + " CDS.\n")

#Building index from HCMV
reference_transcriptome = "kallisto index -i index.idx cds.fna"
os.system(reference_transcriptome)

quantification_one_2dpi = "kallisto quant -i index.idx -o quantification_results_one_2dpi -b 10 -t 2 SRR5660033_1.fastq SRR5660033_2.fastq"
os.system(quantification_one_2dpi)
os.chdir("quantification_results_one_2dpi")
with open("abundance.tsv","r")as f:
	lines = f.readlines()[1:]
	base_results = [(i.split("\t")[4]) for i in lines]
	results = [float(x.strip()) for x in base_results]
def find_median(results):
    results.sort()
    list_length = len(results)
    if list_length % 2 == 0:
        mid1 = results[list_length // 2 - 1]
        mid2 = results[list_length // 2]
        median = (mid1 + mid2) / 2
    else:
        median = results[list_length // 2]
    return median
min_1_2dpi = min(results)
med_1_2dpi = find_median(results) 
sum_results = sum(results)
length_results = len(results)
print(sum_results)
print(length_results)
mean_1_2dpi = sum_results/length_results
max_1_2dpi = max(results)
os.chdir("..")
with open("PipelineProject.log", "a+") as f:
    f.write("sample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n")
    f.write("Donor 1\t" + "2dpi\t" + str(min_1_2dpi) + "\t" + str(med_1_2dpi) + "\t" + str(mean_1_2dpi) + "\t" + str(max_1_2dpi) + "\n")
