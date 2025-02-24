#Import Necessary Dependencies
import os
import subprocess
#Create a directory for output of all files
os.system("mkdir PipelineProject_John_Floros")
os.system("cp *.fastq PipelineProject_John_Floros") 
os.system("cp sleuth_command.R PipelineProject_John_Floros")
os.chdir("PipelineProject_John_Floros")

#Create an file for the output of each command in the pipeline
os.system("echo > PipelineProject.log")

#Create  a genome for HCMV using NCBI command datasets
#HCMV_cds = os.system("datasets download virus genome accession NC_006273.2 --include cds")
HCMV_cds = os.system("datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report")
#Unzip the file and go inside to analyze
os.system("unzip ncbi_dataset.zip")
os.system("cp ncbi_dataset/data/cds.fna .")
#Find total number of CDS for first output to .log file using subprocess module, output if captured in genomic_CDS_num
genomic_CDS_num_command = subprocess.run("grep -c '>' cds.fna", shell=True, capture_output=True, text=True)
genomic_CDS_num = genomic_CDS_num_command.stdout.strip()
#Write first output of .log file
with open("PipelineProject.log","a+") as f:
	f.write("The HCMV genome (NC_006273.2) has " +str( genomic_CDS_num) + " CDS.\n")
	f.write("\n")

#Building index from HCMV
reference_transcriptome = "kallisto index -i index.idx cds.fna"
os.system(reference_transcriptome)

#First quantification
#Kallisto quantification is used for analysis. The previously assembled index is used as well as the paired FASTQ files that were previously downloaded
quantification_one_2dpi = "kallisto quant -i index.idx -o quantification_results_one_2dpi -b 10 -t 2 SRR5660030_1.fastq SRR5660030_2.fastq"
#The system then runs this command and changes to the output directory to read results from the abduncance tsv file
os.system(quantification_one_2dpi)
os.chdir("quantification_results_one_2dpi")
with open("abundance.tsv","r")as f:
#Lines are read and then each value is split into a list, newline characters are removed and stastitiscal analysis is done
	lines = f.readlines()[1:]
	base_results = [(i.split("\t")[4]) for i in lines]
	results = [float(x.strip()) for x in base_results]
#Function to calculate median result of list
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
#Finding the minimum, median, mean, and maximum values of the TPMs
minimum = min(results)
median = find_median(results) 
sum_results = sum(results)
length_results = len(results)
mean = sum_results/length_results
maximum = max(results)
#Return to main directory and write to the output file
os.chdir("..")
with open("PipelineProject.log", "a+") as f:
    f.write("sample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n")
    f.write("Donor 1\t" + "2dpi\t" + str(minimum) + "\t" + str(median) + "\t" + str(mean) + "\t" + str(maximum) + "\n")

#Second Quantification
quantification_one_6dpi = "kallisto quant -i index.idx -o quantification_results_one_6dpi -b 10 -t 2 SRR5660033_1.fastq SRR5660033_2.fastq"
os.system(quantification_one_6dpi)
os.chdir("quantification_results_one_6dpi")
with open("abundance.tsv","r")as f:
        lines = f.readlines()[1:]
        base_results = [(i.split("\t")[4]) for i in lines]
        results = [float(x.strip()) for x in base_results]
minimum = min(results)
median = find_median(results) 
sum_results = sum(results)
length_results = len(results)
mean = sum_results/length_results
maximum = max(results)
os.chdir("..")
with open("PipelineProject.log", "a+") as f:
    f.write("Donor 1\t" + "6dpi\t" + str(minimum) + "\t" + str(median) + "\t" + str(mean) + "\t" + str(maximum) + "\n")

#Third Quantification
quantification_three_2dpi = "kallisto quant -i index.idx -o quantification_results_three_2dpi -b 10 -t 2 SRR5660044_1.fastq SRR5660044_2.fastq"
os.system(quantification_three_2dpi)
os.chdir("quantification_results_three_2dpi")
with open("abundance.tsv","r")as f:
        lines = f.readlines()[1:]
        base_results = [(i.split("\t")[4]) for i in lines]
        results = [float(x.strip()) for x in base_results]
minimum = min(results)
median = find_median(results) 
sum_results = sum(results)
length_results = len(results)
mean = sum_results/length_results
maximum = max(results)
os.chdir("..")
with open("PipelineProject.log", "a+") as f:
    f.write("Donor 3\t" + "2dpi\t" + str(minimum) + "\t" + str(median) + "\t" + str(mean) + "\t" + str(maximum) + "\n")

#Fourth Quantification
quantification_three_6dpi = "kallisto quant -i index.idx -o quantification_results_three_6dpi -b 10 -t 2 SRR5660045_1.fastq SRR5660045_2.fastq"
os.system(quantification_three_6dpi)
os.chdir("quantification_results_three_6dpi")
with open("abundance.tsv","r")as f:
        lines = f.readlines()[1:]
        base_results = [(i.split("\t")[4]) for i in lines]
        results = [float(x.strip()) for x in base_results]
minimum = min(results)
median = find_median(results) 
sum_results = sum(results)
length_results = len(results)
mean = sum_results/length_results
maximum = max(results)
os.chdir("..")
with open("PipelineProject.log", "a+") as f:
	 f.write("Donor 3\t" + "6dpi\t" + str(minimum) + "\t" + str(median) + "\t" + str(mean) + "\t" + str(maximum) + "\n")


#Creating tsv file for reading for sleuth input
with open("sleuth_input.tsv","a+") as f:
	f.write("sample\tdonor\tcondition\tpath\n")
	f.write("SRR5660030\t1\t2dpi\tquantification_results_one_2dpi\n")
	f.write("SRR5660033\t1\t6dpi\tquantification_results_one_6dpi\n")
	f.write("SRR5660044\t3\t2dpi\tquantification_results_three_2dpi\n")
	f.write("SRR5660045\t3\t6dpi\tquantification_results_three_6dpi\n")

#Run Rscript command, output will be a tsv file with significant values
os.system('Rscript sleuth_command.R')

#Extract data from Rscript output
with open("significant_values.tsv", "r") as f:
	lines = f.readlines()[1:]
	target_id = [i.split("\t")[1] for i in lines]
	test_stat_raw = [(i.split("\t")[2]) for i in lines]
	test_stat = [float(x.strip()) for x in test_stat_raw]
	pval_raw = [i.split("\t")[3] for i in lines]
	pval = [float(x.strip()) for x in pval_raw]
	qval_raw = [i.split("\t")[4] for i in lines]
	qval = [float(x.strip()) for x in qval_raw]

#Write output to log file
with open("PipelineProject.log", "a+") as f:
    f.write("target_id\ttest_stat\tpval\tqval\n")
    for i in range(len(target_id)):
   	 f.write(target_id[i] + "\t" + test_stat[i] + "\t" + pval[i] + "\t" + qval[i] + "\n")

