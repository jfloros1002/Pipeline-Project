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
os.system("cp ncbi_dataset/data/GCF_000845245.1/cds_from_genomic.fna .")
#Find total number of CDS for first output to .log file using subprocess module, output is captured in genomic_CDS_num
genomic_CDS_num_command = subprocess.run("grep -c '>' cds_from_genomic.fna", shell=True, capture_output=True, text=True)
genomic_CDS_num = genomic_CDS_num_command.stdout.strip()
#Write first output of .log file
with open("PipelineProject.log","a+") as f:
	f.write("The HCMV genome (NC_006273.2) has " +str( genomic_CDS_num) + " CDS.\n")
	f.write("\n")

#Building index from HCMV
reference_transcriptome = "kallisto index -i index.idx cds_from_genomic.fna"
os.system(reference_transcriptome)

#First quantification
#Kallisto quantification is used for analysis. The previously assembled index is used as well as the paired FASTQ files that were previously downloaded
quantification_one_2dpi = "kallisto quant -i index.idx -o quantification_results_one_2dpi -b 10 -t 2 *SRR5660030_1.fastq *SRR5660030_2.fastq"
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
quantification_one_6dpi = "kallisto quant -i index.idx -o quantification_results_one_6dpi -b 10 -t 2 *SRR5660033_1.fastq *SRR5660033_2.fastq"
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
quantification_three_2dpi = "kallisto quant -i index.idx -o quantification_results_three_2dpi -b 10 -t 2 *SRR5660044_1.fastq *SRR5660044_2.fastq"
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
quantification_three_6dpi = "kallisto quant -i index.idx -o quantification_results_three_6dpi -b 10 -t 2 *SRR5660045_1.fastq *SRR5660045_2.fastq"
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
	f.write("\n")

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
	target_id = [i.split("\t")[0] for i in lines]
	test_stat_raw = [(i.split("\t")[1]) for i in lines]
	test_stat = [float(x.strip()) for x in test_stat_raw]
	pval_raw = [i.split("\t")[2] for i in lines]
	pval = [float(x.strip()) for x in pval_raw]
	qval_raw = [i.split("\t")[3] for i in lines]
	qval = [float(x.strip()) for x in qval_raw]

#Write output to log file
with open("PipelineProject.log", "a+") as f:
    f.write("target_id\ttest_stat\tpval\tqval\n")
    for i in range(len(target_id)):
   	 f.write(str(target_id[i]) + "\t" +str(test_stat[i]) + "\t" +str(pval[i]) + "\t" +str(qval[i]) + "\n")
    f.write("\n")

#Build Index for Bowtie2
os.system("bowtie2-build ncbi_dataset/data/GCF_000845245.1/cds_from_genomic.fna HCMV")

#Running Bowtie2 on each paired read
os.system("bowtie2 --quiet -x HCMV -1 *SRR5660030_1.fastq -2 *SRR5660030_2.fastq -s SRR5660030_aligned.sam --al-conc-gz SRR5660030_mapped_%.fq.gz")
os.system("bowtie2 --quiet -x HCMV -1 *SRR5660033_1.fastq -2 *SRR5660033_2.fastq -s SRR5660033_aligned.sam --al-conc-gz SRR5660033_mapped_%.fq.gz")
os.system("bowtie2 --quiet -x HCMV -1 *SRR5660044_1.fastq -2 *SRR5660044_2.fastq -s SRR5660044_aligned.sam --al-conc-gz SRR5660044_mapped_%.fq.gz")
os.system("bowtie2 --quiet -x HCMV -1 *SRR5660045_1.fastq -2 *SRR5660045_2.fastq -s SRR5660045_aligned.sam --al-conc-gz SRR5660045_mapped_%.fq.gz")

#Unzip all reads
os.system("gunzip SRR5660030_mapped_1.fq")
os.system("gunzip SRR5660030_mapped_2.fq")
os.system("gunzip SRR5660033_mapped_1.fq")
os.system("gunzip SRR5660033_mapped_2.fq")
os.system("gunzip SRR5660044_mapped_1.fq")
os.system("gunzip SRR5660044_mapped_2.fq")
os.system("gunzip SRR5660045_mapped_1.fq")
os.system("gunzip SRR5660045_mapped_2.fq")

#Each set of paired reads has only one + to represent the line with ASCII scores, this means that we can count the number of "+" to determined how many paired reads there are
#Find the output of a grep command for raw reads and mapped reads and store to variables for each sample/conditions
length_Donor1_2dpi_all_raw = subprocess.run("grep -c '+' SRR5660030_1.fastq", shell=True, capture_output=True, text=True)
length_Donor1_2dpi_aligned_raw = subprocess.run("grep -c '+' SRR5660030_mapped_1.fq", shell=True, capture_output=True, text=True)
length_Donor1_2dpi_all = length_Donor1_2dpi_all_raw.stdout.strip()
length_Donor1_2dpi_aligned = length_Donor1_2dpi_aligned_raw.stdout.strip()

length_Donor1_6dpi_all_raw = subprocess.run("grep -c '+' SRR5660033_1.fastq", shell=True, capture_output=True, text=True)
length_Donor1_6dpi_aligned_raw = subprocess.run("grep -c '+' SRR5660033_mapped_1.fq", shell=True, capture_output=True, text=True)
length_Donor1_6dpi_all = length_Donor1_6dpi_all_raw.stdout.strip()
length_Donor1_6dpi_aligned = length_Donor1_6dpi_aligned_raw.stdout.strip()

length_Donor3_2dpi_all_raw = subprocess.run("grep -c '+' SRR5660044_1.fastq", shell=True, capture_output=True, text=True)
length_Donor3_2dpi_aligned_raw = subprocess.run("grep -c '+' SRR5660044_mapped_1.fq", shell=True, capture_output=True, text=True)
length_Donor3_2dpi_all = length_Donor3_2dpi_all_raw.stdout.strip()
length_Donor3_2dpi_aligned = length_Donor3_2dpi_aligned_raw.stdout.strip()

length_Donor3_6dpi_all_raw = subprocess.run("grep -c '+' SRR5660045_1.fastq", shell=True, capture_output=True, text=True)
length_Donor3_6dpi_aligned_raw = subprocess.run("grep -c '+' SRR5660045_mapped_1.fq", shell=True, capture_output=True, text=True)
length_Donor3_6dpi_all = length_Donor3_6dpi_all_raw.stdout.strip()
length_Donor3_6dpi_aligned = length_Donor3_6dpi_aligned_raw.stdout.strip()

#Output found values to PipelineProject.log
with open("PipelineProject.log","a+") as f:
	f.write("Donor 1 (2dpi) had " + str(length_Donor1_2dpi_all) + " read pairs before Bowtie2 filtering and " + str(length_Donor1_2dpi_aligned) + " read pairs after.\n")
	f.write("Donor 1 (6dpi) had " + str(length_Donor1_6dpi_all) + " read pairs before Bowtie2 filtering and " + str(length_Donor1_6dpi_aligned) + " read pairs after.\n")
	f.write("Donor 3 (2dpi) had " + str(length_Donor3_2dpi_all) + " read pairs before Bowtie2 filtering and " + str(length_Donor3_2dpi_aligned) + " read pairs after.\n")
	f.write("Donor 3 (6dpi) had " + str(length_Donor3_6dpi_all) + " read pairs before Bowtie2 filtering and " + str(length_Donor3_6dpi_aligned) + " read pairs after.\n")
	f.write("\n")

#Commands for running spades with each set of Donor reads. A k size of 77 was requested, so this is what k is set to
Donor1_spades_command = "spades.py -k 77 -t 2 --only-assembler --pe-1 1 SRR5660030_mapped_1.fq --pe-2 1 SRR5660030_mapped_2.fq --pe-1 2 SRR5660033_mapped_1.fq --pe-2 2 SRR5660033_mapped_2.fq -o Donor1_assembly/"
Donor3_spades_command = "spades.py -k 77 -t 2 --only-assembler --pe-1 1 SRR5660044_mapped_1.fq --pe-2 1 SRR5660044_mapped_2.fq --pe-1 2 SRR5660045_mapped_1.fq --pe-2 2 SRR5660045_mapped_2.fq -o Donor3_assmebly/"

#Running spades commands
os.system(Donor1_spades_command)
os.system(Donor3_spades_command)

#Writing commands to PipelineProject.log
with open("PipelineProject.log","a+") as f:
	f.write("SPAdes command for Donor 1: " + str(Donor1_spades_command) + "\n")
	f.write("SPAdes command for Donor 3: " + str(Donor3_spades_command) + "\n")
	f.write("\n")

os.chdir("Donor1_assembly")
with open("scaffolds.fasta","r+") as f:
	lines = f.readlines()
	total_strand = ""
	current_strand = ""
	while ">" not in current_strand:
		for i in range(1,len(lines)+1):
			if ">" in current_strand:
				break
			else:
				current_strand = lines[i]
				current_strand = current_strand.strip("\n")
				print(current_strand)
				total_strand += str(current_strand)

print(total_strand)
