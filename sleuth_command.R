#Load library
library(sleuth)
stab = read.table("sample_table.txt",header=TRUE)
so = sleuth_prep(stab)
