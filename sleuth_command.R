#Load libraries
library(dplyr)
library(sleuth)
#Read in table from four reads, ignoring header
stab <- read.table('sleuth_input.tsv', header = TRUE)

#Create object for sleuth, using full model
so <- sleuth_prep(stab)

#Fitting the data to different models, includes error and reduced fit
so <- sleuth_fit(so, ~condition, 'condition.fit')
so <- sleuth_fit(so, ~1, 'condition.re')

#Use likehood ratio test for differential analysis for gene expression
so <- sleuth_lrt(so, 'condition.re', 'condition.fit')

#Find significant results
sleuth_table <- sleuth_results(so, 'condition.re:condition.fit', 'lrt', show_all = FALSE)

#Filter to only results with qval less than 0.05, arrange p values in ascending order
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval)

#Extract target_id, test_stat, pval, and qval for output
sleuth_output <- dplyr::select(sleuth_significant, target_id, test_stat, pval, qval)

#Wrties values to table for extraction to python pipeline log file using tab as a delimeter
write.table(sleuth_output, 'significant_values.tsv', sep = '\t', row.names = FALSE)
