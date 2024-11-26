#!/usr/bin/env Rscript

# system specifications
## Architecture:                    x86_64
## CPU op-mode(s):                  32-bit, 64-bit
## CPU(s):                          8
## Thread(s) per core:              2
## Core(s) per socket:              4
## Model name:                      Intel(R) Core(TM) i7-10510U CPU @ 1.80GHz
## R:                               4.3.0
## RStudio:                         2022.02.3 Build 492

# goals: calculate metrics comparing variants (SNPs or InDels) comparing reference and workflow files
## recall=sensitivity=TP/(TP+FN)
## precision=TP/(TP+FP)
## specificity=TN/(FP+TN)
## accuracy=(TP+TN)/(TP+TN+FP+FN)
## F-score (2*TP)/((2+TP)+FP+FN)

# controlled samples (input_Reference_Variants.tsv and input_Workflow_Variants.tsv)
## Reference	   Workflow	       Results    Positions (S1439.014)
## genotype	     genotype	       TP         254
## genotype	     other genotype  FP         653
## undetected	   genotype	       FP         748
## missing	     genotype	       MD         25014
## genotype	     undetected	     FN         35124
## genotype	     missing	       MD         7589
## undetected	   undetected	     TN         47000
## undetected	   missing	       MD         5012
## missing	     missing         MD         43000
## missing	     undetected	     MD         45000

# expected outcomes (Rscript VariantsMetricsReference:1.0.R -r input_Reference_Variants.tsv -w input_Workflow_Variants.tsv -g 50000 -o output_Reference_)
#    samples MDr MDw MD TP    TN FP FN  recall precision  specificity accuracy  Fscore
#1 S1239.014   3   3  5  1 49991  2  1 0.50000   0.33333      0.99996  0.99994 0.40000
#2 S1339.010   3   3  5  1 49991  2  1 0.50000   0.33333      0.99996  0.99994 0.40000
#3 S1339.020   3   2  3  4 49993  0  0 1.00000   1.00000      1.00000  1.00000 1.00000
#4 S1442.017   3   2  3  4 49993  0  0 1.00000   1.00000      1.00000  1.00000 1.00000

# skip lines related to installation of libraries because there are supposed to be already installed
skip_instalation <- scan(what="character", quiet = TRUE)
# install libraries
install.packages("remotes") # version 2.4.2
require(remotes)
install_version("argparse", version = "2.2.3", repos = "https://cloud.r-project.org")

# load packages avoiding warning messages
suppressPackageStartupMessages(library(argparse))
# clean environment
rm(list=ls())
# clean graphical device
graphics.off()
# use English language
invisible(capture.output(Sys.setlocale("LC_TIME", "C")))
# keep in mind start time
start.time <- Sys.time()
# keep in mind start time as human readable
start.time.readable <- format(Sys.time(), "%X %a %b %d %Y")

# create parser object
parser <- ArgumentParser(description = "This script computes performance metrics (i.e. TP, TN, FP, FN, recall, precision, specificity, accuracy and Fscore) comparing genotypes ('undetected' stands for 'undetected genotype' and 'missing' stands for 'missing data') from expected (i.e. reference file) and predicted (i.e. workflow file) variants (i.e. SNP and InDels).")
# add in parser all the arguments sorted by order of appearance
parser$add_argument("-r", "--reference", type="character", metavar="CHARACTERS", required=TRUE, dest="reference", 
                    help="Reference input file with an absolute or relative path (tab-separated values). First column: positions of variants (header: whatever). Second column: profiles of genotypes (header: whatever). [MANDATORY]")
parser$add_argument("-w", "--workflow", type="character", metavar="CHARACTERS", required=TRUE, dest="workflow", 
                    help="Workflow input file with an absolute or relative path (tab-separated values). First column: positions of variants (header: whatever). Other columns: profiles of genotypes (header: sample identifiers). [MANDATORY]")
parser$add_argument("-g", "--genome", type="integer", metavar="INTEGER", required=TRUE, dest="genome", 
                    help="Size of the reference genome (bases). [MANDATORY]")
parser$add_argument("-o", "--prefix", type="character", default="output_", metavar="CHARACTERS", dest="prefix", 
                    help="Absolute or relative output path with or without output file prefix. [OPTIONAL, default %(default)s]")
parser$add_argument("-b", "--backup", type="logical", default=FALSE, metavar="LOGICAL", dest="backup", 
                    help="Save an external representation of R objects (i.e. saved_data.RData) and a short-cut of the current workspace (i.e. saved_images.RData). [OPTIONAL, default %(default)s]")
# get command line options, if help option encountered print help and exit, otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

# retrieve genome size
genome <- args$genome
# read dataframes
df.reference <- read.table(args$reference, dec = ".", header=TRUE, sep = "\t", quote = "", check.names = FALSE)
df.workflow.samples <- read.table(args$workflow, dec = ".", header=TRUE, sep = "\t", quote = "", check.names = FALSE)

# R console tests
#rm(list=ls())
#setwd("/home/IZSNT/n.radomski/Documents/RstudioWorkingDirectory/VariantsMetrics-20241112")
#genome <- 50000
#df.reference <- read.table("input_Reference_Variants.tsv", dec = ".", header=TRUE, sep = "\t", quote = "", check.names = FALSE)
#df.workflow.samples <- read.table("input_Workflow_Variants.tsv", dec = ".", header=TRUE, sep = "\t", quote = "", check.names = FALSE)

# prepare help messages
help1 <- "Help: Rscript VariantsMetricsReference.R -h"
help2 <- "Help: Rscript VariantsMetricsReference.R --help"

# format checking of mandatory input files
## test the number of columns of the reference input file
if ((ncol(df.reference) != 2) == TRUE) {
  stop("The reference input file must harbor two columns", "\n", help1, "\n", help2)
  }
## test the number of columns of the workflow input file
if ((ncol(df.workflow.samples) < 2) == TRUE) {
  stop("The workflow input file must harbor at least two columns", "\n", help1, "\n", help2)
  }
## test the number of rows of the reference input file
if ((nrow(df.reference) >= 1) == FALSE) {
  stop("The reference input file must harbor at least one row", "\n", help1, "\n", help2)
  }
## test the number of rows of the workflow input file
if ((nrow(df.workflow.samples) >= 1) == FALSE) {
  stop("The workflow input file must harbor at least one row", "\n", help1, "\n", help2)
  }
## test the absence of empty cells in the first column of the reference input file
if (any((is.na(df.reference[,1]))) == TRUE) {
  stop("The first column of the reference input file must not harbor empty cells", "\n", help1, "\n", help2)
  }
## test the absence of empty cells in the first column of the workflow input file
if (any((is.na(df.workflow.samples[,1]))) == TRUE) {
  stop("The first column of the workflow input file must not harbor empty cells", "\n", help1, "\n", help2)
  }
## test the absence of empty cells in the headers of the reference input file
if (("" %in% colnames(df.reference)) == TRUE) {
  stop("The header of the reference input file must not harbor empty cells", "\n", help1, "\n", help2)
}
## test the absence of empty cells in the headers of the workflow input file
if (("" %in% colnames(df.workflow.samples)) == TRUE) {
  stop("The header of the the workflow input file must not harbor empty cells", "\n", help1, "\n", help2)
}
## test the presence of integers in the first column of the reference input file
if ((is.integer(df.reference[,1])) == FALSE) {
  stop("The first column of the reference input file must only harbor integers", "\n", help1, "\n", help2)
  }
## test the presence of integers in the first column of the workflow input file
if ((is.integer(df.workflow.samples[,1])) == FALSE) {
  stop("The first column of the workflow input file must only harbor integers", "\n", help1, "\n", help2)
  }

# checking of positions
## test the absence of duplicated positions in the reference input file
if ((any((duplicated(df.reference[,1])) == TRUE)) == FALSE) {
  cat("The absence of duplicated positions from the reference input file was correctly controlled", "\n", sep = "")
} else {
  stop("The positions from the reference input file must not harbor duplicates", "\n", help1, "\n", help2)
}
## test the absence of duplicated positions in the workflow input file
if ((any((duplicated(df.workflow.samples[,1])) == TRUE)) == FALSE) {
  cat("The absence of duplicated positions from the workflow input file was correctly controlled", "\n", sep = "")
} else {
  stop("The positions from the workflow input file must not harbor duplicates", "\n", help1, "\n", help2)
}

# checking of samples
## retrieve samples names
samples <- names(df.workflow.samples)[2:ncol(df.workflow.samples)]
## test the absence of sample identifier duplicates in the workflow input file
if ((any((duplicated(samples)) == TRUE)) == FALSE) {
  cat("The absence of duplicated samples from the workflow input file was correctly controlled", "\n", sep = "")
} else {
  stop("The samples from the workflow input file must not harbor duplicates", "\n", help1, "\n", help2)
}

# checking of reference genome size
## test the reference genome size compared to first column of the reference input file
if (genome >= (max(df.reference[,1]))) {
  cat("The reference genome size was correctly controlled as higher than, or equal to, the highest position of the reference input file", "\n", sep = "")
} else {
  stop("The reference genome size must be higher than, or equal to, the highest position of the reference input file", "\n", help1, "\n", help2)
}
## test the reference genome size compared to first column of the workflow input file
if (genome >= (max(df.workflow.samples[,1]))) {
  cat("The reference genome size was correctly controlled as higher than, or equal to, the highest position of the workflow input file", "\n", sep = "")
} else {
  stop("The reference genome size must be higher than, or equal to, the highest position of the workflow input file", "\n", help1, "\n", help2)
}

# harmonize the reference and workflow dataframes
## rename columns
colnames(df.reference) <- c("position","reference")
colnames(df.workflow.samples)[1] = "position"
## join
df <- merge(df.reference, df.workflow.samples, by = "position", all = TRUE, sort = FALSE)
## replace <NA> by "undetected"
df[is.na(df)] = "undetected"
## replace empty cells by "missing"
df[df == ''] <- "missing"
## sort by positions
df <- df[order(df$position, decreasing = FALSE), ]
## retrieve sorted positions
positions <- df$position
## retrieve reference genotype
r <- df$reference

# TP: The variant genotypes from the workflow and reference are not undetected (i.e. undetected), not missing (i.e. missing) and identical
## Reference	   Workflow	       Results    Positions (2024.EXT.0126.1439.014)
## genotype	   genotype	         TP         254
## compare workflow (w) and reference (r)
### create an empty boolean vector
vTP <- logical()
### loop over columns
for (w in df[, 3:ncol(df)]) {
  output <- ifelse(((w != "undetected") & 
                      (w != "missing") & 
                      (r != "undetected") & 
                      (r != "missing") &
                      (r == w)), TRUE, FALSE)
  vTP <- c(vTP, output) # add the output vectors into the empty vector
} 
### create a matrix adding the vector by row independently of the sample amount
mat.TP <- matrix(vTP, nrow = nrow(df), ncol = (ncol(df)-2), byrow = FALSE)
### transform as dataframe
df.TP <- as.data.frame(mat.TP)
### add samples names
colnames(df.TP) <- samples
### add positions to check results
df.TP <- cbind(positions, df.TP)
## count TRUE
### create an empty integer vector
TP <- integer()
### loop over columns
for (s in df.TP[, 2:ncol(df.TP)]) {
  output <- sum(s, na.rm=TRUE)
  TP <- c(TP, output)
} 
## transform as dataframe
df.TP.count <- as.data.frame(TP)
## add samples names
df.TP.count <- cbind(samples, df.TP.count)

# FN: The variant genotype from the workflow is undetected (i.e. undetected) but not missing (i.e. missing), while it is detected into the reference
## Reference	   Workflow	       Results    Positions (2024.EXT.0126.1439.014)
## genotype	   undetected	       FN         35124
## compare workflow (w) and reference (r)
### create an empty boolean vector
vFN <- logical()
### loop over columns
for (w in df[, 3:ncol(df)]) {
  output <- ifelse((w == "undetected") & ((r != "undetected") & (r != "missing")), TRUE, FALSE)
  vFN <- c(vFN, output) # add the output vectors into the empty vector
} 
### create a matrix adding the vector by row independently of the sample amount
mat.FN <- matrix(vFN, nrow = nrow(df), ncol = (ncol(df)-2), byrow = FALSE)
### transform as dataframe
df.FN <- as.data.frame(mat.FN)
### add samples names
colnames(df.FN) <- samples
### add positions to check results
df.FN <- cbind(positions, df.FN)
## count TRUE
### create an empty integer vector
FN <- integer()
### loop over columns
for (s in df.FN[, 2:ncol(df.FN)]) {
  output <- sum(s, na.rm=TRUE)
  FN <- c(FN, output)
} 
## transform as dataframe
df.FN.count <- as.data.frame(FN)
## add samples names
df.FN.count <- cbind(samples, df.FN.count)

# FP: The variant genotype from the workflow is detected (i.e. not "undetected" and not "missing"), while it is undetected (i.e. undetected) or different, but not missing (i.e. missing) into the reference
## Reference	   Workflow	       Results    Positions (2024.EXT.0126.1439.014)
## genotype	     other genotype	   FP         653
## undetected 	 genotype	         FP         748
## compare workflow (w) and reference (r)
### create an empty boolean vector
vFP <- logical()
### loop over columns
for (w in df[, 3:ncol(df)]) {
  output <- ifelse((w != "undetected") & (w != "missing") & (r != "undetected") & (r != "missing") & (r != w) | 
                     (w != "undetected") & (w != "missing") & (r == "undetected"), TRUE, FALSE)
  vFP <- c(vFP, output) # add the output vectors into the empty vector
} 
### create a matrix adding the vector by row independently of the sample amount
mat.FP <- matrix(vFP, nrow = nrow(df), ncol = (ncol(df)-2), byrow = FALSE)
### transform as dataframe
df.FP <- as.data.frame(mat.FP)
### add samples names
colnames(df.FP) <- samples
### add positions to check results
df.FP <- cbind(positions, df.FP)
## count TRUE
### create an empty integer vector
FP <- integer()
### loop over columns
for (s in df.FP[, 2:ncol(df.FP)]) {
  output <- sum(s, na.rm=TRUE)
  FP <- c(FP, output)
} 
## transform as dataframe
df.FP.count <- as.data.frame(FP)
## add samples names
df.FP.count <- cbind(samples, df.FP.count)

# MD (missing data): The variant genotype from the reference and/or workflow is missing (i.e. "missing")
## Reference	   Workflow	       Results    Positions (2024.EXT.0126.1439.014)
## genotype	     missing	       MD         7589
## missing	     genotype	       MD         25014
## undetected	   missing	       MD         5012
## missing	     missing         MD         43000
## missing	     undetected	     MD         45000
### create an empty boolean vector
vMD <- logical()
### loop over columns
for (w in df[, 3:ncol(df)]) {
  output <- ifelse(((r == "missing") & (w != "missing")) | 
                     ((r != "missing") & (w == "missing")) | 
                     ((r == "missing") & (w == "missing")) , TRUE, FALSE)
  vMD <- c(vMD, output) # add the output vectors into the empty vector
} 
### create a matrix adding the vector by row independently of the sample amount
mat.MD <- matrix(vMD, nrow = nrow(df), ncol = (ncol(df)-2), byrow = FALSE)
### transform as dataframe
df.MD <- as.data.frame(mat.MD)
### add samples names
colnames(df.MD) <- samples
### add positions to check results
df.MD <- cbind(positions, df.MD)
## count TRUE
### create an empty integer vector
MD <- integer()
### loop over columns
for (s in df.MD[, 2:ncol(df.MD)]) {
  output <- sum(s, na.rm=TRUE)
  MD <- c(MD, output)
} 
## transform as dataframe
df.MD.count <- as.data.frame(MD)
## add samples names
df.MD.count <- cbind(samples, df.MD.count)

# MDr (missing data in reference): The variant genotype from the reference is missing (i.e. "missing")
## Reference	   Workflow	       Results    Positions (2024.EXT.0126.1439.014)
## missing	     genotype	       MD         25014
## missing	     missing         MD         43000
## missing	     undetected	     MD         45000
### create an empty boolean vector
vMDr <- logical()
### loop over columns
for (w in df[, 3:ncol(df)]) {
  output <- ifelse(((r == "missing") & (w != "missing")) | 
                     ((r == "missing") & (w == "missing")) , TRUE, FALSE)
  vMDr <- c(vMDr, output) # add the output vectors into the empty vector
} 
### create a matrix adding the vector by row independently of the sample amount
mat.MDr <- matrix(vMDr, nrow = nrow(df), ncol = (ncol(df)-2), byrow = FALSE)
### transform as dataframe
df.MDr <- as.data.frame(mat.MDr)
### add samples names
colnames(df.MDr) <- samples
### add positions to check results
df.MDr <- cbind(positions, df.MDr)
## count TRUE
### create an empty integer vector
MDr <- integer()
### loop over columns
for (s in df.MDr[, 2:ncol(df.MDr)]) {
  output <- sum(s, na.rm=TRUE)
  MDr <- c(MDr, output)
} 
## transform as dataframe
df.MDr.count <- as.data.frame(MDr)
## add samples names
df.MDr.count <- cbind(samples, df.MDr.count)

# MDw (missing data): The variant genotype from the workflow is missing (i.e. "missing")
## Reference	   Workflow	       Results    Positions (2024.EXT.0126.1439.014)
## genotype	     missing	       MD         7589
## undetected	   missing	       MD         5012
## missing	     missing         MD         43000
### create an empty boolean vector
vMDw <- logical()
### loop over columns
for (w in df[, 3:ncol(df)]) {
  output <- ifelse(((r != "missing") & (w == "missing")) | 
                     ((r == "missing") & (w == "missing")) , TRUE, FALSE)
  vMDw <- c(vMDw, output) # add the output vectors into the empty vector
} 
### create a matrix adding the vector by row independently of the sample amount
mat.MDw <- matrix(vMDw, nrow = nrow(df), ncol = (ncol(df)-2), byrow = FALSE)
### transform as dataframe
df.MDw <- as.data.frame(mat.MDw)
### add samples names
colnames(df.MDw) <- samples
### add positions to check results
df.MDw <- cbind(positions, df.MDw)
## count TRUE
### create an empty integer vector
MDw <- integer()
### loop over columns
for (s in df.MDw[, 2:ncol(df.MDw)]) {
  output <- sum(s, na.rm=TRUE)
  MDw <- c(MDw, output)
} 
## transform as dataframe
df.MDw.count <- as.data.frame(MDw)
## add samples names
df.MDw.count <- cbind(samples, df.MDw.count)

# TN: The variant genotype is not detected neither in the reference (i.e. "undetected") not in the workflow (i.e. "undetected")
## Reference	   Workflow	       Results    Positions (2024.EXT.0126.1439.014)
## undetected	   undetected	       TN         47000
## merge count from TP, FN, FP and MD
df.TP.FN.count <- merge(df.TP.count, df.FN.count, by="samples", all = FALSE)
df.TP.FN.FP.count <- merge(df.TP.FN.count, df.FP.count, by="samples", all = FALSE)
df.TP.FN.FP.MD.count <- merge(df.TP.FN.FP.count, df.MD.count, by="samples", all = FALSE)
df.TP.FN.FP.MD.MDr.count <- merge(df.TP.FN.FP.MD.count, df.MDr.count, by="samples", all = FALSE)
df.results <- merge(df.TP.FN.FP.MD.MDr.count, df.MDw.count, by="samples", all = FALSE)
## calculate TN
df.results$TN <- genome-df.results$TP-df.results$FN-df.results$FP-df.results$MD
## reorder columns
df.results <- df.results[, c("samples","MDr", "MDw", "MD", "TP", "TN", "FP", "FN")]

# calculate metrics
## recall=sensitivity=TP/(TP+FN)
df.results$recall <- df.results$TP/(df.results$TP+df.results$FN)
## precision=TP/(TP+FP)
df.results$precision <- df.results$TP/(df.results$TP+df.results$FP)
## specificity=TN/(FP+TN)
df.results$specificity <- df.results$TN/(df.results$FP+df.results$TN)
## accuracy=(TP+TN)/(TP+TN+FP+FN)
df.results$accuracy <- (df.results$TP+df.results$TN)/(df.results$TP+df.results$TN+df.results$FP+df.results$FN)
## F-score (2*TP)/((2+TP)+FP+FN)
df.results$Fscore <- (2*df.results$TP)/((2*df.results$TP)+df.results$FP+df.results$FN)
# transform as numerical values
df.results[,6:10] <- sapply(df.results[,6:10], as.numeric)
# round
options(digits=6)
df.rounded.results <- df.results
df.rounded.results$recall <- format(round(df.rounded.results$recall, 5), nsmall = 5)
df.rounded.results$precision <- format(round(df.rounded.results$precision, 5), nsmall = 5)
df.rounded.results$specificity <- format(round(df.rounded.results$specificity, 5), nsmall = 5)
df.rounded.results$accuracy <- format(round(df.rounded.results$accuracy, 5), nsmall = 5)
df.rounded.results$Fscore <- format(round(df.rounded.results$Fscore, 5), nsmall = 5)

# keep in mind end time
end.time <- Sys.time()

# calculate execution time
time.taken <- difftime(end.time, start.time, units="secs")

# output
## dataframes
write.table(df.TP, file = paste(args$prefix, "dfTP.tsv", sep = ""), append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(df.FN, file = paste(args$prefix, "dfFN.tsv", sep = ""), append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(df.FP, file = paste(args$prefix, "dfFP.tsv", sep = ""), append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(df.MDr, file = paste(args$prefix, "dfMDr.tsv", sep = ""), append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(df.MDw, file = paste(args$prefix, "dfMDw.tsv", sep = ""), append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(df.MD, file = paste(args$prefix, "dfMD.tsv", sep = ""), append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
## metrics
write.table(df.rounded.results, file = paste(args$prefix, "metrix.tsv", sep = ""), append = FALSE, quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
## RData
if (isTRUE(args$backup)){
  ## to load with load("output_saved_data.RData")
  save(list = ls(), file = paste(args$prefix, "saved_data.RData", sep = ""))
  ## to load with load("output_saved_images.RData")
  save.image(file = paste(args$prefix, "saved_images.RData", sep = ""))
}
## log
### open
sink(paste(args$prefix, "log.txt", sep = ""))
### print
cat("Start time:", start.time.readable, "\n")
cat("R version:", strsplit(version[['version.string']], ' ')[[1]][3], "\n") # 4.4.0
cat("remotes version:", getNamespaceVersion("remotes"), "\n") # 2.5.0
cat("argparse version:", getNamespaceVersion("argparse"), "\n") # 2.2.3
cat("The absence of duplicated positions from the reference input file was correctly controlled", "\n", sep = "")
cat("The absence of duplicated samples from the workflow input file was correctly controlled", "\n", sep = "")
cat("The reference genome size was correctly controlled as higher than, or equal to, the highest position of the reference input file", "\n", sep = "")
cat("The reference genome size was correctly controlled as higher than, or equal to, the highest position of the workflow input file", "\n", sep = "")
cat("MDr stands for missing data in the reference", "\n", sep = "")
cat("MDw stands for missing data in the workflow", "\n", sep = "")
cat("MD stands for missing data in the reference and/or workflow", "\n", sep = "")
cat("TP stands for true positive", "\n", sep = "")
cat("TN stands for true negative", "\n", sep = "")
cat("FP stands for false positive", "\n", sep = "")
cat("FN stands for false negative", "\n", sep = "")
cat("Running time (seconds):", time.taken, "\n")
cat("Outcomes: ", args$prefix,"\n", sep = "")
cat("Developped by Nicolas Radomski since October 2024", "\n")
cat("GitHub: https://github.com/Nicolas-Radomski/VariantsMetrics", "\n")
cat("Docker: https://hub.docker.com/r/nicolasradomski/variantsmetricsreference", "\n")
### close
sink()

# add messages
cat("Running time (seconds):", time.taken, "\n")
cat("Outcomes: ", args$prefix,"\n", sep = "")
cat("Developped by Nicolas Radomski since October 2024", "\n")
cat("GitHub: https://github.com/Nicolas-Radomski/VariantsMetrics", "\n")
cat("Docker: https://hub.docker.com/r/nicolasradomski/variantsmetricsreference", "\n")
