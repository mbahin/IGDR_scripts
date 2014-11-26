# Mathieu Bahin, 07/04/14

# Script to produce a Venn diagram.
# The input is a csv file with 0 or 1 for each feature.
# The outout is a 'png' image showing the input data into a Venn diagram.

# Loading library
suppressPackageStartupMessages(library(limma))

# Getting paramter back
args <- commandArgs(trailingOnly = TRUE)
input <- args[1]

# Producing the Venn diagram
input <- read.csv(input,header=T,row.names=1)
distribution <- vennCounts(input)
png('vennDiag.png',width=800,height=800)
vennDiagram(distribution,names=colnames(input),circle.col="darkgreen",counts.col="blue")
sink('/dev/null')
dev.off()