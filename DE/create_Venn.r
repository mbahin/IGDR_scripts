# Mathieu Bahin, 07/04/14

# Script to produce a Venn diagram.
# The input is a matrix.
# The outout is a 'png' image showing the input data into a Venn diagram.

# Loading library
suppressPackageStartupMessages(library(limma))

# Producing the Venn diagram
input <- read.csv('file.venn.csv',header=T,row.names=1)
distribution <- vennCounts(input)
png('vennDiag.png',width=800,height=800)
vennDiagram(distribution,names=colnames(input),circle.col="darkgreen",counts.col="blue")
sink('/dev/null')
dev.off()