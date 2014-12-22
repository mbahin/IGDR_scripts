# Loading library
print("Loading libraries...")
suppressPackageStartupMessages(library("DESeq2"))

# Getting parameters back
args <- commandArgs(trailingOnly = TRUE)
directory <- args[1]
htsfilter <- (args[2] == TRUE)
paired <- (args[3] == TRUE)

# Finding the data files in the directory
print("Loading count data...")
sampleFiles <- as.character(read.csv('../metadata.csv',header=T)[,1])

# Finding the factors from the file 'conditions.csv' in the directory
sampleCondition <- as.character(read.csv('../metadata.csv',header=TRUE)$condition)
if (paired) {
	samplePatient <- as.character(read.csv('../metadata.csv',header=T)$patient)
}

# Creating a data frame linking the count files and the fatcors
if (paired) {
	sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleCondition, patient = samplePatient)
} else {
	sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleCondition)
}
# Loading count data from HTSeq-count
if (paired) {
	ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design= ~ patient + condition)
} else {
	ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design= ~ condition)
}
ddsHTSeq

# Making sure 'control' is the first factor (alphabetical order by default)
#colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,levels=c("DOG","WOLF"))

# Applying DESeq differential analysis steps to the dataset
suppressMessages(dds <- DESeq(ddsHTSeq))

# Writing the normalized counts in 'file.normalized.counts'
write.csv(as.data.frame(counts(dds,normalized=TRUE)), file='file.normalized.counts')

# Extracting results from the DESeq analysis
if (htsfilter) {
	suppressPackageStartupMessages(library("HTSFilter"))
	print("HTSFilter normalisation and filtering...")
	filter <- HTSFilter(dds, s.len=25, plot=FALSE)$filteredData
	print("Computing count data differential analysis...")
	res <- results(filter, independentFiltering=FALSE)
} else {
	print("DESeq2 normalisation and filtering...")
	print("Computing count data differential analysis...")
	res <- results(dds)
}

# Ordering the results 	after the adjusted p-value
res <- res[order(res$padj),]
print("Overview of the output:")
head(res)

# Writing results in a file
print("Writing the output file...")
write.csv(as.data.frame(res), file='file.DESeq2.output.csv')

# Producing labelled PCA plot
print("Producing the labelled PCA plot... (PCA.png)")
	# Redefining plotPCA function to label the PCA plot
customPlotPCA = function (x, intgroup = "condition", ntop = 500) 
{
	library("lattice")
	suppressPackageStartupMessages(library("genefilter"))
    rv = rowVars(assay(x))
    select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca = prcomp(t(assay(x)[select, ]))
    fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop = FALSE]), 
        1, paste, collapse = " : "))
    names = sub("file.(.*).PE.picard.reorder.count","\\1",sampleFiles) # added row
    if (nlevels(fac) >= 3) 
        colours = brewer.pal(nlevels(fac), "Paired")
    else colours = c("lightgreen", "dodgerblue")
    xyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$x), 
        pch = 16, cex = 2, 
        
        panel=function(x, y, ...) {
			panel.xyplot(x, y, ...);
			ltext(x=x, y=y, labels=names, pos=1, offset=1, cex=0.8)
		} #added_row
        ,
        
        aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours), 
            text = list(levels(fac)), rep = FALSE)))
}

	# Computing plot
suppressMessages(rld <- rlogTransformation(dds, blind=TRUE))
png(filename='../PCA.png',width=600,height=600)
print(customPlotPCA(rld, intgroup=c("condition")))
sink('/dev/null')
dev.off()
sink()

# Producing sample-to-sample distance heatmap
print("Producing the sample-to-sample distance heatmap... (heatmap.png)")
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("gplots"))
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL) 
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
png(filename='../heatmap.png',width=800,height=800)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(20, 20), cexRow=1.5, cexCol=1.5)
sink('/dev/null')
dev.off()
