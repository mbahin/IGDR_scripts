# Loading library
print("Loading libraries...")
suppressPackageStartupMessages(library(edgeR))

# Getting parameters back
args <- commandArgs(trailingOnly = TRUE)
htsfilter <- (args[1] == TRUE)
glm <- (args[2] == TRUE)
paired <- (args[3] == TRUE)

# Finding the data files in the directory
sampleFiles <- list.files(path='../Counts/')
counts <- readDGE(sampleFiles,path='../Counts/',header=F)

# Finding the conditions from the file 'conditions.csv' in the directory
sampleCondition <- as.character(read.csv('../conditions.csv',header=FALSE)$V1)

# Creating a DGEList
y <- DGEList(counts=counts,group=sampleCondition)

# Calculating normalization factors to scale the raw library sizes
print("Normalisation and filtering...")
y <- calcNormFactors(y)

if (!glm) {
	# Estimating the dataset dispersion
	print('Estimating dataset dispersion...')
	y <- estimateCommonDisp(y)
	y <- estimateTagwiseDisp(y)
	
	# Writing the normalized counts in 'file.normalized.counts'
	write.csv(as.data.frame(y$pseudo.counts), file='file.normalized.counts')

	# Computing the differences between counts
	print("Computing count data differential analysis...")
	et <- exactTest(y)
	if (htsfilter) {
		print('Applying HTSFilter...')
		suppressPackageStartupMessages(library("HTSFilter"))
		et <- HTSFilter(et, DGEList=y, s.len=25, plot=FALSE)$filteredData
		}
	
	# Producing the MDS plot
	print("Producing the MDS plot... (MDS.png)")
	png(filename='../MDS.png',width=600,height=600)
	plotMDS(y, method="bcv")
	sink('/dev/null')
	dev.off()
	sink()
} else {
	# Creating a design matrix
	if (!paired) {
		design <- model.matrix(~sampleCondition)
	} else {
		samplePatient <- as.character(read.csv('../conditions.csv',header=FALSE)$V2)
		design <- model.matrix(~samplePatient+sampleCondition)
		}

	# Estimating dataset dispersion
	print('Estimating dataset dispersion...')
	y <- estimateGLMCommonDisp(y,design)
	y <- suppressMessages(estimateGLMTrendedDisp(y,design))
	y <- estimateGLMTagwiseDisp(y,design)
	
	# Producing the MDS plot
	print("Producing the MDS plot... (MDS.png)")
	png(filename='../MDS.png',width=600,height=600)
	plotMDS(y, method="bcv")
	sink('/dev/null')
	dev.off()
	sink()

	# Fitting to the model
	print("Computing count data differential analysis...")
	fit <- glmFit(y,design)
	et <- glmLRT(fit)
	if (htsfilter) {
		print('Applying HTSFilter...')
		suppressPackageStartupMessages(library("HTSFilter"))
		et <- HTSFilter(et, DGEGLM=y, s.len=25, plot=FALSE)$filteredData
		}
	}

# Writing results in a file
print("Overview of the output:")
print(topTags(et))
print("Writing the output file...")
write.csv(as.data.frame(topTags(et,n=dim(et$table)[1])), file='file.output.csv')