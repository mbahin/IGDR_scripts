#! /bin/bash

# Mathieu Bahin, 31/03/13

# Sourcing environment(s)
. /local/env/envR-3.1.0.sh

# Config
#DESeq2_launch='/home/genouest/umr6061/recomgen/dog/mbahin/DE/scripts/DESeq2_launch.r'
#edgeR_launch='/home/genouest/umr6061/recomgen/dog/mbahin/DE/scripts/edgeR_launch.r'
#merge_analyses_script='/home/genouest/umr6061/recomgen/dog/mbahin/DE/scripts/merge_analyses.py'
#create_Venn='/home/genouest/umr6061/recomgen/dog/mbahin/DE/scripts/create_Venn.r'
DESeq2_launch='/home/genouest/genouest/mbahin/DE/DESeq2_launch.r'
edgeR_launch='/home/genouest/genouest/mbahin/DE/edgeR_launch.r'
merge_analyses_script='/home/genouest/genouest/mbahin/DE/merge_analyses.py'
create_Venn='/home/genouest/genouest/mbahin/DE/create_Venn.r'

# Getting options back
mode='de'
directory='multiDE'
paired_design=FALSE
glm=FALSE
htsfilter=FALSE
threshold=0.05
while getopts "m:d:c:t:pgh" OPTION
do
	case $OPTION in
    	m) mode=$OPTARG;;
    	d) directory=$OPTARG;;
    	c) metadata=$OPTARG;;
    	t) threshold=$OPTARG;;
    	p) paired_design=TRUE;;
    	g) glm=TRUE;;
    	h) htsfilter=TRUE;;
    esac
done

# Checking parameters
if [[ -z $metadata ]]; then
	echo "Metadata file (option '-c') is mandatory. Please provide it. Aborting."
	exit 1
fi
if [[ ! ("$mode" != 'de' || "$mode" != 'ed' || "$mode" != 'e' || "$mode" != 'd') ]]; then
	echo $mode
	echo "The mode should be 'de' (or 'ed'), 'd' or 'e'. Aborting."
	exit 1
fi
if [[ ! ("$metadata" =~ ^/) ]]; then
	echo "The metadata file (options '-c') must specified with absolute path. Aborting."
	exit 1
fi
if [[ ! ("$threshold" =~ 0.[0-9]*) ]]; then
	echo "The threshold must be numeric. Aborting."
	exit 1
fi
if [[ ("$paired_design" == TRUE) && ("$mode" =~ 'e') && ("$glm" == FALSE) ]]; then
	echo "Paired mode and edgeR analysis were chosen but GLM approach wasn't. Aborting."
	exit 1
fi
if [[ -n $(awk -F, 'NR>1 && NF<2' $metadata) ]]; then
	echo "At least one line in the factor file, beside the header, doesn't have the 2 minimum fields (filename and condition). Aborting."
fi
if [[ (-n $(awk -F, 'NR>1 && NF<3' $metadata)) && ("$paired_design" == TRUE) ]]; then
	echo "Paired mode chosen but pairing information are not properly found in the factor file. Aborting."
	exit 1
fi
if [[ ("$paired_design" == TRUE) && ("$htsfilter" == TRUE) ]]; then
	echo "Warning: Paired mode and htsfilter chosen. HTSFilter need at least two replicates, it won't be used."
	htsfilter=FALSE
fi

if [[ -d "$directory" ]]; then
	echo "Directory $directory already exists. Aborting."
	exit 1
else
	mkdir $directory
	cd $directory
	mkdir Counts
	for file in $(sed 1d $metadata | cut -f1 -d,); do
		cp $(head -1 $metadata | cut -f1 -d,)/$file Counts/
	done
	cp $metadata metadata.csv
fi

# Resizing the count filenames (removing the longest common prefix and suffix)
cd Counts
prefix_size=0
prefix=''
suffix_size=0
suffix=''
first_file=$(ls | head -1)

# Determining the longest common prefix
while true; do
    prefix_size=$((prefix_size+1))
    prefix=${first_file:0:prefix_size}
    for filename in $(ls); do
    	if [[ ${filename:0:$prefix_size} != "$prefix" ]]; then
    		prefix_size=$((prefix_size-1))
    		prefix=${first_file:0:prefix_size}
    		break 2
    	fi
    done
done

# Determining the longest common suffix
while true; do
    suffix_size=$((suffix_size+1))
    suffix=${first_file:${#first_file}-suffix_size:${#first_file}}
    for filename in $(ls); do
    	if [[ ${filename:${#filename}-suffix_size:${#filename}} != "$suffix" ]]; then
    		suffix_size=$((suffix_size-1))
    		suffix=${first_file:${#first_file}-suffix_size:${#first_file}}
    		break 2
    	fi
    done
done

# Modifying count filenames
for filename in $(ls); do
	if [[ $filename != ${filename#$prefix} ]]; then
		mv $filename ${filename#$prefix}
		sed -i "s/$filename/${filename#$prefix}/g" ../metadata.csv
	fi
done

for filename in $(ls); do
	if [[ $filename != ${filename#$suffix} ]]; then
		mv $filename ${filename%$suffix}
		sed -i "s/$filename/${filename%$suffix}/g" ../metadata.csv
	fi
done
cd ..

# Launching the DE analyses
if [[ "$mode" =~ 'd' ]]; then
	DESeq2_mode='DESeq2'
	if [[ "$paired_design" == TRUE ]]; then
		DESeq2_mode=$DESeq2_mode'_paired'
	fi
	if [[ "$htsfilter" == TRUE ]]; then
		DESeq2_mode=$DESeq2_mode'_HTSF'
	fi
	echo "====== Launching a $DESeq2_mode analysis on the input data..."
	mkdir $DESeq2_mode
	cd $DESeq2_mode
	Rscript $DESeq2_launch "../Counts" $htsfilter $paired_design
	cd ..
fi
if [[ "$mode" =~ 'e' ]]; then
	edgeR_mode='edgeR'
	if [[ "$glm" == TRUE ]]; then
		edgeR_mode=$edgeR_mode'_GLM'
	fi
	if [[ "$paired_design" == TRUE ]]; then
		edgeR_mode=$edgeR_mode'_paired'
	fi
	if [[ "$htsfilter" == TRUE ]]; then
		edgeR_mode=$edgeR_mode'_HTSF'
	fi
	echo "===== Launching a $edgeR_mode analysis on the input data..."
	mkdir $edgeR_mode
	cd $edgeR_mode
	Rscript $edgeR_launch "../Counts" $htsfilter $glm $paired_design
	cd ..
fi

# Intersecting the results if both DESeq2 and edgeR analyses are requested
if [[ ${#mode} == 2 ]]; then
	echo "===== Intersecting DESeq2 and edgeR analysis results..."
	/local/python/2.7/bin/python $merge_analyses_script -t $threshold -d $DESeq2_mode/file.DESeq2.output.csv -e $edgeR_mode/file.edgeR.output.csv

	# Creating a Venn diagram
	Rscript $create_Venn 'file.venn.csv'
fi
echo "Done."

# Cleaning
rm -r Counts metadata.csv
