#!/usr/bin/env bash

#$ -q sgi144g.q
#$ -M mbahin@univ-rennes1.fr
#$ -m bea
#$ -cwd

# Mathieu Bahin, 10/02/14

# Script to launch CRAC.

# Template command to launch the script
# Paired:
# ~/Fusion_genes/CRAC/CRAC_script.sh -r /omaha-beach/JurasHic/FASTQ.dir/LUPA.dir/LUPA.trim_paired.dir/LUPA13_R1.paired.trim.fastq.gz -s /omaha-beach/JurasHic/FASTQ.dir/LUPA.dir/LUPA.trim_paired.dir/LUPA13_R2.paired.trim.fastq.gz -t 5
# To work on multi-threads
# -pe mpichtc 5

# Sourcing python environment
. /local/env/envcrac-1.3.2.sh
. /local/env/envsamtools.sh

# Getting options back
index=/omaha-beach/mbahin/CRAC/Index/dogIndex
kmer=22
paired_end=TRUE
output=pairs.sam
stringent_chimera=FALSE
stranded=FALSE;
detailed_sam=FALSE
while getopts "i:k:r:s:po:uvdt:" OPTION
do
        case $OPTION in
                i) index=$OPTARG;;
                k) kmer=$OPTARG;;
                r) reads1=$OPTARG;;
                s) reads2=$OPTARG;;
                p) paired_end=FALSE;;
                o) output=$OPTARG;;
                u) stringent_chimera=TRUE;;
                v) stranded=TRUE;;
                d) detailed_sam=TRUE;;
                t) threads=$OPTARG;;
        esac
done

# Checking parameters
if [[ "$paired_end" == FALSE && -n $reads2 ]]; then
	echo "Single-end mode chosen (option '-e') but second read file given (option '-s'). Aborting."
	exit 1
fi
if [[ "$paired_end" == TRUE && ( -z $reads1 || -z $reads2 ) ]]; then
	echo "Paired-end mode chosen (option '-e') but first or second read file missing. Aborting."
	exit 1
fi

# Creating a working directory
name=$(basename $reads1 '_R1.paired.trim.fastq.gz')
rep=$name.dir
if [[ -d $rep ]]; then
	echo -e "Failed to create a directory named '$rep' because there is already one existing. Aborting."
	exit 1
fi
mkdir $rep
cd $rep

# Printing script metadata
log=file.log
echo -e "Date: "$(date)"\n" > file.log
echo -e "Original command line:\nCRAC_script.sh -i $index -k $kmer -r $reads1 $reads2 -o $output --summary summary.output --chimera chimera.output\n" >> file.log
echo -e "Index:\n$index" >> file.log
if [[ "$paired_end" == TRUE ]]; then
	echo -e "\nRead file (R1):\n$reads1\n" >> file.log
	echo -e "Read file (R2):\n$reads2" >> file.log
fi

# Building and launching the command
command="crac -i $index -k $kmer -r $reads1"
if [[ "$paired_end" == TRUE ]]; then
	command=$command" $reads2 --paired-end-chimera pe.chimera.output"
fi
if [[ "$stringent" == TRUE ]]; then
	command=$command' --stringent-chimera'
fi
if [[ "$stranded" == TRUE ]]; then
	command=$command' --stranded'
fi
if [[ "$detailed_sam" == TRUE ]]; then
	command=$command' --detailed-sam'
fi
if [[ -n $threads ]]; then
	command=$command" --nb-threads $threads"
fi
command=$command" -o $output --summary summary.output --chimera chimera.output"
$command

# Creating a BAM output if necessary
samtools view -bSh $output > $(basename $output .sam).bam 2> /dev/null
rm $output