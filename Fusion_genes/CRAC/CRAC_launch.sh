#!/usr/bin/env bash

#$ -q sgi144g.q
#$ -M mbahin@univ-rennes1.fr
#$ -m bea
#$ -cwd
#$ -pe mpichtc 5

# Mathieu Bahin, 10/02/14

# Script to launch CRAC.

# Template command to launch the script
# Paired:
# qsub ~/Fusion_genes/CRAC/CRAC_script.sh -r /home/SAVE-OMAHA/JurasHic/FASTQ.dir/LUPA.dir/LUPA.trim_paired.dir/LUPA13_R1.paired.trim.fastq.gz -s /home/SAVE-OMAHA/JurasHic/FASTQ.dir/LUPA.dir/LUPA.trim_paired.dir/LUPA13_R2.paired.trim.fastq.gz -t 5
# To work on multi-threads
# -pe mpichtc 5

# Sourcing python environment
. /local/env/envcrac-1.5.0.sh
. /local/env/envsamtools.sh

# Getting options back
#directory='CRAC_results'
index=/home/genouest/genouest/mbahin/Fusion_genes/CRAC/Index/dogIndex
kmer=22
#paired_end=TRUE
output=pairs.bam
stringent_chimera=FALSE
noAmbiguity=FALSE
stranded=FALSE;
detailed_sam=FALSE
#while getopts "c:i:k:r:s:po:uvdt:" OPTION
while getopts "i:k:r:s:unvdt:" OPTION
do
	case $OPTION in
		#d) directory=$OPTARG;;
        i) index=$OPTARG;;
        k) kmer=$OPTARG;;
        r) reads1=$OPTARG;;
        s) reads2=$OPTARG;;
        #p) paired_end=FALSE;;
        u) stringent_chimera=TRUE;;
        n) noAmbiguity=TRUE;;
        v) stranded=TRUE;;
        d) detailed_sam=TRUE;;
        t) threads=$OPTARG;;
    esac
done

# Checking parameters
#if [[ "$paired_end" == FALSE && -n $reads2 ]]; then
#	echo "Single-end mode chosen (option '-e') but second read file given (option '-s'). Aborting."
#	exit 1
#fi
#if [[ "$paired_end" == TRUE && ( -z $reads1 || -z $reads2 ) ]]; then
#	echo "Paired-end mode chosen (option '-e') but first or second read file missing. Aborting."
#	exit 1
#fi

# Creating a directory for the job
if [[ "$stringent_chimera" == TRUE ]]; then
	rep=$(basename $reads1 '_R1.trim.fastq.gz')_stringent_CRAC
else
	rep=$(basename $reads1 '_R1.trim.fastq.gz')_CRAC
fi
if [[ ! -d "$rep.dir" ]]; then
	mkdir $rep.dir
	cd $rep.dir
else
	echo "Directory $rep already exists. Aborting."
	exit 1
fi

# Printing script metadata
log=file.log
echo -e "Date: "$(date)"\n" > file.log
echo -e "Index:\n$index" >> file.log
#if [[ "$paired_end" == TRUE ]]; then
#	echo -e "\nRead file (R1):\n$reads1\n" >> file.log
#	echo -e "Read file (R2):\n$reads2" >> file.log
#fi
echo -e "\nRead file (R1):\n$reads1" >> file.log
if [[ -n "$reads2" ]]; then
	echo -e "Read file (R2):\n$reads2" >> file.log
fi

# Building and launching the command
command="crac -i $index -k $kmer -r $reads1"
#if [[ "$paired_end" == TRUE ]]; then
#	command=$command" $reads2 --paired-end-chimera pe.chimera.output"
#fi
if [[ -n "$reads2" ]]; then
	command=$command" $reads2"
fi
if [[ "$stringent_chimera" == TRUE ]]; then
	command=$command' --stringent-chimera'
fi
if [[ "$noAmbiguity" == TRUE ]]; then
	command=$command' --no-ambiguity'
fi
if [[ "$stranded" == TRUE ]]; then
	command=$command' --stranded'
fi
if [[ "$detailed_sam" == TRUE ]]; then
	command=$command' --detailed-sam'
fi
if [[ -n "$threads" ]]; then
	command=$command" --nb-threads $threads"
fi
command=$command" -o- --summary summary.output"
echo -e "\nOriginal command line:\n$command" >> file.log
$command | samtools view -Sbh -  > $output

# Creating the BAM and BEDPE file with only the primary alignments (used by the script to get the paired-end reads around a breakpoint)
samtools view -hb -F 0x800 $output > pairs.primary_alignment.bam
bamToBed -bedpe -i pairs.primary_alignment.bam > pairs.primary_alignment.bedpe