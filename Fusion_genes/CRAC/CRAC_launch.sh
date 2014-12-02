#!/usr/bin/env bash

# -q sgi144g.q
#$ -q sgi512g.q
#$ -M mbahin@univ-rennes1.fr
#$ -m bea
#$ -cwd
#$ -pe mpichtc 5
# -pe make 5

# Mathieu Bahin, 10/02/14

# Script to launch CRAC.

# Template command to launch the script
# qsub ~/Fusion_genes/CRAC/CRAC_launch.sh -r $(pwd)/GIGA-HYS-2012_T01_ROTW_R1.trim.fastq.gz -s $(pwd)/GIGA-HYS-2012_T01_ROTW_R2.trim.fastq.gz -t 5

# Sourcing python environment
. /local/env/envcrac-1.5.0.sh
. /local/env/envsamtools.sh

# Getting options back
index=/home/genouest/genouest/mbahin/Fusion_genes/CRAC/Index/dogIndex
kmer=22
output=pairs.bam
stringent_chimera=FALSE
noAmbiguity=FALSE
stranded=FALSE;
detailed_sam=FALSE
while getopts "i:k:r:s:unvdt:" OPTION
do
	case $OPTION in
        i) index=$OPTARG;;
        k) kmer=$OPTARG;;
        r) reads1=$OPTARG;;
        s) reads2=$OPTARG;;
        u) stringent_chimera=TRUE;;
        n) noAmbiguity=TRUE;;
        v) stranded=TRUE;;
        d) detailed_sam=TRUE;;
        t) threads=$OPTARG;;
    esac
done

# Checking parameters
if [[ ! ("$reads1" =~ ^/) || ((-n "$reads2") && ! ("$reads2" =~ ^/)) ]]; then
	echo "The FASTQ input file(s) (options '-r' and '-s') must be specified with an absolute path. Aborting."
	exit 1
fi

# Creating a directory for the job (from the first reads file and chopping some extensions)
rep=$(basename $reads1 '.fastq.gz')
rep=${rep%.trim}
rep=${rep%_R1}
if [[ "$stringent_chimera" == TRUE && "$noAmbiguity" == TRUE ]]; then
	rep=${rep}_stringent_noAmbig_CRAC
elif [[ "$stringent_chimera" == TRUE ]]; then
	rep=${rep}_stringent_CRAC
elif [[ "$noAmbiguity" == TRUE ]]; then
	rep=${rep}_noAmbig_CRAC
else
	rep=${rep}_CRAC
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
$command | samtools view -@ $threads -Sbh - > $output

# Creating the sorted BAM and index file without secondary and supplementary alignments (used by the script to get the paired-end reads around a breakpoint)
samtools view -@ $threads -hb -F 0x800 -F 0x100 $output > pairs.clean.bam
samtools sort -@ $threads -o pairs.clean.bam sorting > pairs.clean.sort.bam
rm pairs.clean.bam
samtools index pairs.clean.sort.bam

# Checking the results (produced BAM line count almost half the input FASTQ line count)
nb_bam=$(samtools view -c $output)
nb_fastq=$(zcat $reads1 | wc -l)
check=$(echo "$nb_fastq/$nb_bam" | bc -l)
echo -e "\n\nFile sizes:\n\t$output: $nb_bam\n\t$reads1: $nb_fastq\n\tCheck test (should be greater than 1.9): $check" >> file.log
if [[ -n "$reads2" ]]; then
	if [[ $(echo "$check <= 1.9" | bc -l) -eq 1 ]]; then
		echo -e "\nWarning: Check failed !!" >> file.log
		echo "Warning: File sizes check failed !!"
	fi
else
	if [[ $(echo "$check <= 1.9" | bc -l) -eq 1 ]]; then
		echo -e "\nWarning: Check failed !!" >> file.log
		echo "Warning: File sizes check failed !!"
	fi
fi