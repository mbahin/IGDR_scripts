#!/bin/bash

# Mathieu Bahin, 12/06/14

# Script to get the sequences involved in a breakpoint from a fasta file produced by chimCT
# Inputs are a fasta output file from chimCT, a breakpoint coordinates (e.g. "(15:60179225,1 / X:35772912,1 # 50,62)") and an output filename.
# The output is fasta file.

# Getting options back
directory=Spanning_reads
total=FALSE
while getopts "d:f:a:b:g:h:t" OPTION
do
	case $OPTION in
    	d) directory=$OPTARG;;
    	f) fasta_file=$OPTARG;;
    	a) alignments_dir=$OPTARG;;
    	b) bkpt=$OPTARG;;
    	g) feat1=$OPTARG;;
    	h) feat2=$OPTARG;;
    	t) total=TRUE;;
    esac
done

GFFv3='/home/genouest/umr6061/recomgen/tderrien/dogomaha/DATA/canFam3/annotation/MasterAnnotation/BROADmRNA_lncRNA_antis.Ens75.gtfclean.06-02-2014.gff3'

# Checking the parameters
if [[ -z $fasta_file || -z $alignments_dir ||-z $bkpt ]]; then
	echo "Fasta file (option '-f'), alignments file directory (option '-a') and the breakpoint (option '-b') are mandatory. Please provide them. Aborting."
	exit 1
fi

# Creating a directory for the job
if [[ "$total" == FALSE ]]; then
	dir_name=$directory.bkpt.dir
else
	dir_name=$directory.full_feat.dir
fi
if [[ ! -d "$dir_name" ]]; then
	mkdir $dir_name
	cd $dir_name
else
	echo "Directory $dir_name already exists. Aborting."
	exit 1
fi

# Getting the breakpoint features
chr1=$(echo $bkpt | cut -f2 -d'(' | cut -f1 -d':')
end1=$(echo $bkpt | cut -f2 -d':' | cut -f1 -d',')
strand1=$(echo $bkpt | cut -f2 -d',' | cut -f1 -d' ')
chr2=$(echo $bkpt | cut -f2 -d':' | cut -f3 -d' ')
start2=$(echo $bkpt | cut -f3 -d':' | cut -f1 -d',')
strand2=$(echo $bkpt | cut -f3 -d',' | cut -f1 -d' ')

#####
# Getting the spanning reads
#####

echo "##### Getting the spanning reads"

# Building the spanning read file
echo "Building the spanning read file..."
bkpt_part1=$chr1'@-?1@'$end1
bkpt_part2=$chr2'@-?1@'$start2
egrep -A 1 "$bkpt_part1" $fasta_file | egrep -A 1 "$bkpt_part2" >> spanning_split_reads.fasta

# Modifying the fasta file headers
i=1
while read line; do
	if [[ $line =~ '^>' ]]; then
		echo ">Spanning_read"$i >> file.tmp; i=$(($i + 1))
	else echo $line >> file.tmp
	fi
done < spanning_split_reads.fasta
mv file.tmp spanning_split_reads.fasta

#####
# Getting the paired-end reads around the breakpoint
#####

# Getting the matched features terminals
echo "Getting the matched features terminals..."
rloc1=$(grep $feat1 $GFFv3 | awk '$3 == "mRNA"' | cut -f9 | cut -f2 -d';' | cut -f2 -d'=' | sort | uniq)
rloc2=$(grep $feat2 $GFFv3 | awk '$3 == "mRNA"' | cut -f9 | cut -f2 -d';' | cut -f2 -d'=' | sort | uniq)
feat1_beg=$(grep $rloc1 $GFFv3 | awk '$3 == "gene"'| cut -f4)
feat1_end=$(grep $rloc1 $GFFv3 | awk '$3 == "gene"'| cut -f5)
feat2_beg=$(grep $rloc2 $GFFv3 | awk '$3 == "gene"'| cut -f4)
feat2_end=$(grep $rloc2 $GFFv3 | awk '$3 == "gene"'| cut -f5)

if [[ "$total" == FALSE ]]; then
	echo "##### Getting the paired-end reads around the breakpoint"

	# Creating a BEDPE file for the extended area around the breakpoint within the two features matched
	if [[ "$strand1" == "1" ]]; then
		one_line_bedpe="$chr1\t$feat1_beg\t$end1\t$chr2\t"
	else
		one_line_bedpe="$chr1\t$end1\t$feat1_end\t$chr2\t"
	fi
	if [[ "$strand2" == "1" ]]; then
		one_line_bedpe=$one_line_bedpe"$start2\t$feat2_end\tfusion_bkpt\t0\t.\t."
	else
		one_line_bedpe=$one_line_bedpe"$feat2_beg\t$start2\tfusion_bkpt\t0\t.\t."
	fi
else
	one_line_bedpe="$chr1\t$feat1_beg\t$feat1_end\t$chr2\t$feat2_beg\t$feat2_end\tfusion_bkpt\t0\t.\t."
fi

echo -e $one_line_bedpe > file.bedpe

# Reducing the alignment file only on the concerned chromosomes
awk -v chr1=$chr1 -v chr2=$chr2 '$1 == chr1 && $4 == chr2' $alignments_dir/*.bedpe > file.primary_alignment.good_chr.bedpe
# Small for test
#awk -v chr1=$chr1 -v chr2=$chr2 '$1 == chr1 && $4 == chr2' /home/genouest/umr6061/recomgen/dog/mbahin/Fusion-genes/CRAC/Get_reads/Mini_test/small.bedpe > file.primary_alignment.good_chr.bedpe

# Doing a pairToPair between the BEDPE file about the breakpoint and the alignment file
echo "Intersecting the read pairs files..."
pairToPair -is -a file.bedpe -b file.primary_alignment.good_chr.bedpe > file.pairToPair.output

# Getting read IDs
cut -f17 file.pairToPair.output > file.readID.txt
pat=''
while read line; do
	pat=$pat"$line|"
done < file.readID.txt
pat=$(echo $pat | sed 's/^\(.*\)|$/\1/g')

# Getting the BAM lines
echo "Getting the matching BAM lines..."
samtools view $alignments_dir/pairs.primary_alignment.bam | egrep "$pat" > spanning_PE.sam
# Small for test
#samtools view /home/genouest/umr6061/recomgen/dog/mbahin/Fusion-genes/CRAC/Get_reads/Mini_test/small.bam | egrep "$pat" > spanning_PE.sam

# Outputting a fasta file
echo "Formatting the output fasta file..."
i=0
while read line; do
	echo $(echo $line | cut -f1 -d' ')"/"$(($i % 2 + 1)) >> spanning_PE_reads.fasta
	echo $(echo $line | cut -f10 -d' ') >> spanning_PE_reads.fasta
	i=$(($i + 1))
done < spanning_PE.sam

# Cleaning directory
rm file.bedpe file.pairToPair.output file.primary_alignment.good_chr.bedpe file.readID.txt spanning_PE.sam
echo "Done."