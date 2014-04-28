#$ -M mbahin@univ-rennes1.fr
#$ -m bea
#$ -cwd

# Sourcing environment(s)
. /local/env/envsamtools.sh
. /local/env/envpython-2.7.sh

# Defining function HTSeq-count_launch
function launch {
	command="htseq-count $1 $gtf"

	if [[ "$nonstranded" == TRUE ]]; then
		command=$command' -s no'
	elif [[ "$reverse" == TRUE ]]; then
		command=$command' -s reverse'
	fi
	if [[ "$samout_option" == TRUE ]]; then
		samout=$1'out'
		command=$command" -o $samout"
	fi
	command=$command' -q'
	
	echo -e "Command launched:\n$command" >> $log
	$command > $2
}

# Getting options back
already_sorted=FALSE
nonstranded=FALSE
samout_option=FALSE
reverse=FALSE
while getopts "b:g:anro" OPTION
do
	case $OPTION in
		b) bam=$OPTARG;;
		g) gtf=$OPTARG;;
		a) already_sorted=TRUE;;
		n) nonstranded=TRUE;;
		r) reverse=TRUE;;
		o) samout_option=TRUE;;
	esac
done

# Checking parameters validity
if [[ -z $bam || -z $gtf ]]; then
	echo "BAM (option '-b') and GTF (option '-g') files are mandatory. Please provide them. Aborting."
	exit 1
fi
if [[ ! ("$gtf" =~ ^/) ]]; then
	echo "The input GTF file path must be absolute. Aborting."
	exit 1
fi
if [[ ("$reverse" == TRUE) && ("$nonstranded" == TRUE) ]]; then
	echo 'The incompatible non-stranded and reverse modes were chosen. Aborting.'
	exit 1
fi

# Creating a directory for the job
rep=$(basename $bam .readname_order.bam)
if [[ ! -d "$rep.dir" ]]; then
	mkdir $rep.dir
	cd $rep.dir
else
	echo "Directory $rep already exists. Aborting."
	exit
fi

# Printing script metadata
log=file.log
echo -e "----- Metadata -----\n" > $log
echo -e "BAM file processed:\n$bam\n" >> $log
echo -e "GTF file used:\n$gtf\n" >> $log

# Changing the BAM input file into a SAM one
sam=file.sorted.sam
header=file.header
samtools view -H $bam > $header
if [[ "$already_sorted" == TRUE ]]; then
	samtools view $bam > $sam
else
	samtools sort -on $bam sorting | samtools view - > $sam
fi

# Splitting the input dataset into a paired-end file and a non-paired-end file
awk '!and($2,0x1)' $sam > file.tmp
if [[ -s file.tmp ]]; then
	sam_SE=file.sorted.SE.sam
	cat $header file.tmp > $sam_SE
fi
rm file.tmp
awk 'and($2,0x1)' $sam > file.tmp
if [[ -s file.tmp ]]; then
	sam_PE=file.sorted.PE.sam
	cat $header file.tmp > $sam_PE
fi

# Launching HTSeq-count
output=file.$rep.count
if [[ (-e $sam_SE) && (-e $sam_PE) ]]; then
	count_SE=file.SE.count
	count_PE=file.PE.count
	
	# Executing HTSeq-count
	echo -e "Paired and non-paired reads processed.\n" >> $log
	launch $sam_SE $count_SE
	launch $sam_PE $count_PE
	
	# Merging the results
	paste $count_SE $count_PE | awk 'BEGIN{OFS="\t"}{print $1,$2+$4}' > $output
else
	if [[ -e $sam_SE ]]; then
		sam=$sam_SE
	else
		sam=$sam_PE
	fi
	# Executing HTSeq-count
	echo -e "Homogeneous reads (paired or non-paired) processed.\n" >> $log
	#HTSeq-count_launch $sam $output
	launch $sam $output
fi

# Cleaning files
#rm $sam $sam_SE $sam_PE $header file.tmp $count_SE $count_PE