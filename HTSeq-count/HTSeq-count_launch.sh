#$ -M mbahin@univ-rennes1.fr
#$ -m bea
#$ -cwd

# Template command to launch the script
#qsub ~/HTSeq-count/HTSeq-count_launch.sh -i file.bam -g file.gtf -afso

# Sourcing environment(s)
. /local/env/envsamtools.sh
. /local/env/envpython-2.7.sh

# Defining the function launching HTSeq-count
function launch {
    if [[ "$transcript_mode" == TRUE ]]; then
    	command="htseq-count -i transcript_id $1 $gtf"
    else
    	command="htseq-count $1 $gtf"
    fi

	# TEMPORARY
	#if [[ "$order_pos" == TRUE ]]; then
	#	command=$command' -r pos'
	#fi
	if [[ "$nonstranded" == TRUE ]]; then
		command=$command' -s no'
	elif [[ "$stranded_reverse" == TRUE ]]; then
		command=$command' -s reverse'
	fi
	if [[ "$samout_option" == TRUE ]]; then
		samout=file.$rep.samout
		command=$command" -o $samout"
	fi
	command=$command' -q'
	
	echo "$command" >> $log
	$command > $2
}

# Getting options back
already_sorted=FALSE # TEMPORARY
format_bam=FALSE
#order_pos=FALSE # TEMPORARY
transcript_mode=FALSE
nonstranded=FALSE
stranded_reverse=FALSE
samout_option=FALSE
while getopts "i:g:aftnso" OPTION
do
	case $OPTION in
		i) input=$OPTARG;;
		g) gtf=$OPTARG;;
		a) already_sorted=TRUE;;	
		f) format_bam=TRUE;;
		#r) order_pos=TRUE;; # TEMPORARY
		t) transcript_mode=TRUE;;
		n) nonstranded=TRUE;;
		s) stranded_reverse=TRUE;;
		o) samout_option=TRUE;;
	esac
done

# Checking parameters validity
if [[ -z $input || -z $gtf ]]; then
	echo "Input (option '-i') and GTF (option '-g') files are mandatory. Please provide them. Aborting."
	exit 1
fi
if [[ ! ("$gtf" =~ ^/) || ! ("$input" =~ ^/) ]]; then
	echo "The input GTF and BAM/SAM file paths must be absolute. Aborting."
	exit 1
fi
if [[ ("$stranded_reverse" == TRUE) && ("$nonstranded" == TRUE) ]]; then
	echo 'The incompatible non-stranded and stranded reverse modes were chosen. Aborting.'
	exit 1
fi

# Creating a directory for the job
rep=$(basename $(echo $input | sed 's/readname_order.//g') .bam)
if [[ ! -d "${rep}_count.dir" ]]; then
	mkdir ${rep}_count.dir
	cd ${rep}_count.dir
else
	echo "Directory ${rep}_count already exists. Aborting."
	exit
fi

# Printing script metadata
log=file.log
echo -e "----- Metadata -----\n" > $log
echo -e "BAM file processed:\n$input\n" >> $log
echo -e "GTF file used:\n$gtf\n" >> $log

# TEMPORARY: sorting file by name, otherwise, "Maximum alignment buffer size exceeded while pairing SAM alignments" error may occur.
if [[ "$format_bam" == TRUE ]]; then
	if [[ "$already_sorted" == FALSE ]]; then
		samtools sort -on $input sorting > file.sorted
	fi
else
	if [[ "$already_sorted" == FALSE ]]; then
		samtools sort -on $input sorting | samtools view - > file.sorted
	fi
fi

# Splitting the input dataset into a paired-end reads file and a single-end reads file and changing into SAM files if necessary
if [[ "$already_sorted" == FALSE ]]; then
	input=file.sorted
fi
if [[ "$format_bam" == TRUE ]]; then
	samtools view -h -F 0x1 $input > file.tmp
	if [[ $(awk '! /^@/' file.tmp | wc -l) -ne 0 ]]; then
		input_SE=file.SE.bam
		mv file.tmp $input_SE
	fi
	samtools view -h -f 0x1 $input > file.tmp
	if [[ $(awk '! /^@/' file.tmp | wc -l) -ne 0 ]]; then
		input_PE=file.PE.bam
		mv file.tmp $input_PE
	fi
else
	samtools view -Sh -F 0x1 $input > file.tmp
	if [[ $(awk '! /^@/' file.tmp | wc -l) -ne 0 ]]; then
		input_SE=file.SE.sam
		mv file.tmp $input_SE
	fi
	samtools view -Sh -f 0x1 $input > file.tmp
	if [[ $(awk '! /^@/' file.tmp | wc -l) -ne 0 ]]; then
		input_PE=file.PE.sam
		mv file.tmp $input_PE
	fi
fi

# Launching HTSeq-count
output=file.$rep.count
if [[ (-e $input_SE) && (-e $input_PE) ]]; then
	count_SE=file.SE.count
	count_PE=file.PE.count
	
	# Executing HTSeq-count
	echo -e "Paired and non-paired reads processed.\n" >> $log
	echo 'Command(s) launched:' >> $log
	launch $input_SE $count_SE
	launch $input_PE $count_PE
	
	# Merging the results
	paste $count_SE $count_PE | awk 'BEGIN{OFS="\t"}{print $1,$2+$4}' > $output
else
	if [[ -e $input_SE ]]; then
		input=$input_SE
	else
		input=$input_PE
	fi
	# Executing HTSeq-count
	echo -e "Homogeneous reads (paired or non-paired) processed.\n" >> $log
	launch $input $output
fi

# Cleaning files
rm -f file.sorted $input_SE $input_PE file.tmp $count_SE $count_PE
