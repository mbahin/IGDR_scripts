#$ -M mbahin@univ-rennes1.fr
#$ -m bea
#$ -cwd

# Mathieu Bahin, 15/04/14

# Template command to launch the script
#qsub ~/Fusion_genes/CRAC/launch_chimCT.sh -c <CRAC_output_path> -n <sample_name>

# Configuring the environment
. /local/env/envsamtools.sh
PERL5LIB=$PERL5LIB:$HOME/Fusion_genes/CRAC/chimCT/lib/perl5
export PERL5LIB
chimCT=/home/genouest/genouest/mbahin/Fusion_genes/CRAC/chimCT/bin/chimCT

# Getting options back
stranded=FALSE
single_end=FALSE
conf='/home/genouest/genouest/mbahin/Fusion_genes/CRAC/CracTools.cfg'
while getopts "d:c:n:stg:" OPTION
do
	case $OPTION in
		d) output_dir=$OPTARG;;
		c) crac_dir=$OPTARG;;
		n) sample_name=$OPTARG;;
		s) stranded=TRUE;;
		t) single_end=TRUE;;
		g) gff=$OPTARG;;
	esac
done

# Checking parameters
if [[ -z $crac_dir || -z $sample_name ]]; then
	echo "Input CRAC output directory (option '-c') and sample name (option '-n') are mandatory. Please provide them. Aborting."
	exit 1
fi
if [[ ! ("$crac_dir" =~ ^/) ]]; then
	echo "The input CRAC output directory path must be absolute. Aborting."
	exit 1
fi

# Creating the output directory
if [[ -n $output_dir ]]; then
	directory=$output_dir
else
	directory=${sample_name}_chimCT.dir
fi
if [[ -d "$directory" ]]; then
	echo "Directory $directory already exists. Aborting."
	exit 1
else
	mkdir $directory
	cd $directory
fi

# Printing script metadata
log=file.log
echo -e "Date: "$(date)"\n" > $log
        # If GFF specified, replacing it in the config file
if [[ -n $gff ]]; then
        head -6 $conf > CracTools.cfg
        conf=CracTools.cfg
        echo "ANNOTATION_GFF '$gff'" >> CracTools.cfg
fi
echo -e "GFF file: "$(grep ^ANNOTATION_GFF $conf | sed "s/ANNOTATION_GFF '\(.*\)'/\1/g") >> $log

# Command building
command="$chimCT -s $crac_dir/pairs.bam -n $sample_name --summary summary.txt --spanning-reads $sample_name --conf $conf"
echo -e "\nOriginal command line:\n$command" >> $log

# Executing chimCT command
$command > ${sample_name}.chimCT.txt
# To get the standard error output (processing steps) to debug: $command > ${sample_name}.chimCT.txt 2> stderr.log

# Processing chimCT output
if [[ "$single_end" == FALSE ]]; then
	~mbahin/Fusion_genes/CRAC/process_chimCT_output.py -n $sample_name -c $crac_dir
else
	~mbahin/Fusion_genes/CRAC/process_chimCT_output.py -n $sample_name -c $crac_dir -s
fi
