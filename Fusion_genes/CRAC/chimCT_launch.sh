#$ -M mbahin@univ-rennes1.fr
#$ -m bea
#$ -cwd

# Mathieu Bahin, 15/04/14

# Template command to launch the script
#qsub ~/Fusion_genes/CRAC/launch_chimCT.sh -c <CRAC_output_path> -n <sample_name>

# Configuring the environment
. /local/env/envsamtools.sh
PERL5LIB=${PERL5LIB}:$HOME/Fusion_genes/CRAC/chimCT/lib/perl5
export PERL5LIB
chimCT=/home/genouest/genouest/mbahin/Fusion_genes/CRAC/chimCT/bin/chimCT

# Getting options back
stranded=FALSE
keep_ig=FALSE
single_end=FALSE
conf='/home/genouest/genouest/mbahin/Fusion_genes/CRAC/CracTools.cfg'
while getopts "c:n:sktg:c:" OPTION
do
	case $OPTION in
    	c) crac_dir=$OPTARG;;
    	n) sample_name=$OPTARG;;
    	s) stranded=TRUE;;
    	k) keep_ig=TRUE;;
    	t) single_end=TRUE;;
    	g) gff=$OPTARG;;
    	c) conf=$OPTARG;;
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
directory=${sample_name}_chimCT.dir
if [[ -d "$directory" ]]; then
	echo "Directory $directory already exists. Aborting."
	exit 1
else
	mkdir $directory
	cd $directory
fi

# Command building
command="$chimCT -s $crac_dir/pairs.bam -n $sample_name --summary summary.txt --spanning-reads $sample_name --conf $conf"

if [[ "$keep_ig" == TRUE ]]; then
	command=$command' --keep-ig'
fi

# Executing chimCT command
$command > ${sample_name}.chimCT.txt
# To get the standard error output (processing steps) to debug: $command > ${sample_name}.chimCT.txt 2> stderr.log

# Processing chimCT output
if [[ "$single_end" == FALSE ]]; then
	~/process_chimCT_output.py -n $sample_name -c $crac_dir
else
	~/process_chimCT_output.py -n $sample_name -c $crac_dir -s
fi
