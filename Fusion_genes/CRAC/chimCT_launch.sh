#$ -M mbahin@univ-rennes1.fr
#$ -m bea
#$ -cwd

# Mathieu Bahin, 15/04/14

# Template command to launch the script
#qsub ~/Fusion_genes/CRAC/launch_chimCT.sh -i pairs.bam -n LUPA13 -r

# Configuring the environment
. /local/env/envsamtools.sh
PERL5LIB=${PERL5LIB}:$HOME/Fusion_genes/CRAC/chimCT/lib/perl5
export PERL5LIB
chimCT=/home/genouest/genouest/mbahin/Fusion_genes/CRAC/chimCT/bin/chimCT

# Getting options back
stranded=FALSE
keep_ig=FALSE
conf='/home/genouest/genouest/mbahin/Fusion_genes/CRAC/CracTools.cfg'
while getopts "g:i:n:tkc:" OPTION
do
	case $OPTION in
    	g) gff=$OPTARG;;
    	i) input=$OPTARG;;
    	n) sample_name=$OPTARG;;
    	t) stranded=TRUE;;
    	k) keep_ig=TRUE;;
    	c) conf=$OPTARG;;
    esac
done

# Checking parameters
if [[ -z $input || -z $sample_name ]]; then
	echo "Input BAM/SAM file (option '-i') and sample name (option '-n') are mandatory. Please provide them. Aborting."
	exit 1
fi
if [[ ! ("$input" =~ ^/) ]]; then
	echo "The input BAM/SAM file path must be absolute. Aborting."
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
command="$chimCT -s $input -n $sample_name --summary summary.txt --spanning-reads $sample_name --conf $conf"

if [[ "$keep_ig" == TRUE ]]; then
	command=$command' --keep-ig'
fi

# Executing command
$command > ${sample_name}.chimCT.txt
# To get the standard error output (processing steps) to debug: $command > ${sample_name}.chimCT.txt 2> stderr.log
