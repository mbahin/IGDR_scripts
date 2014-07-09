#$ -M mbahin@univ-rennes1.fr
#$ -m bea
#$ -cwd

# Mathieu Bahin, 15/04/14

# Template command to launch the script
#qsub ~/Fusion_genes/CRAC/launch_chimCT.sh -g /home/genouest/umr6061/recomgen/tderrien/DATA/canFam3/annotation/MasterAnnotation/BROADmRNA_lncRNA_antis.Ens75.gtfclean.05-06-2014.gff3 -s pairs.sam -n LUPA13 -r

# Config
. ~/.bash_profile
chimCT=/home/genouest/genouest/mbahin/Fusion_genes/CRAC/chimCT/bin/chimCT

# Getting options back
#directory='chimCT_results'
stranded=FALSE
keep_ig=FALSE
spanning_reads=FALSE
conf='/home/genouest/genouest/mbahin/Fusion_genes/CRAC/CracTools.cfg'
while getopts "g:s:n:tkrc:" OPTION
do
	case $OPTION in
    	#d) directory=$OPTARG;;
    	g) gff=$OPTARG;;
    	s) sam=$OPTARG;;
    	n) sample_name=$OPTARG;;
    	t) stranded=TRUE;;
    	k) keep_ig=TRUE;;
    	r) spanning_reads=TRUE;;
    	c) conf=$OPTARG;;
    esac
done

# Checking parameters
if [[ -z $sam || -z $sample_name ]]; then
	echo "SAM file (option '-s') and sample name (option '-n') are mandatory. Please provide them. Aborting."
	exit 1
fi
if [[ ! ("$sam" =~ ^/) ]]; then
	echo "The input SAM file path must be absolute. Aborting."
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
command="$chimCT -s $sam -n $sample_name --summary summary.txt --conf $conf"

if [[ "$spanning_reads" == TRUE ]]; then
	command=$command' --spanning-reads '$sample_name
fi
if [[ "$keep_ig" == TRUE ]]; then
	command=$command' --keep-ig'
fi

# Executing command
$command > file.chim.txt 2> stderr.log