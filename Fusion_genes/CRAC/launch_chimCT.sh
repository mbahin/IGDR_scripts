#$ -M mbahin@univ-rennes1.fr
#$ -m bea
#$ -cwd

# Mathieu Bahin, 15/04/14

# Script to launch a chimCT study

# Config
. ~/.bash_profile
chimCT=/home/genouest/genouest/mbahin/Fusion_genes/CRAC/chimCT/bin/chimCT

#sam=/home/genouest/umr6061/recomgen/dog/mbahin/Fusion-genes/chimCT/file.withXP.sam
#gff3=/home/genouest/umr6061/recomgen/dog/mbahin/chimCT/Res/file.gff3
#gff3_withID=/home/genouest/umr6061/recomgen/dog/mbahin/chimCT/Res/file.withID.gff3

# Getting options back
directory='chimCT_results'
stranded=FALSE
keep_ig=FALSE
spanning_reads=FALSE
while getopts "d:g:s:n:tkr" OPTION
do
	case $OPTION in
    	d) directory=$OPTARG;;
    	g) gff=$OPTARG;;
    	s) sam=$OPTARG;;
    	n) sample_name=$OPTARG;;
    	t) stranded=TRUE;;
    	k) keep_ig=TRUE;;
    	r) spanning_reads=TRUE;;
    esac
done

# Checking parameters
if [[ -z $sam || -z $sample_name ]]; then
	echo "Mode (option '-m'), counts file (option '-c') and factors file (option '-f') are mandatory. Please provide them. Aborting."
	# TODO
	exit 1
fi
if [[ -d "$directory" ]]; then
	echo "Directory $directory already exists. Aborting."
	exit 1
else
	mkdir $directory
	cd $directory
fi

# Command building
command="$chimCT -s $sam -n $sample_name --summary summary.txt"

if [[ -n $spanning_reads ]]; then
	command=$command' --spanning-reads spanR'
fi

# Executing command
#chimCT -g $gff3 -s $sam -n LUPA11_chimCT_test --primers primers.txt --summary summary.txt --spanning-reads spanR
#chimCT -g $gff3 -s $sam -n LUPA11_chimCT_test --primers --summary summary.txt --spanning-reads spanR
#chimCT -g $gff3 -s $sam -n LUPA11_chimCT_test --summary summary.txt --spanning-reads spanR
#chimCT -g $gff3_withID -s $sam -n LUPA11_chimCT_test --primers --summary summary.txt --spanning-reads spanR
#chimCT -g $gff3_withID -s $sam -n LUPA11_chimCT_test --summary summary.txt --spanning-reads spanR
#chimCT -s $sam -n LUPA11_chimCT_test --keep-ig --summary summary.txt --spanning-reads spanR
log=file.log
echo -e "----- Inputs -----\n" > $log
echo -e "SAM file processed:\n$sam\n" >> $log
echo -e "GFF file used:" >> $log
if [[ -n $gff ]]; then
	echo -e "$gff\n" >> $log
else
	echo -e "GFF file from config file used.\n" >> $log
fi
echo -e "Command launched:\n$command" >> $log
$command > file.chim.txt 2> stderr.log

# Copy the config file in the directory and delete it