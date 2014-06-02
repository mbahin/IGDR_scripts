#$ -S /bin/bash
#$ -cwd

# Template command to launch the script in a loop
#for i in $(ls <BAM files directory>/*.bam); do qsub <script> -b $i -g <file.gtf> -c -o -s; done

# Sourcing python environment
. /local/env/envpython-2.7.sh
. /local/env/envsamtools.sh

# Getting options back
while getopts "b:g:soc" OPTION
do
        case $OPTION in
                b) bam=$OPTARG;;
                g) gtf=$OPTARG;;
                s) stranded="no";; # To launch with '-s no'
                o) samout_option="yes";;
                c) custom_stat_option="yes";; # To launch with custom_stat

        esac
done

if [[ -z $bam || -z $gtf ]]; then
	echo "BAM (option '-b') and GTF (option '-g') files are mandatory. Please provide them. Aborting."
	exit 1
fi

# Defining output files
output=file.output.count
bamout=file.output.bam
log=file.log

# Creating a working directory
#rep=$(basename $1 .PE.picard.reorder.bam).dir
rep=$(basename $bam .bam).dir
mkdir $rep
cd $rep

# Working file
sam=file.sorted.sam
header=file.header
samout_head=file.feature.h.sam

# Printing script metadata
echo -e "----- Metadata -----\n" > $log
echo -e "BAM file processed:\n$bam\n" >> $log
echo -e "GTF file used:\n$gtf\n" >> $log

# Building command-line to count reads
command="/home/genouest/genouest/mbahin/HTSeq/HTSeq_retouch/HTSeq-count.retouch.py $sam $gtf -q"
if [ "$stranded" == "no" ]; then
	command=$command" -s no"
fi
if [ "$samout_option" == "yes" ]; then
	samout=file.feature.sam
	command=$command" -o $samout"
fi
if [ "$custom_stat_option" == "yes" ]; then
	stat=file.custom_stat.txt
	command=$command" -c $stat"
fi
#command=$command" > $output 2> err.txt"
echo -e "Command-line:\n$command" >> $log

# Sorting BAM file and converting to SAM
samtools sort -on $bam sorting | samtools view -h - > $sam 2> /dev/null

# Counting reads
$command > $output 2> err.txt

# Adding header to the annotated SAM file
samtools view -SH $sam > $header 2> /dev/null
cat $header $samout > $samout_head

# Converting SAM to BAM
samtools view -bSh $samout_head > $bamout 2> /dev/null

# Cleaning files
rm $sam $header $samout $samout_head err.txt

# Going back to the original directory
cd - > /dev/null