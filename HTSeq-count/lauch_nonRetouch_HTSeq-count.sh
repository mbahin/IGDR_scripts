#$ -M mbahin@univ-rennes1.fr
#$ -m bea
#$ -cwd

# Sourcing environments
. ~/.bash_profile

#bam=/home/SAVE-OMAHA/JurasHic/FASTQ.dir/MUCOSA.dir/MUCOSA.trim_paired.dir/MUCOSA-H-PDL_BRREFAC.tophat2.0.10.unstranded.UCSC_ROUND2.dir/BAM.dir/MUCOSA-H3-PDL_BRREFAC.readname_order.bam
bam=$1
file=$(basename $bam)

samtools view -h $bam > $file.sam

htseq-count $file.sam /home/genouest/genouest/mbahin/merged.withStrand.gtf -s no > $file.count

rm $file.sam