#!/bin/bash
#
#SBATCH --job-name=grepfinch
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --output=output.%j.grep_finch
#SBATCH --time=72:00:00
#SBATCH --account=clayton
#SBATCH --partition=notchpeak
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sabrina.mcnew@gmail.com

module purge
module load samtools
module load bwa 
module load bedtools

export WORK_DIR=/scratch/general/nfs1/u0813967/insilico_dig
export PATH=/scratch/general/nfs1/u0813967/insilico_dig/bioawk/:$PATH

echo "Job started on `date`"
echo "hostname:`hostname`"

cd $WORK_DIR

fragmatic.pl -i GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna -r "CTGCA^G" -f -o finch
fragmatic.pl -i Mmel_Chromosemble_FinalVersion_Nov2018_SeqFold80.fasta -r "CTGCA^G" -f -o gamo

bioawk -c fastx 'length($seq) > 200{ print ">"$name; print $seq }' gamo.CTGCAG-CTGCAG.fasta | bioawk -c fastx 'length($seq) < 800{ print ">"$name; print $seq }' > gamo.frags.fa

bioawk -c fastx 'length($seq) > 200{ print ">"$name; print $seq }' zebrafinch.CTGCAG-CTGCAG.fasta | bioawk -c fastx 'length($seq) < 800{ print ">"$name; print $seq }' > zfinch.frags.fa

bwa index GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna
bwa index Mmel_Chromosemble_FinalVersion_Nov2018_SeqFold80.fasta

bwa mem GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna gamo.frags.fa > gamo.map.sam
bwa mem Mmel_Chromosemble_FinalVersion_Nov2018_SeqFold80.fasta zfinch.frags.fa > zfinch.map.sam

samtools flagstat gamo.map.sam
samtools flagstat zfinch.map.sam

# count up actual CpGs

cat gamo.frags.fa | grep -oh "CG" | wc -l 
cat zfinch.frags.fa | grep -oh "CG" | wc -l

cat gamo.frags.fa | grep ">" | wc -l 
cat zfinch.frags.fa | grep ">" | wc -l 

# now look at which genes were digested in the mockingbird specifically 
bwa mem Mmel_Chromosemble_FinalVersion_Nov2018_SeqFold80.fasta gamo.frags.fa | samtools view -S -b > gamo.align.gamo.bam
bedtools sort -i mmel_ns_up_onlygene.ipr.gff > mmel.sorted.gff
bedtools intersect -a gamo.align.gamo.bam -b  mmel.sorted.gff -wb -bed > intersect.gamo.bed


intersect.gamo.bed | cut -f21 | grep "Similar to" | sed 's/^.*\(Similar to .*:\).*$/\1/' | awk '{print $3}'| sed 's/://'|sort -u > insilico.genes.gamo

# how many cpg sites in each genome? 
cat Mmel_Chromosemble_FinalVersion_Nov2018_SeqFold80.fasta | grep -oh "CG" | wc -l #7 million or so
cat GCF_008822105.2_bTaeGut2.pat.W.v2_genomic.fna | grep -oh "CG" | wc -l # 7.8 million

echo "Job ended on `date`"
