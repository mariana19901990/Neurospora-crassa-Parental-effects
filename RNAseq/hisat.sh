#!/bin/bash
#SBATCH --job-name=hisat_align
#SBATCH --account=project_2000350
#SBATCH --time=00:10:00
#SBATCH --ntasks=1 
#SBATCH --mem-per-cpu=1G

module load biokit
module load parallel/20200122
module load bioconda

#Stop at any error
set -ue
#Alignent base RNA-seq analysis
S1PATH=/path_to_rawdata/raw_data
S2PATH=/path_to_aligned_files/rna-seq/bam
IDX=/path_to_indexed_genome/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata_ERCC.fasta
ids=/path_to_sample_ids/ids
#Use hisat to align reads to the reference, first build the reference and index it
hisat2-build $IDX $IDX
# Index the reference genome with samtools.
samtools faidx $IDX

# Align the FASTQ files to the reference genome.
cat $ids | parallel --compress "hisat2 --max-intronlen 2500 -x ${IDX} -1 $S1PATH/{}/{}_1.fq.gz -2 $S1PATH/{}/{}_2.fq.gz | samtools sort > $S2PATH/{}.bam"

# Index each BAM file.
cat $ids | parallel --compress "samtools index $S2PATH/{}.bam"

# To see the alignments in IGV is necessart to turn them into bigWig format:
# Turn each BAM file into bedGraph coverage. The files will have the .bg extension.
cat $ids | parallel --compress "bedtools genomecov -ibam $S2PATH/{}.bam -split -bg | sort -k1,1 -k2,2n > $S2PATH/{}.sort.bg"

source activate bed2wig
# Convert each bedGraph coverage into bigWig coverage. The files will have the .bw extension.
cat $ids | parallel "bedGraphToBigWig $S2PATH/{}.sort.bg  $IDX.fai $S2PATH/{}.bw"
conda deactivate

