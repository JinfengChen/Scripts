#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --time=4:00:00
#SBATCH --output=From_bam_to_local_assembly.sh.stdout
#SBATCH -p intel
#SBATCH --chdir=./


CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=$1
fi

module load ncbi-blast/2.2.26
module load samtools
module load muscle
module load bwa
module load picard

velvet=/rhome/cjinfeng/BigData/software/Velvet/velvet
seqtk=/rhome/cjinfeng/BigData/software/seqtk-master/seqtk


bam=GN158.mPingSV.bam
region=Chr1:29318601-29323757
prefix=GN158_Chr1_29M
outdir=$prefix

if [ ! -e $prefix\.fa ]; then
    echo "Prepare reference sequence"
    echo $region | sed 's/:/\t/' | sed 's/-/\t/' > $prefix\.bed
    $seqtk subseq MSU_r7.fa $prefix\.bed > $prefix\.fa
    formatdb -i $prefix\.fa -p F 
fi

if [ ! -e $outdir/$prefix\_1.fq ]; then
    echo "Extract fastq from bam: $bam"
    mkdir $outdir
    samtools view -b $bam $region -o $prefix\.bam 
    #samtools fastq $bam -0 $outdir/$prefix\_unpaired.fq -1 $outdir/$prefix\_1.fq -2 $outdir/$prefix\_2.fq
    samtools sort -n $prefix\.bam -o $prefix\.sorted.bam
    bedtools bamtofastq -i $prefix\.sorted.bam -fq $outdir/$prefix\_1.fq -fq2 $outdir/$prefix\_2.fq
fi

if [ ! -e $prefix\.assembly ]; then
    echo "Local assembly"
    $velvet/velveth $prefix\.assembly 31 -shortPaired -fastq -separate $outdir/$prefix\_1.fq $outdir/$prefix\_2.fq
    $velvet/velvetg $prefix\.assembly -ins_length 200 -exp_cov 50 -min_contig_lgth 200 -scaffolding yes
    $seqtk seq -r $prefix\.assembly/contigs.fa > $prefix\.assembly/contigs.rec.fa
    cat $prefix\.assembly/contigs.fa $prefix\.assembly/contigs.rec.fa $prefix\.fa > $prefix\.compare.fa
fi

if [ ! -e $prefix\.assembly/contigs.rec.bam ]; then
    echo "Map read back to assembly"
    # reverse complement
    bwa index $prefix\.assembly/contigs.rec.fa
    picard CreateSequenceDictionary R=$prefix\.assembly/contigs.rec.fa O=$prefix\.assembly/contigs.rec.dict
    bwa mem $prefix\.assembly/contigs.rec.fa $outdir/$prefix\_1.fq $outdir/$prefix\_2.fq | samtools view -b - -o $prefix\.assembly/contigs.rec.raw.bam
    samtools sort $prefix\.assembly/contigs.rec.raw.bam -o $prefix\.assembly/contigs.rec.bam
    samtools index $prefix\.assembly/contigs.rec.bam
    # raw
    bwa index $prefix\.assembly/contigs.fa
    picard CreateSequenceDictionary R=$prefix\.assembly/contigs.fa O=$prefix\.assembly/contigs.dict
    bwa index $prefix\.assembly/contigs.fa
    bwa mem $prefix\.assembly/contigs.fa $outdir/$prefix\_1.fq $outdir/$prefix\_2.fq | samtools view -b - -o $prefix\.assembly/contigs.raw.bam
    samtools sort $prefix\.assembly/contigs.raw.bam -o $prefix\.assembly/contigs.bam
    samtools index $prefix\.assembly/contigs.bam
fi

if [ ! -e $prefix\.blast.output ]; then
    blastall -p blastn -i $prefix\.assembly/contigs.fa -d $prefix\.fa -e 1e-10 -o $prefix\.blast.output
    muscle -in $prefix\.compare.fa -out $prefix\.compare.msa
    muscle -clw -in $prefix\.compare.fa -out $prefix\.compare.clw
fi

echo "Done"
