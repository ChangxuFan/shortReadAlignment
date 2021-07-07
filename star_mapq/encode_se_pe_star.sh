#!/bin/bash -e
if [ $# -ne 7 ]; then
    echo "usage v1: lrna_align_star_pe_se.sh <star_index_directory> <read1.fq.gz><read2.fq.gz> <library_id> <ncpus><out_put_dir><root_name_for_bam_file>"
    echo "Align single-end or pair-end reads with STAR. adopted from encode https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/dnanexus/align-star-se/resources/usr/bin/lrna_align_star_se.sh."
    echo "also incorparated: https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/DAC/STAR_RSEM.sh"
    echo "if sample is single-end, read2 should be 'single'"
    echo "i you don't know what to use as library_id, just leave SRR number there"
    exit -1; 
fi
star_index=$1  # STAR Index archive.
read1_fq_gz=$(readlink --canonicalize $2)    # gzipped fastq of of single-end reads.
if [ $3 != "single" ]
then
    read2_fq_gz=$(readlink --canonicalize $3)
else
    read2_fq_gz=""
fi
library_id=$4      # Library identifier which will be added to bam header.
ncpus=$5            # Number of cpus available.
bam_root="$7_star" # root name for output bam (e.g. "out_bam" will create "out_bam_star_genome.bam" and "out_bam_star_anno.bam") 
outdir=$6
echo "-- Alignments file will be: '${bam_root}_genome.bam' and '${bam_root}_anno.bam'"

mkdir -p $outdir
cd $outdir
echo "-- Set up headers..."
set -x
libraryComment="@CO\tLIBID:${library_id}"
echo -e ${libraryComment} > COfile.txt
# cat $outdir/*_bamCommentLines.txt >> COfile.txt
echo `cat COfile.txt`
set +x

echo "-- Map reads..."
set -x
STAR --genomeDir $star_index --readFilesIn $read1_fq_gz $read2_fq_gz                                \
    --readFilesCommand zcat --runThreadN $ncpus --genomeLoad LoadAndKeep     \
    --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1    \
    --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04              \
    --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000         \
    --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate   \
    --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
    --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate                \
    --quantMode TranscriptomeSAM --sjdbScore 1 --limitBAMsortRAM 60000000000  \
#    --outSAMprimaryFlag AllBestScore # note I disabled this one intentionally so that there will be no reads 
# reported as primary more than once.
# I was going to allBestScore but ended up deciding against it because I will filter out reads with low mapq anyway.
mv Aligned.sortedByCoord.out.bam ${bam_root}_genome.bam
mv Log.final.out ${bam_root}_Log.final.out
set +x
ls -l ${bam_root}_genome.bam

echo "-- Sorting annotation bam..."
set -x
cat <( samtools view -H Aligned.toTranscriptome.out.bam ) \
    <( samtools view -@ $ncpus Aligned.toTranscriptome.out.bam | sort -S 60G -T ./ ) | \
    samtools view -@ $ncpus -bS - > ${bam_root}_anno.bam
set +x
ls -l ${bam_root}_anno.bam

echo "-- Collect bam flagstats..."
set -x
samtools flagstat ${bam_root}_genome.bam > ${bam_root}_genome_flagstat.txt
samtools flagstat ${bam_root}_anno.bam > ${bam_root}_anno_flagstat.txt
set +x

echo "-- The results..."
ls -l ${bam_root}*

