HOMEDIR=/data/fleming/fousteri-te

# STAR index generation
STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeDir $HOME/reference
  --genomeFastaFiles $HOME/reference/genome.fa \
  --sjdbGTFfile $HOME/reference/genes.gtf

# STAR alignment allowing multimapping for TE analysis
for FILE in `ls $HOMEDIR/fastq/*.fastq.gz`
do
  SAMPLE=`basename $FILE | sed s/\.fastq\.gz//`
  STAR \
    --runThreadN 16 \
    --genomeDir $HOMEDIR/reference \
    --runMode alignReads \
    --readFilesIn $FILE \
    --outSAMtype BAM SortedByCoordinate \
    --winAnchorMultimapNmax 100 \
    --outFilterMultimapNmax 100 \
    --readFilesCommand zcat \
    --outFileNamePrefix $HOMEDIR/bam/$SAMPLE
done

# Run TE transcripts with BAM files (and replicates) for each comparison
$HOMEDIR/TEtranscripts/bin/TEtranscripts \
  --treatment $HOMEDIR/bam/CSB1NOUVAligned.sortedByCoord.out.bam $HOMEDIR/bam/CSB3NOUVAligned.sortedByCoord.out.bam \
  --control $HOMEDIR/bam/VH103NOUVAligned.sortedByCoord.out.bam $HOMEDIR/bam/VH10DmNOUVAligned.sortedByCoord.out.bam \
  --GTF $HOMEDIR/reference/genes.gtf \
  --TE $HOMEDIR/teanalysis/hg19_rmsk_TE.gtf \
  --format BAM \
  --mode multi \
  --sortByPos

$HOMEDIR/TEtranscripts/bin/TEtranscripts \
  --treatment $HOMEDIR/bam/CSBprimary1Aligned.sortedByCoord.out.bam $HOMEDIR/bam/CSBprimary2Aligned.sortedByCoord.out.bam $HOMEDIR/bam/CSBprimary3Aligned.sortedByCoord.out.bam \
  --control $HOMEDIR/bam/VH103NOUVAligned.sortedByCoord.out.bam $HOMEDIR/bam/VH10DmNOUVAligned.sortedByCoord.out.bam \
  --GTF $HOMEDIR/reference/genes.gtf \
  --TE $HOMEDIR/teanalysis/hg19_rmsk_TE.gtf \
  --format BAM \
  --mode multi \
  --sortByPos

$HOMEDIR/TEtranscripts/bin/TEtranscripts \
  --treatment $HOMEDIR/bam/CSBprimary1Aligned.sortedByCoord.out.bam $HOMEDIR/bam/CSBprimary2Aligned.sortedByCoord.out.bam $HOMEDIR/bam/CSBprimary3Aligned.sortedByCoord.out.bam \
  --control $HOMEDIR/bam/CSB1NOUVAligned.sortedByCoord.out.bam $HOMEDIR/bam/CSB3NOUVAligned.sortedByCoord.out.bam \
  --GTF $HOMEDIR/reference/genes.gtf \
  --TE $HOMEDIR/teanalysis/hg19_rmsk_TE.gtf \
  --format BAM \
  --mode multi \
  --sortByPos

