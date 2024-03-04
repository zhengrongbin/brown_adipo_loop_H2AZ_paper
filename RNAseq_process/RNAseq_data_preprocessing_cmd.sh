###
# This a example for bash command lines for RNA-seq data processing
# including read mapping, read counting, generating bigwig, calculating FPKM and TPM for normalized expression index
# the count matrix will be used for differential expression analysis in DESeq2 (a R package)
###

echo '++++ activating biogrids'
source /programs/biogrids.shrc

## ==== set input ======
jobName=shNT_plusCL_rep1
commpath=/project/RC_Cardio-Chen-e2/
fq1=Fastq/RNA2_R1.fq.gz
fq2=Fastq/RNA2_R2.fq.gz
species=mm10

## ==== set output dir =====
workpath=/project/RC_Cardio-Chen-e2//rongbinzheng/BrownAdipo/analysis/shNT_plusCL_rep1
mkdir -p /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1

mkdir -p /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/QC

mkdir -p /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis

echo '++++++ STAR version'
STAR --version

### ========== index path ==========

mainpath=/project/RC_Cardio-Chen-e2//rongbinzheng/Genome/mm10
fapath=/project/RC_Cardio-Chen-e2//rongbinzheng/Genome/mm10/GRCm39.primary_assembly.genome.fa
gtfpath=/project/RC_Cardio-Chen-e2//rongbinzheng/Genome/mm10/gencode.vM23.annotation.gtf
indexpath=/project/RC_Cardio-Chen-e2//rongbinzheng/Genome/mm10/STAR_mm10
rsem_ref=/project/RC_Cardio-Chen-e2//rongbinzheng/Genome/mm10/rsem_ref/rsem_ref


## =========== process ============
echo '++++ automatic trimming'
trim_galore --version >> /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/version.control

trim_galore --paired --retain_unpaired --dont_gzip -o /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/QC --fastqc_args '-d /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/QC' Fastq/RNA2_R1.fq.gz Fastq/RNA2_R2.fq.gz



echo '++++ star runnning'
echo '\n\nSTAR version' >> /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/version.control
STAR --version >> /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/version.control

STAR --runThreadN 8 --genomeDir /project/RC_Cardio-Chen-e2//rongbinzheng/Genome/mm10/STAR_mm10 --readFilesIn /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/QC/RNA2_R1_val_1.fq /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/QC/RNA2_R2_val_2.fq --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --quantMode TranscriptomeSAM GeneCounts --sjdbGTFfile /project/RC_Cardio-Chen-e2//rongbinzheng/Genome/mm10/gencode.vM23.annotation.gtf --outFileNamePrefix /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/STAR_output

## ======= TPM, FPKM using RSEM
echo '+++++++ RSEM'
rsem-calculate-expression --version  >> /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/version.control

rsem-calculate-expression --bam --no-bam-output -p 8     --paired-end --forward-prob 0.5     /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/STAR_outputAligned.toTranscriptome.out.bam     /project/RC_Cardio-Chen-e2//rongbinzheng/Genome/mm10/rsem_ref/rsem_ref /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/rsem

## ======= make bigwig file
echo '++++++++ generating bigwig'
bamCoverage --version >> /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/version.control
samtools --version|head -1 >> /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/version.control

samtools index /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/STAR_outputAligned.sortedByCoord.out.bam

bamCoverage --bam /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/STAR_outputAligned.sortedByCoord.out.bam --outFileName /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/shNT_plusCL_rep1.bw --normalizeUsing RPKM --outFileFormat bigwig

## ======= change name and move ======
echo '++++++++ move, remove, change result file name'

mkdir -p /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/result
mv /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/shNT_plusCL_rep1.bw /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/result/shNT_plusCL_rep1.bw
mv /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/rsem.genes.results /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/result/shNT_plusCL_rep1.gene.rsem.txt
mv /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/STAR_outputReadsPerGene.out.tab /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/result/shNT_plusCL_rep1.gene_count.txt
mv /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/rsem.isoforms.results /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/result/shNT_plusCL_rep1.isoforms.rsem.txt

rm -rf /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/STAR_output_STARgenome
rm -rf /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/STAR_output_STARtmp
rm -rf /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/QC/val*fq
rm -rf /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/analysis/*.toTranscriptome.out.bam 

## ======= annotate genes =====
python rnaseq_output_merge.py ann_gene -t /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/result/shNT_plusCL_rep1.gene.rsem.txt -s mm10
python rnaseq_output_merge.py ann_gene -c /project/RC_Cardio-Chen-e2//rongbinzheng/DataProcess/Yang_shH2AZ_RNA/shNT_plusCL_rep1/result/shNT_plusCL_rep1.gene_count.txt -s mm10

echo '++++++ finished'

