###
# This a example for bash command lines for Micro-C data processing
# including read mapping, read counting, generating BAM file and read pair file
# the read pair file will be used in the downstream analysis for generating hic/cool file, 
# the TAD/compartment/loop analysis will be based on hic/cool file
###


source /programs/biogrids.shrc

echo 'running'

fqname='MicroC-16_S2'
name='KD_plusCL_rep4'

mkdir -p $name
cat ${fqname}_L001_R1_001.fastq.gz ${fqname}_L002_R1_001.fastq.gz ${fqname}_L003_R1_001.fastq.gz ${fqname}_L004_R1_001.fastq.gz > ${name}/${name}.R1.fq.gz
cat ${fqname}_L001_R2_001.fastq.gz ${fqname}_L002_R2_001.fastq.gz ${fqname}_L003_R2_001.fastq.gz ${fqname}_L004_R2_001.fastq.gz > ${name}/${name}.R2.fq.gz

## trim adaptors
trim_galore --paired -o ${name} ${name}/${name}.R1.fq.gz ${name}/${name}.R2.fq.gz

## the path for adaptor trimmed fastq files
prefix=${name}/${name}
fq1=/temp_work/ch228298/MicroC/${name}/${name}.R1_val_1.fq.gz
fq2=/temp_work/ch228298/MicroC/${name}/${name}.R2_val_2.fq.gz


nThreads=16
## read mapping using bwa mem software
ref='/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/bwa_GRCm38/GRCm38.primary_assembly.genome.fa'
genome_file='/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/GRCm38.primary_assembly.genome.fa.Genome'

bwa mem -t $nThreads -SP5M $ref $fq1 $fq2 | samtools view -Shb - > ${prefix}.bam

source /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/activate
## generate BAM and pair file using pairtools (a python-based package)
mkdir -p ${name}/temp1
pairtools parse --add-columns mapq --walks-policy 5unique --max-inter-align-gap 30 --nproc-in $nThreads --nproc-out $nThreads --chroms-path $genome_file ${prefix}.bam | pairtools sort --tmpdir=${name}/temp1/ --nproc 8 |pairtools dedup --nproc-in $nThreads --nproc-out $nThreads --mark-dups --output-stats ${prefix}_stats.txt|pairtools split --nproc-in $nThreads --nproc-out $nThreads --output-pairs ${prefix}.pairs

