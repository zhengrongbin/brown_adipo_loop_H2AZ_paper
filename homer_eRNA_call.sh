#!/bin/bash
#SBATCH --partition=cbp-compute # queue to be used
#SBATCH --time=15:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=pro # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # send and email when the job begins, ends or fails
#SBATCH --mail-user=your_email_address # Email address to send the job status
#SBATCH --output=output_%A_%a.txt # Name of the output file
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=4 # Number of cpu cores on one node
#SBATCH --mem=20G

echo "++++ activating anaconda"

## NTplusCL, replicate a sample, forward
sample=NTplusCL
strand=forward 
rep=a
ID=1600

# ### NTplusCL, replicate a sample, reverse
# sample=NTplusCL
# strand=reverse 
# rep=a 
# ID=1600

# ### NTplusCL, replicate b sample, forward
# sample=NTplusCL
# strand=forward 
# rep=b
# ID=1603

# ### NTplusCL, replicate b sample, reverse
# sample=NTplusCL
# strand=reverse 
# rep=b 
# ID=1603

# ### KD_plusCL, replicate a sample, forward
# sample=KDplusCL
# strand=forward 
# rep=a
# ID=1601

# ### KD_plusCL, replicate a sample, reverse
# sample=KDplusCL
# strand=reverse 
# rep=a 
# ID=1601

# ### KD_plusCL, replicate b sample, forward
# sample=KDplusCL
# strand=forward 
# rep=b
# ID=1604

# ### KD_plusCL, replicate b sample, reverse
# sample=KDplusCL
# strand=reverse 
# rep=b 
# ID=1604

# ### NTminusCL, replicate a sample, forward
# sample=NTminusCL
# strand=forward 
# rep=a
# ID=1599

# ### NTminusCL, replicate a sample, reverse
# sample=NTminusCL
# strand=reverse 
# rep=a 
# ID=1599

# ### NTminusCL, replicate b sample, forward
# sample=NTminusCL
# strand=forward 
# rep=b
# ID=1602

# ### NTminusCL, replicate b sample, reverse
# sample=NTminusCL
# strand=reverse 
# rep=b 
# ID=1602

## ===== run ====
source /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/activate

conda activate chips

bedgraph=/temp_work/ch228298/Coll_Yang/Nature_revise_Yang/PROseq_102024/bedgraphs/Tseng001_PR${ID}_${sample}_${rep}_mm10_dedup_5pr_${strand}.bedGraph

prefix=${sample}_${strand}_${rep}
## bdg to bed
awk '{print $1"\t"$2"\t"$3"\t.\t"$4"\t+"}' ${bedgraph} |tail -n +3 > ${prefix}_bdg_to.bed

makeTagDirectory ${prefix} ${prefix}_bdg_to.bed -format bed -force5th

findPeaks ${prefix} -style factor -center -fragLength 150 -size 200 -minDist 400 -tbp 0 -L 3 -localSize 10000 -fdr 0.001 -o ${prefix}_homer_peaks.bed

## intergenic
bedtools subtract -a ${prefix}_bdg_to.bed -b /lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/gencode.vM23.annotation.gene_annotation.bed > ${prefix}_bdg_to.bed_intergenic.bed

prefix1=${prefix}_intergenic

makeTagDirectory ${prefix1} ${prefix}_bdg_to.bed_intergenic.bed -format bed -force5th

findPeaks ${prefix1} -style factor -center -fragLength 150 -size 200 -minDist 400 -tbp 0 -L 3 -localSize 10000 -fdr 0.001 -o ${prefix}_homer_peaks_intergenic_call.bed


echo "+++++ Finished"



