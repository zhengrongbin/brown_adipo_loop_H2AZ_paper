#sample=NTplusCL
#strand=forward

sample=$1
strand=$2
rep=$3
ID=$4

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


