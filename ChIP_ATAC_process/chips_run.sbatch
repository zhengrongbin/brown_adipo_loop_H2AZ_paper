## an example command line for excute CHIPs pipeline for ChIP-seq or ATAC-seq data proprocessing

echo "++++ activating anaconda"
source /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/activate
conda activate chips

echo "+++++ run CHIPs"

cd /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/H2AZ_2022_1st

snakemake -s CHIPS/chips.snakefile --rerun-incomplete -j 8 > run.out
#snakemake -s CHIPS/chips.snakefile -j 8 -npr --configfile config.yaml 

echo "+++++ Finished"



