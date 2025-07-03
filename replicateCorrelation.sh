#!/bin/bash
#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=100:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=hicrep # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # send and email when the job begins, ends or fails
#SBATCH --mail-user=your_email_address # Email address to send the job status
#SBATCH --output=output_%A_%a.txt # Name of the output file
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=2 # Number of cpu cores on one node
#SBATCH --mem=50G

echo "++++ activating"
source /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/activate
conda activate H2AZ

## python hicrep at https://github.com/dejunlin/hicrep
r=5000

cool1=../../../DataProcess/Human_MicroC/merge_pairs_mapq5/cool/NT_plusFSK_rep1.pairs_${r}.cool
cool2=../../../DataProcess/Human_MicroC/merge_pairs_mapq5/cool/NT_plusFSK_rep2.pairs_${r}.cool
cool3=../../../DataProcess/Human_MicroC/merge_pairs_mapq5/cool/NT_plusFSK_rep3.pairs_${r}.cool
cool4=../../../DataProcess/Human_MicroC/merge_pairs_mapq5/cool/NT_plusFSK_rep4.pairs_${r}.cool

hicrep $cool1 $cool2 ./replicate_correlation/shNT_plusFSK_rep1_rep2_resolu_${r}.txt --h 3 --dBPMax 500000
hicrep $cool1 $cool3 ./replicate_correlation/shNT_plusFSK_rep1_rep3_resolu_${r}.txt --h 3 --dBPMax 500000
hicrep $cool1 $cool4 ./replicate_correlation/shNT_plusFSK_rep1_rep4_resolu_${r}.txt --h 3 --dBPMax 500000
hicrep $cool2 $cool3 ./replicate_correlation/shNT_plusFSK_rep2_rep3_resolu_${r}.txt --h 3 --dBPMax 500000
hicrep $cool2 $cool4 ./replicate_correlation/shNT_plusFSK_rep2_rep4_resolu_${r}.txt --h 3 --dBPMax 500000
hicrep $cool3 $cool4 ./replicate_correlation/shNT_plusFSK_rep3_rep4_resolu_${r}.txt --h 3 --dBPMax 500000


cool1=../../../DataProcess/Human_MicroC/merge_pairs_mapq5/cool/NT_minusFSK_rep1.pairs_${r}.cool
cool2=../../../DataProcess/Human_MicroC/merge_pairs_mapq5/cool/NT_minusFSK_rep2.pairs_${r}.cool
cool3=../../../DataProcess/Human_MicroC/merge_pairs_mapq5/cool/NT_minusFSK_rep3.pairs_${r}.cool
cool4=../../../DataProcess/Human_MicroC/merge_pairs_mapq5/cool/NT_minusFSK_rep4.pairs_${r}.cool

hicrep $cool1 $cool2 ./replicate_correlation/shNT_minusFSK_rep1_rep2_resolu_${r}.txt --h 3 --dBPMax 500000
hicrep $cool1 $cool3 ./replicate_correlation/shNT_minusFSK_rep1_rep3_resolu_${r}.txt --h 3 --dBPMax 500000
hicrep $cool1 $cool4 ./replicate_correlation/shNT_minusFSK_rep1_rep4_resolu_${r}.txt --h 3 --dBPMax 500000
hicrep $cool2 $cool3 ./replicate_correlation/shNT_minusFSK_rep2_rep3_resolu_${r}.txt --h 3 --dBPMax 500000
hicrep $cool2 $cool4 ./replicate_correlation/shNT_minusFSK_rep2_rep4_resolu_${r}.txt --h 3 --dBPMax 500000
hicrep $cool3 $cool4 ./replicate_correlation/shNT_minusFSK_rep3_rep4_resolu_${r}.txt --h 3 --dBPMax 500000


cool1=../../../DataProcess/Human_MicroC/merge_pairs_mapq5/cool/KD_plusFSK_rep1.pairs_${r}.cool
cool2=../../../DataProcess/Human_MicroC/merge_pairs_mapq5/cool/KD_plusFSK_rep2.pairs_${r}.cool
cool3=../../../DataProcess/Human_MicroC/merge_pairs_mapq5/cool/KD_plusFSK_rep3.pairs_${r}.cool
cool4=../../../DataProcess/Human_MicroC/merge_pairs_mapq5/cool/KD_plusFSK_rep4.pairs_${r}.cool

hicrep $cool1 $cool2 ./replicate_correlation/KD_plusFSK_rep1_rep2_resolu_${r}.txt --h 3 --dBPMax 500000
hicrep $cool1 $cool3 ./replicate_correlation/KD_plusFSK_rep1_rep3_resolu_${r}.txt --h 3 --dBPMax 500000
hicrep $cool1 $cool4 ./replicate_correlation/KD_plusFSK_rep1_rep4_resolu_${r}.txt --h 3 --dBPMax 500000
hicrep $cool2 $cool3 ./replicate_correlation/KD_plusFSK_rep2_rep3_resolu_${r}.txt --h 3 --dBPMax 500000
hicrep $cool2 $cool4 ./replicate_correlation/KD_plusFSK_rep2_rep4_resolu_${r}.txt --h 3 --dBPMax 500000
hicrep $cool3 $cool4 ./replicate_correlation/KD_plusFSK_rep3_rep4_resolu_${r}.txt --h 3 --dBPMax 500000


echo "+++++ Finished"