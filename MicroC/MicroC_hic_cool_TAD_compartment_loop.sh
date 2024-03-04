###
# This a example for bash command lines for Micro-C data processing based on pair files generated by data preprocessing
# including filtering out low-quality read pairs, generating hic/cool files, calling TAD domans and compartment, loop identification
# the read pair file will be used in the downstream analysis for generating hic/cool file, 
# the TAD/compartment/loop analysis will be based on hic/cool file
###

echo 'running'

pair_new='shNT_plusCL_mapq5_merge.pairs'
## get high quality reads: remove mapping quality scores less than 5, and remove read pair distance less than 100bp due to potential random ligation
awk '{
    if ( $0 ~ /^#chromsize/ ) {
        if ( $0 ~ /chr[0-9]|chrX|chrY|chrM/ ) {
            print
        }
    }else if ( $0 ~ /^#samheader/ && $0 ~ /SQ/) {
        if ( $0 ~ /chr[0-9]|chrX|chrY|chrM/ ) {
            print
        }
    }else if ( $0 ~ /^#/ ){
        print
    }else{
        if ( $2 ~ /chr[0-9]|chrX|chrY|chrM/ && $4 ~ /chr[0-9]|chrX|chrY|chrM/ && $9 > 5 && $10 > 5 && (($3 - $5) > 100 || ($3 - $5) < -100)){
            print
        }
    }
}' ../${pair_new} > ${pair_new} 


source /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/activate

## convet to .hic and mcool
echo "+++convert pair file to hic file using juicer tool for different resolution"

jarfile='/lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/MicroC/20220926_preadip/chromap/juicer_tools_1.22.01.jar'
chromsize=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/mm10.chromSize
/lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/java -Xmx48000m  -Djava.awt.headless=true \
 -jar ${jarfile} pre -q 5 --threads 16 -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000 \
 ${pair_new} ${pair_new}.hic ${chromsize}

echo "+++convert pair file to cool file using cooler for different resolution"
for r in 2500000 1000000 500000 250000 100000 50000 25000 10000 5000 2000 1000
do
    echo ++++${r}
    cooler_path=cool/${pair_new}_${r}.cool
    cooler cload pairs --assembly mm10 -c1 2 -p1 3 -c2 4 -p2 5 ${chromsize}:${r} ${pair_new} ${cooler_path}
    echo ++++balance
    cooler balance ${cooler_path}
done

# call TAD domains using arrowhead algorith in juicer tools
/lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/java -Xmx10g -jar ${jarfile} arrowhead -k KR -r 1000000 ${pair_new}.hic ${pair_new}_domain --ignore-sparsity
/lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/java -Xmx10g -jar ${jarfile} arrowhead -k KR -r 100000 ${pair_new}.hic ${pair_new}_domain --ignore-sparsity
/lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/java -Xmx10g -jar ${jarfile} arrowhead -k KR -r 10000 ${pair_new}.hic ${pair_new}_domain --ignore-sparsity
/lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/java -Xmx10g -jar ${jarfile} arrowhead -k KR -r 5000 ${pair_new}.hic ${pair_new}_domain --ignore-sparsity

### compartment calling using cooltool, a python package

r=5000  ## we do multiple resolution, here is an example
cond=shNT_plusCL
python compartment.py cool/${cond}_mapq5_merge.pairs_${r}.cool /lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/GRCm38.primary_assembly.genome.fa ./cooltool_comp/${cond}_merge_mapq5_$r 5


### loop/dots identification using peakachu software at 50kb resolution, using plusCL as an example
source /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/activate

peakachu score_genome -r 5000 --balance -p ../../cool/shNT_plusCL_mapq5_merge.pairs_5000.cool -O shNT_plusCL_mapq5_merge.pairs_5000.bedpe -m ./models/high-confidence.250million.5kb.w6.pkl --minimum-prob 0

peakachu pool -r 5000 -i shNT_plusCL_mapq5_merge.pairs_5000.bedpe -o shNT_plusCL_mapq5_merge.pairs_5000_loops.bedpe -t 0.95

echo "Finish"
#done


