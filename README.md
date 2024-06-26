# scripts and notebooks for Micro-C, ChIP-seq, ATAC-seq, and RNA-seq data analysis
### dependency and version control
<li>Python==3.8, R==4.1</li>
<li>Micro-C analysis: Galore (version 0.6.6), BWA (version 0.7.17-r1188), pairtools (version 1.0.0), COOLER (version 0.9.1), JUICER (version 1.22.01), peakachu (downloaded at July 2023), BedTools (version 2.31.0), pybedtools (version 0.9.1)</li>
<li>RNA-seq analysis: Galore (version 0.6.6), STAR (version 2.7.9a), deepTools (version 3.3.1), SamTools (version 1.16.1), DESeq2 (R package, version 1.32.0), fgsea (R package, version 1.18.0), gseapy (Python package, version 1.0.6)</li>
<li>ChIP-seq and ATAC-seq analysis: CHIPs pipeline (https://github.com/liulab-dfci/CHIPS), DANPOS2, BedTools (version 2.31.0), pybedtools (Python package, version 0.9.1), pyBigWig (Python package, version 0.3.20)</li>

### Installation
The script can be run on a normal computer with enough RAM. However, we recommend using the high-performance computing server with Linux operating system as it can provide enough computig resource. These software can be installed by either one of the four approaches:
<li>1. anaconda, for example, Galore can be installed from anaconda at https://anaconda.org/bioconda/trim-galore</li>
<li>Python packages can be installed from "pip install" command lines</li>
<li>R packages can be installed Bioconductor using the BiocManager::install() function</li>
<li>[BWA](https://github.com/lh3/bwa), [JUICER](https://github.com/aidenlab/juicer), [CHIPs](https://f1000research.com/articles/10-517), [DANPOS2](https://sites.google.com/site/danposdoc/) can be installed from their home page.</li>

### Downloading the data
<p>The raw data for RNA-seq, Micro-C, ChIP-seq, and ATAC-seq data has been deposited to NCBI GEO repository under accession numbers:</p>
<li>GSE261413 for RNA-seq</li>
<li>GSE261416 for Micro-C</li>
<li>GSE261412 for ChIP-seq</li>
<li>GSE261410 for ATAC-seq</li>

### RNA-seq
<li>RNAseq_process/RNAseq_data_preprocessing_cmd.sh</li>
<p>This a example for bash command lines for RNA-seq data processing, including read mapping, read counting, generating bigwig, calculating FPKM and TPM for normalized expression index. The count matrix will be used for differential expression analysis in DESeq2 (a R package). This task is expected to finish within 1hr per sample</p>

<li>RNAseq_process/diffexp_analysis.R</li>
<p>The R script to perform differential expression analysis using DESeq2, following by gene ontology analysis and gene set enrichment ananlysis (GSEA). The input for this script is the count matrix and comparison design. This task is expected to finish within 20mins per analysis.</p>

<li>RNAseq_process/diffexp_run.sh</li>
<p>a example bash command line to excute diffexp_analysis.R</p>

<li>RNAseq_process/comparison_design.csv</li>
<p>a example of comparison design in our study, the 1-labeled samples versus 0-labeled samples in the comparison</p>


### Micro-C
<li>MicroC/MicroC_data_preprocessing_cmd.sh</li>
<p>This a example for bash command lines for Micro-C data processing, including read mapping, read counting, generating BAM file and read pair file. The read pair file will be used in the downstream analysis for generating hic/cool file, the TAD/compartment/loop analysis will be based on hic/cool file. This task is expected to finish within 2 days per sample.</p>

<li>MicroC/MicroC_hic_cool_TAD_compartment_loop.sh</li>
<p>This a example for bash command lines for Micro-C data processing based on pair files generated by data preprocessing, including filtering out low-quality read pairs, generating hic/cool files, calling TAD domans and compartment, loop identification. The read pair file will be used in the downstream analysis for generating hic/cool file, the TAD/compartment/loop analysis will be based on hic/cool file. This task is expected to finish within a half day per sample.</p>

<li>MicroC/compartment.py</li>
<p>A python script to perform compartment analysis based on cool files using cooltools.</p>

### ChIP-seq and ATAC-seq
<li>ChIP_ATAC_process/chips_run.sbatch</li>
<p>an example command line for excute CHIPs pipeline for ChIP-seq or ATAC-seq data proprocessing. [CHIPs](https://f1000research.com/articles/10-517) is snakemake pipeline for chromatin data processing. This will produce needed outputs for downstream analysis, including mapped reads, ChIP or ATAC signal, quality control results, etc. This task is expected to finish within 1 day per sample.</p>

<li>ChIP_ATAC_process/chips_config_for_ATAC.yaml</li>
<p>An example config file for ATAC-seq data processing, required for CHIPs pipeline</p>

<li>ChIP_ATAC_process/chips_config_for_ChIPseq.yaml</li>
<p>An example config file for ChIP-seq data processing, required for CHIPs pipeline</p>

<li>ChIP_ATAC_process/chips_metasheet_for_ATAC.csv</li>
<p>An example metasheet table for ATAC-seq data processing, required for CHIPs pipeline</p>

<li>ChIP_ATAC_process/chips_metasheet_for_ChIPseq.csv</li>
<p>An example metasheet table for ChIP-seq data processing, required for CHIPs pipeline</p>

<li>ChIP_ATAC_process/danpos_run.sh</li>
<p>An example for running [DANPOS2](https://sites.google.com/site/danposdoc/) software to call peaks for ChIP-seq and ATAC-seq data</p>

### Downstream analysis
<li>Compartment_TAD_analysis_03012024.ipynb or Compartment_TAD_analysis_03012024.py</li>
<p>The jupytor notebook or python script for Compartment and TAD downstream analysis based the outputs from proprocessing steps</p>

<li>MicroC_H2AZ_20240229.ipynb or MicroC_H2AZ_20240229.py</li>
<p>The jupytor notebook or python script for loop analysis and integration with H2AZ ChIP-seq, ATAC-seq, RNA-seq, and GWAS based the outputs from proprocessing steps</p>

<li>exp_analysis_02292024.ipynb or exp_analysis_02292024.py</li>
<p>The jupytor notebook or python script for RNA-seq data downstream analysis based the outputs from proprocessing steps</p>
