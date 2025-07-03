#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os,sys
import numpy as np
import pandas as pd
import collections
import pickle as pk
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import colors as mcolors
import seaborn as sns
from adjustText import adjust_text
from scipy.stats import pearsonr, spearmanr, chisquare, ttest_ind, ranksums, wilcoxon, fisher_exact, mannwhitneyu
from scipy import stats
import pybedtools
import cooler
import cooltools.lib.plotting
import cooltools
from scipy.stats import gaussian_kde
import statsmodels.api as sm
from coolpuppy import coolpup
from coolpuppy.lib import numutils
from coolpuppy.lib.puputils import divide_pups
from coolpuppy import plotpup
import cooler
import bioframe
from cooltools import expected_cis, expected_trans
from skmisc.loess import loess
import pybedtools
import pyBigWig
from matplotlib_venn import venn2, venn3

plt.rcParams.update(plt.rcParamsDefault)
rc={"axes.labelsize": 16, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "figure.titleweight":"bold", #"font.size":14,
    "figure.figsize":(5.5,4.2), "font.weight":"regular", "legend.fontsize":10,
    'axes.labelpad':8, 'figure.dpi':300}
plt.rcParams.update(**rc)



# ### PRO-seq vs RNA-seq

# In[3]:


## load differential expression result and expression TPM matrix
out = open('../../MicroC/20221212_RNA-seq_diffexp_repCorrect.pk', 'rb')
diffexp = pk.load(out)
out.close()


# In[5]:


pro_diffexp = {}
d = pd.read_csv('../../DataProcess/PROseq_102024/HMD_Core/CoreAnalysis/differential_expression/NTplus.vs.NTminus/Tseng001_deseq_NTplus.vs.NTminus_results_table.txt',
                                                   sep = '\t').dropna()
norm_pro_mat = d[['symbol', 'NormCounts_NTminusA','NormCounts_NTplusA','NormCounts_KDplusA','NormCounts_NTminusB', 'NormCounts_NTplusB','NormCounts_KDplusB']]
norm_pro_mat = norm_pro_mat.groupby('symbol').max()

d = d[['symbol', 'log2FoldChange_NTplus.vs.NTminus', 'stat_NTplus.vs.NTminus', 'padj_NTplus.vs.NTminus']]
d.columns = ['symbol', 'log2FoldChange', 'stat', 'padj']
d = d[~d['symbol'].duplicated()]
d.index = d['symbol'].tolist()
pro_diffexp['shNT_plusCL_vs_minusCL'] = d

d = pd.read_csv('../../DataProcess/PROseq_102024/HMD_Core/CoreAnalysis/differential_expression/KDplus.vs.NTplus/Tseng001_deseq_KDplus.vs.NTplus_results_table.txt',
                                                   sep = '\t').dropna()
d = d[['symbol', 'log2FoldChange_KDplus.vs.NTplus', 'stat_KDplus.vs.NTplus', 'padj_KDplus.vs.NTplus']]
d.columns = ['symbol', 'log2FoldChange', 'stat', 'padj']
d = d[~d['symbol'].duplicated()]
d.index = d['symbol'].tolist()
pro_diffexp['plusCL_KD_vs_shNT'] = d


# In[50]:


## plusCL vs minusCL
plot_df = pd.concat([diffexp['shNT_plusCL_vs_minusCL'][['log2FoldChange']],
           pro_diffexp['shNT_plusCL_vs_minusCL'][['log2FoldChange']]], axis = 1).dropna()
plot_df.columns = ['RNA-seq', 'PRO-seq']
plot_df = plot_df[['PRO-seq', 'RNA-seq']]

de_genes = pro_diffexp['shNT_plusCL_vs_minusCL'].query('abs(log2FoldChange) > 0 and padj < 0.05').index.tolist()

fig, ax = plt.subplots(figsize = (4,4))
sns.heatmap(data = plot_df[plot_df.index.isin(de_genes)].sort_values('PRO-seq'), cmap = 'bwr',
            center = 0, vmax = 1.5, vmin = -1.5, yticklabels = False)
plt.tight_layout()
plt.show()
plt.close()

## KD vs shNT
plot_df2 = pd.concat([diffexp['plusCL_KD_vs_shNT'][['log2FoldChange']],
           pro_diffexp['plusCL_KD_vs_shNT'][['log2FoldChange']]], axis = 1).dropna()
plot_df2.columns = ['RNA-seq', 'PRO-seq']
plot_df2 = plot_df2[['PRO-seq', 'RNA-seq']]

de_genes = pro_diffexp['plusCL_KD_vs_shNT'].query('abs(log2FoldChange) > 0 and padj < 0.05').index.tolist()

fig, ax = plt.subplots(figsize = (4,4))
sns.heatmap(data = plot_df2[plot_df2.index.isin(de_genes)].sort_values('PRO-seq'),
            cmap = 'bwr', center = 0, vmax = 1, vmin = -1, yticklabels=False)
plt.tight_layout()
plt.show()
plt.close()


# ### eRNA analysis

# In[ ]:


## gene annotations for TSS, TTS, promoter, and gene body for mouse mm10
mm10_size = '/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/mm10.chromSize'

## TSS
mm10_tss_ann = pd.read_csv('/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/gencode.vM23.annotation.protein_coding.tss.csv',
                          sep = '\t', header = None)
mm10_tss_ann.columns = ['chr', 'start', 'end', 'name', 'label', 'strand']

## TTS
mm10_protein_ann = pd.read_csv('/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/gencode.vM23.annotation.protein_coding.csv')
mm10_protein_ann_pos = mm10_protein_ann.query('strand == "+"')
mm10_tts_ann_pos = mm10_protein_ann_pos[['seqname', 'start', 'end', 'gene_name', 'frame', 'strand']]
mm10_tts_ann_pos['start'] = mm10_tts_ann_pos['end'].copy()
mm10_tts_ann_pos['end'] = mm10_tts_ann_pos['start']+1
mm10_tts_ann_pos.columns = ['chr', 'start', 'end', 'name', 'label', 'strand']

mm10_protein_ann_neg = mm10_protein_ann.query('strand == "-"')
mm10_tts_ann_neg = mm10_protein_ann_neg[['seqname', 'start', 'end', 'gene_name', 'frame', 'strand']]
mm10_tts_ann_neg['end'] = mm10_tts_ann_neg['start']+1
mm10_tts_ann_neg.columns = ['chr', 'start', 'end', 'name', 'label', 'strand']
mm10_tts_ann = pd.concat([mm10_tts_ann_pos, mm10_tts_ann_neg])
mm10_tts_ann['start'] = [x if x > 0 else 0 for x in mm10_tts_ann['start'].tolist()]

d = 2000
## promoter define from mm10 genome
mm10_tss_2k = mm10_tss_ann.copy()
mm10_tss_2k['start'] = mm10_tss_2k['start'] - d
mm10_tss_2k['end'] = mm10_tss_2k['end'] + d
mm10_tss_2k = mm10_tss_2k[mm10_tss_2k['start']>0]
mm10_tss_2k_bed = pybedtools.BedTool.from_dataframe(mm10_tss_2k)

## transcription termination site
## promoter define from mm10 genome
mm10_tts_2k = mm10_tts_ann.copy()
mm10_tts_2k['start'] = mm10_tts_2k['start'] - d
mm10_tts_2k['end'] = mm10_tts_2k['end'] + d
mm10_tts_2k = mm10_tts_2k[mm10_tts_2k['start']>0]
mm10_tts_2k_bed = pybedtools.BedTool.from_dataframe(mm10_tts_2k)


d = 5000
## promoter define from mm10 genome
mm10_tss_5k = mm10_tss_ann.copy()
mm10_tss_5k['start'] = mm10_tss_5k['start'] - d
mm10_tss_5k['end'] = mm10_tss_5k['end'] + d
mm10_tss_5k = mm10_tss_5k[mm10_tss_5k['start']>0]
mm10_tss_2k_bed = pybedtools.BedTool.from_dataframe(mm10_tss_5k)

## transcription termination site
## promoter define from mm10 genome
mm10_tts_5k = mm10_tts_ann.copy()
mm10_tts_5k['start'] = mm10_tts_5k['start'] - d
mm10_tts_5k['end'] = mm10_tts_5k['end'] + d
mm10_tts_5k = mm10_tts_5k[mm10_tts_5k['start']>0]
mm10_tts_5k_bed = pybedtools.BedTool.from_dataframe(mm10_tts_5k)


mm10_gene = pd.read_csv('/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/gencode.vM23.annotation.gene_annotation.csv')
mm10_gene_forward = mm10_gene.query('strand == "+"')[['seqname', 'start', 'end', 'gene_name', 'strand', 'gene_type']]
mm10_gene_reverse = mm10_gene.query('strand == "-"')[['seqname', 'start', 'end', 'gene_name', 'strand', 'gene_type']]


# In[ ]:


## PRO-seq for defining actively transcribled genes 
d = pd.read_csv('../DataProcess/PROseq_102024/HMD_Core/CoreAnalysis/differential_expression/NTplus.vs.NTminus/Tseng001_deseq_NTplus.vs.NTminus_results_table.txt',
                                                   sep = '\t').dropna()
norm_pro_mat = d[['symbol', 'NormCounts_NTminusA','NormCounts_NTplusA','NormCounts_KDplusA','NormCounts_NTminusB', 'NormCounts_NTplusB','NormCounts_KDplusB']]
norm_pro_mat = norm_pro_mat.groupby('symbol').max()

## set the pro-seq cutoff to 100 (norm_cout) for transcribing genes
transcrible_gene = {x.replace('NormCounts_', ''): norm_pro_mat[x][norm_pro_mat[x]>100].index.tolist() for x in norm_pro_mat.columns.tolist()}


# In[ ]:


### run homer to call eRNA peaks for each sample, forward and reverse strand separated
### bash: bash homer_eRNA_call.sh


# In[ ]:


## filters eRNA among homer called peaks
def get_erna_filter(cond, perb, strand, rep, intergenic = False):
    if intergenic:
        erna_all_path = '../DataProcess/PROseq_102024/eRNA/homer_eRNA_peaks/%s%sCL_%s_%s_homer_peaks_intergenic_call.bed'%(perb, cond, strand, rep.lower())
    else:
        erna_all_path = '../DataProcess/PROseq_102024/eRNA/homer_eRNA_peaks/%s%sCL_%s_%s_homer_peaks.bed'%(perb, cond, strand, rep.lower())
    erna_all = [line.rstrip().split('\t')[1:]+[line.rstrip().split('\t')[0]] for line in open(erna_all_path).readlines() if not line.startswith('#')]
    erna_all = pd.DataFrame(erna_all)
    erna_all_bed = pybedtools.BedTool.from_dataframe(erna_all)
    
    ## remove peaks in tss +- 2k regions
    erna_all_bed_filter = erna_all_bed.intersect(mm10_tss_2k_bed, v = True, wa = True, f=1.0)
    ## remove peaks in tts +- 2k regions
    erna_all_bed_filter = erna_all_bed_filter.intersect(mm10_tts_2k_bed, v = True, wa = True, f=1.0)
    ## remove peaks with embeding into the transcription of the same direction
    sample_label = perb+cond+rep
    if strand == "forward":
        gene_bed = pybedtools.BedTool.from_dataframe(mm10_gene_forward[mm10_gene_forward['gene_name'].isin(transcrible_gene[sample_label])])
    else:
        gene_bed = pybedtools.BedTool.from_dataframe(mm10_gene_reverse[mm10_gene_reverse['gene_name'].isin(transcrible_gene[sample_label])])
    
    erna_all_bed_filter = erna_all_bed_filter.intersect(gene_bed, v = True, wa = True, f=1.0)
    
    erna_all_bed_filter = erna_all_bed_filter.to_dataframe().drop_duplicates()
    erna_all_bed_filter.columns = ['chr','start','end','strand','Normalized Tag Count','focus ratio','findPeaks Score', 'Fold Change vs Local','p-value vs Local','Clonal Fold Change', 'PeakID']
    # remove < 0.2 RPM peaks, 'Normalized Tag Count' was per 10M, so set cutoff to 0.02
    erna_all_bed_filter = erna_all_bed_filter[erna_all_bed_filter['Normalized Tag Count'].astype('float') > 0.02]
    
    if intergenic:
        erna_all_bed_filter['PeakID'] = erna_all_bed_filter['PeakID']+'-inter'
    return(erna_all_bed_filter)

## load eRNA peaks in all regions call, apply basic filters
erna_res = {}
erna_res['NTplusA_forward'] = get_erna_filter(cond='plus', perb='NT', strand='forward', rep='A')
erna_res['NTplusB_forward'] = get_erna_filter(cond='plus', perb='NT', strand='forward', rep='B')
erna_res['NTplusA_reverse'] = get_erna_filter(cond='plus', perb='NT', strand='reverse', rep='A')
erna_res['NTplusB_reverse'] = get_erna_filter(cond='plus', perb='NT', strand='reverse', rep='B')

erna_res['NTminusA_forward'] = get_erna_filter(cond='minus', perb='NT', strand='forward', rep='A')
erna_res['NTminusB_forward'] = get_erna_filter(cond='minus', perb='NT', strand='forward', rep='B')
erna_res['NTminusA_reverse'] = get_erna_filter(cond='minus', perb='NT', strand='reverse', rep='A')
erna_res['NTminusB_reverse'] = get_erna_filter(cond='minus', perb='NT', strand='reverse', rep='B')

erna_res['KDplusA_forward'] = get_erna_filter(cond='plus', perb='KD', strand='forward', rep='A')
erna_res['KDplusB_forward'] = get_erna_filter(cond='plus', perb='KD', strand='forward', rep='B')
erna_res['KDplusA_reverse'] = get_erna_filter(cond='plus', perb='KD', strand='reverse', rep='A')
erna_res['KDplusB_reverse'] = get_erna_filter(cond='plus', perb='KD', strand='reverse', rep='B')

## load eRNA peaks in intergenic region specific peak call, apply basic filters
## this is to rescue peaks proximal to promoter
erna_inter_res = {}
erna_inter_res['NTplusA_forward'] = get_erna_filter(cond='plus', perb='NT', strand='forward', rep='A', intergenic=True)
erna_inter_res['NTplusB_forward'] = get_erna_filter(cond='plus', perb='NT', strand='forward', rep='B', intergenic=True)
erna_inter_res['NTplusA_reverse'] = get_erna_filter(cond='plus', perb='NT', strand='reverse', rep='A', intergenic=True)
erna_inter_res['NTplusB_reverse'] = get_erna_filter(cond='plus', perb='NT', strand='reverse', rep='B', intergenic=True)

erna_inter_res['NTminusA_forward'] = get_erna_filter(cond='minus', perb='NT', strand='forward', rep='A', intergenic=True)
erna_inter_res['NTminusB_forward'] = get_erna_filter(cond='minus', perb='NT', strand='forward', rep='B', intergenic=True)
erna_inter_res['NTminusA_reverse'] = get_erna_filter(cond='minus', perb='NT', strand='reverse', rep='A', intergenic=True)
erna_inter_res['NTminusB_reverse'] = get_erna_filter(cond='minus', perb='NT', strand='reverse', rep='B', intergenic=True)

erna_inter_res['KDplusA_forward'] = get_erna_filter(cond='plus', perb='KD', strand='forward', rep='A', intergenic=True)
erna_inter_res['KDplusB_forward'] = get_erna_filter(cond='plus', perb='KD', strand='forward', rep='B', intergenic=True)
erna_inter_res['KDplusA_reverse'] = get_erna_filter(cond='plus', perb='KD', strand='reverse', rep='A', intergenic=True)
erna_inter_res['KDplusB_reverse'] = get_erna_filter(cond='plus', perb='KD', strand='reverse', rep='B', intergenic=True)


# In[ ]:


## merge and deduplicate for peaks in all and intergenic
erna_res_all = {}
for x in erna_res:
    ## remove the same peak in two peak calling 
    tmp1 = erna_res[x].copy() ## all
    tmp1['cord'] = tmp1['chr']+':'+tmp1['start'].astype('str')+'-'+tmp1['end'].astype('str')
    tmp2 = erna_inter_res[x].copy() ## intergeneic call
    tmp2['cord'] = tmp2['chr']+':'+tmp2['start'].astype('str')+'-'+tmp2['end'].astype('str')
    ## intersect
    interpeak = np.intersect1d(tmp1['cord'], tmp2['cord'])
    ## remove from one call
    tmp2 = tmp2[~tmp2['cord'].isin(interpeak)]
    
    erna_res_all[x] = pd.concat([tmp1, tmp2]).drop_duplicates()
    for i in ['start', 'end', 'Normalized Tag Count', 'focus ratio', 'findPeaks Score', 'Fold Change vs Local', 'p-value vs Local','Clonal Fold Change']:
        if i in ['start', 'end']:
            erna_res_all[x][i] = erna_res_all[x][i].astype('int')
        else:
            erna_res_all[x][i] = erna_res_all[x][i].astype('float')
    erna_res_all[x] = erna_res_all[x].sort_values(['chr','start','end'])
    erna_res_all[x]['center'] = erna_res_all[x]['start']+(erna_res_all[x]['end'] - erna_res_all[x]['start']).astype('int')
    erna_res_all[x]['uid'] = ['peak'+str(ii) for ii in range(erna_res_all[x].shape[0])]


# #### identify enhancer location using eRNA peaks

# In[ ]:


## define eRNA regions based on peaks
def _get_erna_region_(cond):
    enhancer_regions = []
    for ch in ['chr'+str(x) for x in range(1, 20)]+['chrX']:
        d = erna_res_all[cond+'_forward'].query('chr == @ch')
        d.index = d['uid'].tolist()
        
        d1 = erna_res_all[cond+'_reverse'].query('chr == @ch')
        d1.index = d1['uid'].tolist()
        
        a=np.array(d['center'])
        b=np.array(d1['center'])
        
        diff = a[np.newaxis,:]-b[:,np.newaxis]
        diff = pd.DataFrame(diff, columns = ('forward_'+d['uid']).tolist(),
                            index = ('reverse_'+d1['uid']).tolist())
        
        ## define enhancer centers and regions
        erna_collect = []
        paired_list = []
        unpaired_list = []
        for f in diff.columns.tolist():
            col = diff[f]
            r = (col[(abs(col) <= 1000)]).index.tolist()
            if len(r) != 0:
                paired_list.extend([f]+r)
                start_erna = d1.loc[[x.replace('reverse_', '') for x in r], ].sort_values('center', ascending = False).iloc[0,:] ## take the max
                end_erna = d.loc[f.replace('forward_', ''), ] # only one should be matched, so no max
                loci = pd.Series([start_erna['center'], end_erna['center']], index = ['reverse_'+start_erna['uid'], 'forward_'+end_erna['uid']]).sort_values()
                start, end = loci.tolist()
                erna_collect.append([ch, start, end, int(start+int(.5*(end-start))), 'paired', '; '.join(loci.index.tolist())]) ## add one erna
            else:
                unpaired_list.append(f) # only for forward, will check reverse later
        ## check unpaired for reverse
        re_all = ('reverse_'+d1['uid'])
        unpaired_list.extend(re_all[~re_all.isin(paired_list)])
        
        ## add erna for unpair, 180bp upstream for center
        for p in unpaired_list:
            if p.startswith('forward_'):
                uid = p.replace('forward_', '')
                c = d.loc[uid, 'center']
                erna_collect.append([ch, c-360, c, int(c-180), 'unpaired', 'forward_'+uid])
            else:
                uid = p.replace('reverse_', '')
                c = d1.loc[uid, 'center']
                erna_collect.append([ch, c, c+360, int(c+180), 'unpaired', 'forward_'+uid])
        enhancer_regions.extend(erna_collect)
        
    enhancer_regions = pd.DataFrame(enhancer_regions).drop_duplicates()
    enhancer_regions['name'] = ['eRNA_'+str(x) for x in range(enhancer_regions.shape[0])]
    enhancer_regions = enhancer_regions[[0,1,2,'name', 3, 4, 5]]
    enhancer_regions.columns = ['chr', 'start', 'end', 'name', 'center', 'pair', 'pair_erna']
    return(enhancer_regions)
    
## enhancer locations (eRNA regions)
erna_finial_df = {}
for cond in ['NTplusA', 'NTplusB', 'NTminusA', 'NTminusB', 'KDplusA', 'KDplusB']:
    erna_finial_df[cond] = _get_erna_region_(cond = cond)

## define reproducibable eRNA regions between replicates
erna_finial_df_reprod = {'NTplus': pybedtools.BedTool.from_dataframe(erna_finial_df['NTplusA']).intersect(pybedtools.BedTool.from_dataframe(erna_finial_df['NTplusB']), wa = True).to_dataframe().drop_duplicates(),
                        'NTminus': pybedtools.BedTool.from_dataframe(erna_finial_df['NTminusA']).intersect(pybedtools.BedTool.from_dataframe(erna_finial_df['NTminusB']), wa = True).to_dataframe().drop_duplicates(),
                        'KDplus': pybedtools.BedTool.from_dataframe(erna_finial_df['KDplusB']).intersect(pybedtools.BedTool.from_dataframe(erna_finial_df['KDplusA']), wa = True).to_dataframe().drop_duplicates()}
for x in erna_finial_df_reprod:
    erna_finial_df_reprod[x].columns = ['chr', 'start', 'end', 'name', 'center', 'pair', 'erna_name']

## union all eRNA from three conditions for the downstream comparison
union_erna = pd.DataFrame()
for x in erna_finial_df_reprod:
    tmp = erna_finial_df_reprod[x][['chr', 'start', 'end', 'name']]
    tmp['name'] = x+'_'+tmp['name']
    union_erna = pd.concat([union_erna, tmp])

union_erna.drop_duplicates(['chr', 'start', 'end']).to_csv('union_reprod_erna.bed', sep = '\t', index = None, header = None)


# #### compute enhancer activity using PROseq signal on eRNA regions

# In[ ]:


## calcualte eRNA levels from 5' bedGraph files, bedgraph converted to bigwig file using bedGraphToBigWig with default setting
## erna transcription level, 1kb around the center
def _erna_level_(d, transcrible_gene_list, pro_bw_op):
    forward_gene_bed = pybedtools.BedTool.from_dataframe(mm10_gene_forward[mm10_gene_forward['gene_name'].isin(transcrible_gene_list)])
    reverse_gene_bed = pybedtools.BedTool.from_dataframe(mm10_gene_reverse[mm10_gene_reverse['gene_name'].isin(transcrible_gene_list)])

    d['center'] = d['start']+(.5*(d['end']-d['start'])).astype('int')
    d['start'] = (d['center']-500).astype('int')
    d['end'] = (d['center']+500).astype('int')
    
    d_bed = pybedtools.BedTool.from_dataframe(d)
    
    ## eRNA in transcrible gene regions
    forward_d = d_bed.intersect(forward_gene_bed, wa=True, f = 1.0).to_dataframe().drop_duplicates()
    forward_d.columns = d.columns
    reverse_d = d_bed.intersect(reverse_gene_bed, wa=True, f = 1.0).to_dataframe().drop_duplicates()
    reverse_d.columns = d.columns
    inter = np.intersect1d(forward_d['name'], reverse_d['name'])
    forward_d = forward_d[~forward_d['name'].isin(inter)]
    reverse_d = reverse_d[~reverse_d['name'].isin(inter)]
    
    ## eRNA in intergenic regions
    other_d = d[~d['name'].isin(forward_d['name'].tolist()+reverse_d['name'].tolist())]
    
    bw_values = []
    # print(forward_d)
    ## eRNA in forward transcrible genes, take reverse reads for eRNA level
    for i, line in forward_d.iterrows():
        c, s, e = line['chr'], line['start'], line['end']
        v = []
        for x in ['NTplus_a', 'NTplus_b', 'NTminus_a', 'NTminus_b', 'KDplus_a', 'KDplus_b']:
            v.append(pro_bw_op[x+'_reverse'].stats(c, s, e, type = 'sum')[0])
        bw_values.append(v)

    ## eRNA in reverse transcrible genes, take forward reads for eRNA level
    for i, line in reverse_d.iterrows():
        c, s, e = line['chr'], line['start'], line['end']
        v = []
        for x in ['NTplus_a', 'NTplus_b', 'NTminus_a', 'NTminus_b', 'KDplus_a', 'KDplus_b']:
            v.append(pro_bw_op[x+'_forward'].stats(c, s, e, type = 'sum')[0])
        bw_values.append(v)

    ## eRNA in intergenic region, take average of forward and reverse reads for eRNA level
    for i, line in other_d.iterrows():
        c, s, e = line['chr'], line['start'], line['end']
        v = []
        for x in ['NTplus_a', 'NTplus_b', 'NTminus_a', 'NTminus_b', 'KDplus_a', 'KDplus_b']:
            v.append(pd.Series([pro_bw_op[x+'_forward'].stats(c, s, e, type = 'sum')[0], pro_bw_op[x+'_reverse'].stats(c, s, e, exact = True)[0]]).dropna().mean())   
        bw_values.append(v)
    bw_df = pd.DataFrame(bw_values, index = forward_d['name'].tolist()+reverse_d['name'].tolist()+other_d['name'].tolist(),
                    columns = ['NTplus_a', 'NTplus_b', 'NTminus_a', 'NTminus_b', 'KDplus_a', 'KDplus_b'])
    bw_df.index = bw_df.index
    return(bw_df)

## extract bw using 5' signal
pro_bw = {'KDplus_a_forward':'../DataProcess/PROseq_102024/HMD_Core/bedgraphs/Tseng001_PR1601_KDplusCL_a_mm10_dedup_5pr_forward.bedGraph.bw',
          'KDplus_a_reverse':'../DataProcess/PROseq_102024/HMD_Core/bedgraphs/Tseng001_PR1601_KDplusCL_a_mm10_dedup_5pr_reverse.bedGraph.bw',

          'KDplus_b_forward':'../DataProcess/PROseq_102024/HMD_Core/bedgraphs/Tseng001_PR1604_KDplusCL_b_mm10_dedup_5pr_reverse.bedGraph.bw',
          'KDplus_b_reverse':'../DataProcess/PROseq_102024/HMD_Core/bedgraphs/Tseng001_PR1604_KDplusCL_b_mm10_dedup_5pr_forward.bedGraph.bw',

          'NTminus_a_forward':'../DataProcess/PROseq_102024/HMD_Core/bedgraphs/Tseng001_PR1599_NTminusCL_a_mm10_dedup_5pr_forward.bedGraph.bw',
          'NTminus_a_reverse':'../DataProcess/PROseq_102024/HMD_Core/bedgraphs/Tseng001_PR1599_NTminusCL_a_mm10_dedup_5pr_reverse.bedGraph.bw',
          
          'NTminus_b_forward':'../DataProcess/PROseq_102024/HMD_Core/bedgraphs/Tseng001_PR1602_NTminusCL_b_mm10_dedup_5pr_forward.bedGraph.bw',
          'NTminus_b_reverse':'../DataProcess/PROseq_102024/HMD_Core/bedgraphs/Tseng001_PR1602_NTminusCL_b_mm10_dedup_5pr_reverse.bedGraph.bw',
          
          'NTplus_a_forward':'../DataProcess/PROseq_102024/HMD_Core/bedgraphs/Tseng001_PR1600_NTplusCL_a_mm10_dedup_5pr_forward.bedGraph.bw',
          'NTplus_a_reverse':'../DataProcess/PROseq_102024/HMD_Core/bedgraphs/Tseng001_PR1600_NTplusCL_a_mm10_dedup_5pr_reverse.bedGraph.bw',
          
          'NTplus_b_forward':'../DataProcess/PROseq_102024/HMD_Core/bedgraphs/Tseng001_PR1603_NTplusCL_b_mm10_dedup_5pr_forward.bedGraph.bw',
          'NTplus_b_reverse':'../DataProcess/PROseq_102024/HMD_Core/bedgraphs/Tseng001_PR1603_NTplusCL_b_mm10_dedup_5pr_reverse.bedGraph.bw',}
## load bw files
pro_bw_op = {x:pyBigWig.open(pro_bw[x]) for x in pro_bw}

## union transcribled genes, which will be used in determing eRNA activity
union_trans_genes = []
for x in transcrible_gene:
    union_trans_genes.extend(transcrible_gene[x])
## compute eRNA activity 
union_erna_level = _erna_level_(union_erna, transcrible_gene_list = union_trans_genes, pro_bw_op=pro_bw_op)
union_erna_level = union_erna_level.fillna(0)
union_erna_level.to_csv('union_erna_level_sum.csv')



# In[ ]:


## normalization to sequence depth and peak length
## total ref reads
total_reads = pd.Series([48991590, 61279698, 58701583, 51734800, 55674279, 55738922],
                       index = ['KDplus_a', 'KDplus_b', 'NTminus_a', 'NTminus_b', 'NTplus_a', 'NTplus_b'])
## eRNA length
erna_len = pd.Series(np.abs(union_erna['end']-union_erna['start']).tolist(), index = union_erna['name'].tolist())

## norm: (x * 1000 / erna_len) * 1e06 / total_reads
union_erna_level_norm = pd.DataFrame(index = union_erna_level.index, columns = union_erna_level.columns)
for x in union_erna_level_norm.columns.tolist():
    union_erna_level_norm[x] = (union_erna_level[x] * 1000 / erna_len.reindex(union_erna_level.index)) * 1e06 / total_reads[x]
union_erna_level_norm = np.log2(union_erna_level_norm+1)

## averaged by replicates
union_erna_level_norm_df = pd.DataFrame(
    [union_erna_level_norm[['NTplus_a', 'NTplus_b']].T.mean(),
    union_erna_level_norm[['NTminus_a', 'NTminus_b']].T.mean(),
    union_erna_level_norm[['KDplus_a', 'KDplus_b']].T.mean()],
    index = ['NTplus', 'NTminus', 'KDplus']
).T
union_erna_level_norm_df.to_csv('./erna_deseq_normed_log2_rep_average_table.csv')

