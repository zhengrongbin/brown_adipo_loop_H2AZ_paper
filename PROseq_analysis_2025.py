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

