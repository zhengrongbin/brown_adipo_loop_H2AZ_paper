#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os,sys
import numpy as np
import pandas as pd
import collections
# plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import colors as mcolors
import seaborn as sns
from adjustText import adjust_text

plt.rcParams.update(plt.rcParamsDefault)
rc={"axes.labelsize": 16, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "figure.titleweight":"bold", #"font.size":14,
    "figure.figsize":(5.5,4.2), "font.weight":"regular", "legend.fontsize":10,
    'axes.labelpad':8, 'figure.dpi':300}
plt.rcParams.update(**rc)





# ## expression GSEA plot
# 
# 
# 

# In[2]:


## read fgsea result for plusCL vs minusCL
shNT_plusCL_vs_minusCL_fgsea_res = pd.read_csv('./diffexp_repCorrect/_shNT_plusCL_vs_minusCL_fgseaRes.tsv', sep = '\t')



# In[3]:


## dot plot to show top pathways
plot_df = shNT_plusCL_vs_minusCL_fgsea_res.query('padj < 0.0025').sort_values('NES')
plot_df['-log10padj'] = -np.log10(plot_df['padj'])

colormap = plt.cm.Reds_r #or any other colormap
normalize = mpl.colors.Normalize(vmin=0, vmax=0.01)
fig, ax = plt.subplots(figsize = (8, 6))
for i in range(plot_df.shape[0]):
    ax.hlines(plot_df.iloc[i]['pathway'], 0, plot_df.iloc[i]['NES'], linewidth = 2, color = 'grey')
sp = ax.scatter(plot_df['NES'], range(plot_df.shape[0]), s = plot_df['size'],
               c = plot_df['padj'], cmap = colormap, norm=normalize, zorder = 100)
ax.set_ylim(-0.5, plot_df.shape[0])
ax.vlines(0, *ax.get_ylim(), color = 'black')
plt.legend(*sp.legend_elements("sizes", num=3), title = 'Gene Count',
          loc = 'center left', bbox_to_anchor=(1, 0.9), frameon = False)
ax.set(xlabel = 'Normalized Enrichment Score')
ax.tick_params(axis = 'y', labelsize = 10)
plt.colorbar(sp, shrink = .4, label = 'padj')
sns.despine(trim = False, left = True)
plt.tight_layout()
fig.savefig('plusCL_vs_minusCL_fgsea_top.pdf')
plt.show()
plt.close()


# In[4]:


## read fgsea result for KD vs shNT plusCL 

plusCL_KD_vs_shNT_fgsea_res = pd.read_csv('./diffexp_repCorrect/_plusCL_KD_vs_shNT_fgseaRes.tsv', sep = '\t')


# In[5]:


comp_pathway = ['mmu00010 Glycolysis %2F Gluconeogenesis',
               'mmu01212 Fatty acid metabolism',
               'mmu04714 Thermogenesis',
               'mmu01200 Carbon metabolism',
               'mmu03320 PPAR signaling pathway',
               'mmu01230 Biosynthesis of amino acids',
               'mmu00500 Starch and sucrose metabolism',
               'mmu00190 Oxidative phosphorylation']



# In[195]:


pplot_df1 = shNT_plusCL_vs_minusCL_fgsea_res[shNT_plusCL_vs_minusCL_fgsea_res['pathway'].isin(comp_pathway)]

colormap = plt.cm.Reds_r #or any other colormap
normalize = mpl.colors.Normalize(vmin=0, vmax=0.01)
fig, ax = plt.subplots(figsize = (7.5, 4))
for i in range(pplot_df1.shape[0]):
    ax.hlines(pplot_df1.iloc[i]['pathway'], 0, pplot_df1.iloc[i]['NES'], linewidth = 2, color = 'grey')
sp = ax.scatter(pplot_df1['NES'], range(pplot_df1.shape[0]), s = pplot_df1['size'],
               c = pplot_df1['padj'], cmap = colormap, norm=normalize, zorder = 100)
ax.set_ylim(-0.5, pplot_df1.shape[0])
ax.vlines(0, *ax.get_ylim(), color = 'black')
plt.legend(*sp.legend_elements("sizes", num=3), title = 'Gene Count',
          loc = 'center left', bbox_to_anchor=(1, 0.9), frameon = False)
ax.set(xlabel = 'Normalized Enrichment Score')
plt.colorbar(sp, shrink = .4, label = 'padj')
sns.despine(trim = False, left = True)
plt.tight_layout()
fig.savefig('plusCL_vs_minusCL_fgsea_top_subset.pdf')
plt.show()
plt.close()


# In[196]:


## dot plot to show top pathways
pplot_df2 = pd.DataFrame()
for i in pplot_df1['pathway'].tolist():
    pplot_df2 = pd.concat([pplot_df2, plusCL_KD_vs_shNT_fgsea_res[plusCL_KD_vs_shNT_fgsea_res['pathway']==i]])

colormap = plt.cm.Reds_r #or any other colormap
normalize = mpl.colors.Normalize(vmin=0, vmax=0.15)
fig, ax = plt.subplots(figsize = (7.5, 4))
for i in range(pplot_df2.shape[0]):
    ax.hlines(pplot_df2.iloc[i]['pathway'], 0, pplot_df2.iloc[i]['NES'], linewidth = 2, color = 'grey')
sp = ax.scatter(pplot_df2['NES'], range(pplot_df2.shape[0]), s = pplot_df2['size'],
               c = pplot_df2['padj'], cmap = colormap, norm=normalize, zorder = 100)
ax.set_ylim(-0.5, pplot_df2.shape[0])
ax.vlines(0, *ax.get_ylim(), color = 'black')
plt.legend(*sp.legend_elements("sizes", num=3), title = 'Gene Count',
          loc = 'center left', bbox_to_anchor=(1, 0.9), frameon = False)
ax.set(xlabel = 'Normalized Enrichment Score')
plt.colorbar(sp, shrink = .4, label = 'padj')
sns.despine(trim = False, left = True)
plt.tight_layout()
fig.savefig('plusCL_KD_vs_WT_fgsea_top_subset.pdf')
plt.show()
plt.close()


# In[ ]:
## load differential expression result and expression TPM matrix
out = open('20221212_RNA-seq_diffexp_repCorrect.pk', 'rb')
diffexp = pk.load(out)
out.close()

kegg_pathway = pd.read_csv('../../CommonData/mouse_KEGG_terms_symbol.txt', sep = '\t', header = None, index_col = 0)
kegg_pathway_dict = {}
for x in kegg_pathway.index.tolist():
    kegg_pathway_dict[x] = list(set(kegg_pathway.loc[x].iloc[0].split(';')) & set(diffexp['plusCL_KD_vs_shNT'].index.tolist()))

cut = np.log2(1.5)
cut2 = -np.log2(1.5)
kegg_pathway_dict['CL_UP'] = diffexp['shNT_plusCL_vs_minusCL'].query('log2FoldChange > @cut and padj < 0.05').index.tolist()
kegg_pathway_dict['CL_DOWN'] = diffexp['shNT_plusCL_vs_minusCL'].query('log2FoldChange < @cut2 and padj < 0.05').index.tolist()

rnk = diffexp['plusCL_KD_vs_shNT'][['log2FoldChange']]
pre_res = gp.prerank(rnk=rnk, # or rnk = rnk,
                     gene_sets=kegg_pathway_dict,
                     threads=4,
                     min_size=5,
                     max_size=5000,
                     permutation_num=1000, # reduce number to speed up testing
                     outdir=None, # don't write to disk
                     seed=6,
                     verbose=True, # see what's going on behind the scenes
                    )
axs = pre_res.plot(terms='CL_UP')
axs = pre_res.plot(terms='mmu04714 Thermogenesis')

## [HUMAN] GSEA for human RNA-seq differential expression analysis
## load data diffexp_batchCorrect

diffexp = {'NT_plusFSK_vs_minusFSK': pd.read_csv('./_plusFSK_vs_minusFSK.protein_gene.csv', index_col = 0),
          'plusFSK_KD_vs_NT': pd.read_csv('./_KD_vs_NT.protein_gene.csv', index_col = 0)}
## load terms
gene_terms = {}
with open('./gene sets for thermogenesis.gmt') as gmt:
    for line in gmt:
        line = line.strip().split('\t')
        gene_terms[line[0]] = line[2:]
with open('./gene sets for lipolysis.gmt') as gmt:
    for line in gmt:
        line = line.strip().split('\t')
        gene_terms[line[0]] = line[2:]
## rank genes by LFC
rnk = diffexp['plusFSK_KD_vs_NT']['log2FoldChange']
pre_res = gp.prerank(rnk=rnk, # or rnk = rnk,
                     gene_sets=gene_terms, #['GO_Biological_Process_2025', 'KEGG_2021_Human'],
                     threads=4,
                     min_size=5,
                     max_size=500,
                     permutation_num=1000, # reduce number to speed up testing
                     outdir=None, # don't write to disk
                     seed=6,
                     verbose=True, # see what's going on behind the scenes
                    )

## plot
axs = pre_res.plot(terms='Regulation of lipolysis in adipocytes', ofname = 'regulation_of_lipolisys__H2AZKD_gsea.pdf')
axs = pre_res.plot(terms='Thermogenesis WP4321', ofname = 'Thermogenesis__H2AZKD_gsea.pdf')

