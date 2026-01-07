#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
plt.rcParams.update(plt.rcParamsDefault)
rc={"axes.labelsize": 16, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "figure.titleweight":"bold", #"font.size":14,
    "figure.figsize":(5.5,4.2), "font.weight":"regular", "legend.fontsize":10,
    'axes.labelpad':8, 'figure.dpi':300}
plt.rcParams.update(**rc)



# In[2]:


mm10_tss_ann = pd.read_csv('/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/gencode.vM23.annotation.protein_coding.tss.csv',
                          sep = '\t', header = None)

mm10_tss_ann.columns = ['chr', 'start', 'end', 'name', 'label', 'strand']


# In[3]:


fruitpunch = sns.blend_palette(['white', 'red'], as_cmap=True)
from matplotlib.ticker import EngFormatter
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import bioframe
import matplotlib.patches as patches

def rectangles_around_dots(dots_df, region, loc="upper", lw=1, ec="cyan", fc="none"):
    """
    yield a series of rectangles around called dots in a given region
    """
    # select dots from the region:
    df_reg = bioframe.select(
        bioframe.select(dots_df, region, cols=("chrom1","start1","end1")),
        region,
        cols=("chrom2","start2","end2"),
    )
    rectangle_kwargs = dict(lw=lw, ec=ec, fc=fc)
    # draw rectangular "boxes" around pixels called as dots in the "region":
    for s1, s2, e1, e2 in df_reg[["start1", "start2", "end1", "end2"]].itertuples(index=False):
        width1 = e1 - s1
        width2 = e2 - s2
        if loc == "upper":
            yield patches.Rectangle((s2, s1), width2, width1, **rectangle_kwargs)
        elif loc == "lower":
            yield patches.Rectangle((s1, s2), width1, width2, **rectangle_kwargs)
        else:
            raise ValueError("loc has to be uppper or lower")
bp_formatter = EngFormatter('b')
def format_ticks(ax, x=True, y=True, rotate=True):
    """format ticks with genomic coordinates as human readable"""
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)
        
def _cont_map_(clr, ch, cond, title='', 
             clr2=None, cond2=None, w = None,
             comp_df=None, comp_df2=None,
             decimal = 0, dots_df = False,
             cmap = 'fall', vmax = 0.01, vmin = 1e-5, figsize = (10, 5),
             dot_loc = 'upper', dot_ec = 'cyan', pdf = False):
    plt.close()
    fig, ax = plt.subplots(figsize=figsize)
    norm = LogNorm(vmax=vmax, vmin=vmin)
    if clr2 is None:
        mat = clr.matrix(balance=True).fetch(ch)
        ylabel = ''
        xlabel = ''
    else:
        mat1 = clr.matrix(balance=True).fetch(ch)
        mat2 = clr2.matrix(balance=True).fetch(ch)
        ## cond1 in lower tri, cond2 in upper tri
        w = 3 if w is None else w
        mat = np.tril(mat1, -w)+np.triu(mat2, w)
        ylabel = '%s'%cond
        xlabel = '%s'%cond2
    
    bins = clr.bins().fetch(ch)
    if ':' in ch:
        c = ch.split(':')[0]
        st, end = ch.split(':')[1].split('-')
        st, end = int(st), int(end)
    else:
        c, st, end = bins.iloc[0,0], bins.iloc[0,1], bins.iloc[-1,2]        

    im = ax.matshow(
        mat,
        norm=norm,
        cmap=cmap,
        extent=(st, end, end, st)
    )
    format_ticks(ax)
    ax.xaxis.tick_top()
    ax.tick_params(axis = 'x', rotation = 45)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im, label='Corrected frequencies', cax = cax, shrink = 0.3);
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.xaxis.set_label_position('top')  
    # format_ticks(ax)
#     plt.tight_layout()
    if type(dots_df) == type(pd.DataFrame()):
        for box in rectangles_around_dots(dots_df, (c, int(st), int(end)), lw=1.5,
                                         loc = dot_loc, ec = dot_ec):
            ax.add_patch(box)
    for x in ['top', 'right', 'bottom', 'left']:
        ax.spines[x].set_visible(False)

    ## add compartment
    if comp_df is not None:
        ax1 = divider.append_axes("bottom", size="20%", pad=0.1, sharex=ax)
        # weights = clr.bins()[:]['weight'].values
        comp1 = comp_df[(comp_df['chrom']==c) & (comp_df['start']>=int(st)) & (comp_df['end']<=int(end))].fillna(0)
        ax1.plot([st, end],[0,0],'k',lw=0.25)
#         ax1.plot(comp1['start'], comp1['E1'].fillna(0))
        ax1.fill_between(
                comp1['start'].tolist(), 0, comp1['E1'].tolist(), 
                where=(comp1['E1']>0),
                edgecolor = 'none',
                facecolor = plt.cm.get_cmap('tab10')(1),
                label = 'A',
                alpha=0.8, 
            )
        ax1.fill_between(
                comp1['start'].tolist(), 0, comp1['E1'].tolist(), 
                where=(comp1['E1']<0),
                facecolor = plt.cm.get_cmap('tab10')(0),
                label = 'B',
                alpha=0.8, 
            )
#         ax1.plot(comp1.query('E1 <= 0')['start'], comp1.query('E1 <= 0')['E1'],
#                 label='B', color = plt.cm.get_cmap('tab10')(0))
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), 
                   title='%s'%cond2, frameon = False)
#         ax1.set_xlabel(cond1)
        ax1.xaxis.set_visible(False)
        ax1.yaxis.set_visible(False)
        for x in ['top', 'right', 'bottom', 'left']:
            ax1.spines[x].set_visible(False)

    if comp_df2 is not None:
        ax2 = divider.append_axes("top", size="20%", pad=1, sharex=ax)
        # weights = clr.bins()[:]['weight'].values
        comp2 = comp_df2[(comp_df2['chrom']==c) & 
                         (comp_df2['start']>=int(st)) & 
                         (comp_df2['end']<=int(end))].fillna(0)
        ax2.plot([st, end],[0,0],'k',lw=0.25)
#         ax1.plot(comp1['start'], comp1['E1'].fillna(0))
        ax2.fill_between(
                comp2['start'].tolist(), 0, comp2['E1'].tolist(), 
                where=(comp2['E1']>0),
                edgecolor = 'none',
                facecolor = plt.cm.get_cmap('tab10')(1),
                label = 'A',
                alpha=0.8, 
            )
        ax2.fill_between(
                comp2['start'].tolist(), 0, comp2['E1'].tolist(), 
                where=(comp2['E1']<0),
                facecolor = plt.cm.get_cmap('tab10')(0),
                label = 'B',
                alpha=0.8, 
            )
#         ax1.plot(comp1.query('E1 <= 0')['start'], comp1.query('E1 <= 0')['E1'],
#                 label='B', color = plt.cm.get_cmap('tab10')(0))
        ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                   title='%s'%cond2, frameon = False)
        ax2.set_xlabel(cond1)
        ax2.xaxis.set_visible(False)
        ax2.yaxis.set_visible(False)
        for x in ['top', 'right', 'bottom', 'left']:
            ax2.spines[x].set_visible(False)
    plt.suptitle(title, fontweight = 'bold') if title else None
    plt.tight_layout()
    pdf.savefig(fig) if pdf else plt.show()
    plt.close()
    
def quantile_norm(X):
    """
    Normalize the columns of X to each have the same distribution.
  
    Given an expression matrix (microarray data, read counts, etc) of M genes
    by N samples, quantile normalization ensures all samples have the same
    spread of data.
  
    The data across each row are averaged to obtain an average column. Each
    column quantile is replaced with the corresponding quantile of the average
    column.
  
    Parameters
    ----------
    X : 2D array of float, shape (M, N)
        The input data, with M rows (genes/features) and N columns (samples).
  
    Returns
    -------
    Xn : 2D array of float, shape (M, N)
        The normalized data.
    """
    # compute the quantiles
    quantiles = np.mean(np.sort(X, axis=0), axis=1)
    # compute the column-wise ranks. Each observation is replaced with its
    # rank in that column: the smallest observation is replaced by 1, the
    # second-smallest by 2, ..., and the largest by M, the number of rows.
    ranks = np.apply_along_axis(stats.rankdata, 0, X)

    # convert ranks to integer indices from 0 to M-1
    rank_indices = ranks.astype(int) - 1

    # index the quantiles for each rank with the ranks matrix
    Xn = quantiles[rank_indices]

    return Xn

### function to draw boxplot for gene expression change of up, down, stable loops
def _plot_loop_exp_(mm10_tss_ann, up_loops, down_loops, stable_loops, de, col='log2FoldChange', d = 10000, figsize = (5,5), pdf = False):
    mm10_promoter = mm10_tss_ann.copy()
    mm10_promoter['start'] = mm10_tss_ann['start'] - d
    mm10_promoter['end'] = mm10_tss_ann['end'] + d
    mm10_promoter = mm10_promoter[mm10_promoter['start']>0]

    up_loops = up_loops.sort_values(['chr', 'start', 'end'])
    py_up_loops = pybedtools.BedTool.from_dataframe(up_loops)
    
    down_loops = down_loops.sort_values(['chr', 'start', 'end'])
    py_down_loops = pybedtools.BedTool.from_dataframe(down_loops)
    
    stable_loops = stable_loops.sort_values(['chr', 'start', 'end'])
    py_stable_loops = pybedtools.BedTool.from_dataframe(stable_loops)

    up_prom_gene = pybedtools.BedTool.from_dataframe(mm10_promoter).intersect(py_up_loops, wa=True).to_dataframe()
    down_prom_gene = pybedtools.BedTool.from_dataframe(mm10_promoter).intersect(py_down_loops, wa=True).to_dataframe()
    stable_prom_gene = pybedtools.BedTool.from_dataframe(mm10_promoter).intersect(py_stable_loops, wa=True).to_dataframe()

    up_prom_gene = np.unique([x.split(':')[-1] for x in up_prom_gene['name'].tolist()])
    down_prom_gene = np.unique([x.split(':')[-1] for x in down_prom_gene['name'].tolist()])
    stable_prom_gene = np.unique([x.split(':')[-1] for x in stable_prom_gene['name'].tolist()])

    plot_df1 = de.loc[de.index.isin(up_prom_gene),]
    plot_df1.loc[:,'label'] = 'Up'
    plot_df2 = de.loc[de.index.isin(down_prom_gene),]
    plot_df2.loc[:,'label'] = 'Down'
    plot_df3 = de.loc[de.index.isin(stable_prom_gene),]
    plot_df3.loc[:,'label'] = 'Stable'
    plot_df = pd.concat([plot_df1, plot_df2, plot_df3])
    
    fig, ax = plt.subplots(figsize = figsize)
    # col = 'log2FoldChange'
    sns.boxplot(data = plot_df, x = 'label', y = col, showfliers = False,
               width = .5)
    s1, p1 = ranksums(plot_df.query('label == "Up"')[col], 
                      plot_df.query('label == "Stable"')[col])
    s2, p2 = ranksums(plot_df.query('label == "Down"')[col], 
                      plot_df.query('label == "Stable"')[col])
    s3, p3 = ranksums(plot_df.query('label == "Up"')[col], 
                      plot_df.query('label == "Down"')[col])
    ax.hlines(0, *ax.get_xlim(), color = 'black', linestyle = 'dashed')
    ax.set(xlabel = 'Looping changes', ylabel = 'mRNA log2FoldChange')
    ax.set_title('Up vs Stable: s=%.2f, p=%.2e\nDown vs Stable: s=%.2f, p=%.2e\nUp vs Down: s=%.2f, p=%.2e'%(s1, p1, s2, p2, s3, p3))
    plt.tight_layout()
    plt.show() if not pdf else pdf.savefig(fig)
    plt.close()
    return(plot_df)


### function to get values for loop change and expression change, lowess smoothing applied
def _get_values_(df, comp = 'shNT_plusCL_vs_minusCL', value = 'log2FoldChange'):
    dx, dy = [], []
    ddf = []
    for i, line in df.iterrows():
        g = line[['r1_gene', 'r2_gene']].dropna().values.tolist()
        for j in g:
            j = j.split(';')
            for jj in j:
                ytmp = diffexp[comp].loc[jj, value] if jj in diffexp[comp].index.tolist() else np.nan
                dy.append(ytmp)
                dx.append(line['score_delta'])
                ddf.append(line.tolist()+[ytmp, jj])

    ddf = pd.DataFrame(ddf, columns = df.columns.tolist()+[value, 'gene'])
#     return(ddf)
    ddf = ddf[~pd.isna(ddf[value]) & ~pd.isna(ddf['score_delta'])]
    ddf = ddf.sort_values('score_delta')
    x, y = ddf['score_delta'], ddf[value]
    l = loess(x,y)
    l.fit()
    pred = l.predict(x, stderror=True)
    conf = pred.confidence()
    lowess = pred.values
    ll = conf.lower
    ul = conf.upper
    ddf['lowess'] = lowess
    ddf['lower_lowess'] = ll
    ddf['upper_lowess'] = ul
    return(ddf)



## define EP, PP, P-other
def _define_EP_PP_(loops, comp = 'shNT_plusCL_vs_minusCL'):
    ## EP
    q1 = '(r1_promoter == True and r2_ELS == True and r2_K27ac == True) or (r2_promoter == True and r1_ELS == True and r1_K27ac == True)'
    ## PP
    q2 = '(r1_promoter == True and r2_promoter == True)'
    ## P-other
    q3 = '(r1_promoter == True and (r2_ELS != True or r2_K27ac != True) and r2_promoter != True) or (r2_promoter == True and (r1_ELS != True or r1_K27ac != True) and r1_promoter != True)'

    ## overlap with promoter
    loops['r1'] = loops['chrom1']+':'+loops['start1'].astype('str')+'-'+loops['end1'].astype('str')
    loops['r2'] = loops['chrom2']+':'+loops['start2'].astype('str')+'-'+loops['end2'].astype('str')

    r1 = pybedtools.BedTool.from_dataframe(loops[['chrom1','start1','end1', 'r1']])
    r2 = pybedtools.BedTool.from_dataframe(loops[['chrom2','start2','end2', 'r2']])
    r1_promoter = r1.intersect(mm10_promoter_bed, wao = True).to_dataframe().drop_duplicates()
    r2_promoter = r2.intersect(mm10_promoter_bed, wao = True).to_dataframe().drop_duplicates()

    ## overlap with ELS
    r1_ELS = r1.intersect(mm10_ELS_bed, wao = True).to_dataframe().drop_duplicates()
    r2_ELS = r2.intersect(mm10_ELS_bed, wao = True).to_dataframe().drop_duplicates()

    ## overlap with K27ac
    r1_K27ac = r1.intersect(k27ac_union_peaks_bed, wao = True).to_dataframe().drop_duplicates()
    r2_K27ac = r2.intersect(k27ac_union_peaks_bed, wao = True).to_dataframe().drop_duplicates()

    ## overlap with ELS_K27ac
    r1_ELS_K27ac = r1.intersect(ELS_k27ac_union_peaks_bed, wao = True).to_dataframe().drop_duplicates()
    r2_ELS_K27ac = r2.intersect(ELS_k27ac_union_peaks_bed, wao = True).to_dataframe().drop_duplicates()

    loops['r1_promoter'] = loops['r1'].isin(r1_promoter.query('score != "."')['name'])
    loops['r2_promoter'] = loops['r2'].isin(r2_promoter.query('score != "."')['name'])
    loops['r1_ELS'] = loops['r1'].isin(r1_ELS.query('score != "."')['name'])
    loops['r2_ELS'] = loops['r2'].isin(r2_ELS.query('score != "."')['name'])
    loops['r1_K27ac'] = loops['r1'].isin(r1_K27ac.query('score != "."')['name'])
    loops['r2_K27ac'] = loops['r2'].isin(r2_K27ac.query('score != "."')['name'])
    loops['r1_ELS_K27ac'] = loops['r1'].isin(r2_ELS_K27ac.query('score != "."')['name'])
    loops['r2_ELS_K27ac'] = loops['r2'].isin(r2_ELS_K27ac.query('score != "."')['name'])
    
    r1_genes = r1_promoter.query('score != "."').groupby('name').apply(lambda df: ';'.join(set([x.split(':')[-1] for x in df['thickEnd'].unique().tolist()]))).to_dict()
    r2_genes = r2_promoter.query('score != "."').groupby('name').apply(lambda df: ';'.join(set([x.split(':')[-1] for x in df['thickEnd'].unique().tolist()]))).to_dict()


    loops['r1_gene'] = [r1_genes.get(x, np.nan) for x in loops['r1'].tolist()]
    loops['r2_gene'] = [r2_genes.get(x, np.nan) for x in loops['r2'].tolist()]

    loops['r1_mRNA_fc'] = [diffexp[comp].reindex(x.split(';'))['log2FoldChange'].mean() if not pd.isna(x) else np.nan for x in loops['r1_gene'].tolist()]
    loops['r2_mRNA_fc'] = [diffexp[comp].reindex(x.split(';'))['log2FoldChange'].mean() if not pd.isna(x) else np.nan for x in loops['r2_gene'].tolist()]

    ## define EP, PP, and PO
    EP_df = loops.query(q1)
    PP_df = loops.query(q2)
    PO_df = loops.query(q3)
    other_df = loops[~loops['label'].isin(EP_df['label'].tolist()+PP_df['label'].tolist()+PO_df['label'].tolist())]
    
    return(EP_df, PP_df, PO_df, other_df)

### identify H2AZ occupied loops
def _plot_H2AZ_overlap_v2_(loops, peaks, title = '', focus = [], flip = False):
    anchor1 = pybedtools.BedTool.from_dataframe(loops[['chrom1','start1','end1']])
    anchor2 = pybedtools.BedTool.from_dataframe(loops[['chrom2','start2','end2']])
    union_loop_an = pd.DataFrame(loops[['chrom1','start1','end1']].values.tolist()+loops[['chrom2','start2','end2']].values.tolist()).drop_duplicates()
    union_loop_an.columns = ['chrom', 'start', 'end']
    union_anchor = pybedtools.BedTool.from_dataframe(union_loop_an)
    H2AZ_overlap = {}
    if not focus:
        peaks_keys = list(peaks.keys())
    else:
        peaks_keys = focus
    for x in peaks_keys:
        tmp1 = anchor1.intersect(peaks[x], wa=True).to_dataframe().drop_duplicates()
        tmp2 = anchor2.intersect(peaks[x], wa=True).to_dataframe().drop_duplicates()
        tmp1['label'] = tmp1['chrom']+':'+tmp1['start'].astype('str')+'-'+tmp1['end'].astype('str')
        tmp2['label'] = tmp2['chrom']+':'+tmp2['start'].astype('str')+'-'+tmp2['end'].astype('str')
        loops_tmp = loops.copy()
        loops_tmp['anchor1'] = loops_tmp['chrom1']+':'+loops_tmp['start1'].astype('str')+'-'+loops_tmp['end1'].astype('str')
        loops_tmp['anchor2'] = loops_tmp['chrom2']+':'+loops_tmp['start2'].astype('str')+'-'+loops_tmp['end2'].astype('str')
        loops_tmp['anchor1_H2AZ'] = loops_tmp['anchor1'].isin(tmp1['label'])
        loops_tmp['anchor2_H2AZ'] = loops_tmp['anchor2'].isin(tmp2['label'])
        tmp_or = loops_tmp[((loops_tmp['anchor1_H2AZ']==True) | (loops_tmp['anchor2_H2AZ']==True)) & ~((loops_tmp['anchor1_H2AZ']==True) & (loops_tmp['anchor2_H2AZ']==True))]
        tmp_and = loops_tmp[(loops_tmp['anchor1_H2AZ']==True) & (loops_tmp['anchor2_H2AZ']==True)]
        H2AZ_overlap[x]={'anchor1': tmp1.shape[0], 'anchor2': tmp2.shape[0],
                                           'both_anchor':tmp_and.shape[0],'either_anchor':tmp_or.shape[0],
                                           'total_loop':loops.shape[0],
                                           'total_H2AZ':peaks[x].to_dataframe().shape[0],
                        'loops':{'both':tmp_and, 'either':tmp_or}}

    df1 = pd.DataFrame(H2AZ_overlap)
    df1 = df1.loc[['anchor1', 'anchor2', 'both_anchor', 'either_anchor']]*100/df1.loc['total_loop']
    df1 = df1.reset_index().melt(id_vars = ['index'])
    
    if flip:
        fig, ax = plt.subplots(figsize = (4, 5))
        sns.barplot(data = df1, x = 'variable', y = 'value', hue = 'index', palette = 'Paired')
        ax.set_ylim(0, 100)
        ax.tick_params(axis = 'x', rotation = 90)
    else:
        fig, ax = plt.subplots()
        sns.barplot(data = df1, y = 'variable', x = 'value', hue = 'index', palette = 'Paired')
        ax.set_xlim(0, 100)
    # ax.tick_params(axis = 'x', rotation = 90)
    ax.set(xlabel = '', ylabel = '')
    ax.set_title(label = title, pad = 20)
#     ax.set_xlim(0, 45)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon = False)
    sns.despine()
    plt.show()
    return(H2AZ_overlap)


### venn diagram for peaks / bed files
def venn_mpl(a, b, colors=None, 
             outfn="out.pdf", labels=None,text=None,
             dpi=300, figsize = (5, 4)):
    """
    *a*, *b*, and *c* are filenames to BED-like files.
    *colors* is a list of matplotlib colors for the Venn diagram circles.
    *outfn* is the resulting output file.  This is passed directly to
    fig.savefig(), so you can supply extensions of .png, .pdf, or whatever your
    matplotlib installation supports.
    *labels* is a list of labels to use for each of the files; by default the
    labels are ['a','b','c']
    
    *dpi* is the dpi setting passed to matplotlib savefig
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle
    except ImportError:
        sys.stderr.write(
            "matplotlib is required to make a Venn diagram with %s\n"
            % os.path.basename(sys.argv[0])
        )
        sys.exit(1)

    a = pybedtools.BedTool(a)
    b = pybedtools.BedTool(b)

    if colors is None:
        colors = ["r", "b"] #if c is not None else ["r", "b"]

    radius = 1.0
    center = 0.0
    offset = radius / 2
    s = sum([a.count(), b.count()]) #if c is None else sum([a.count(), b.count(), c.count()])
    
    if labels is None:
        labels = ["a", "b"] #if c is not None else ["a", "b"]
    aradius = radius
    bradius = aradius * b.count() / a.count()
    ab = (a + b).count()
    
    Oa = ab * aradius / a.count()
    Ob = ab * bradius / b.count()
    
    aoffset = aradius - Oa
    boffset = bradius - Ob 
    
    circle_a = Circle(
        xy=(center - aoffset, center),
        radius=aradius,
        edgecolor=colors[0],
        label=labels[0],
    )
    
    circle_b = Circle(
        xy=(center + boffset, center),
        radius=bradius,
        edgecolor=colors[1],
        label=labels[1],
    )
    
    fig = plt.figure(facecolor="w", figsize = figsize)
    ax = fig.add_subplot(111)

    for circle in (circle_a, circle_b):
        circle.set_facecolor("none")
        circle.set_linewidth(3)
        ax.add_patch(circle)


    ax.axis("tight")
    ax.axis("equal")
    ax.set_axis_off()

    kwargs = dict(horizontalalignment="center")

    
    atextset = aradius - 0.5*Oa
    btextset = bradius - 0.5*Ob 
    if text is True:
        # Unique to A
        t1 = ax.text(center - atextset, center, str((a - b).count()), **kwargs)

        # Unique to B
        t2 = ax.text(center + btextset, center, str((b - a).count()), **kwargs)

        t3 = ax.text(
                center, center, str((a + b).count()), **kwargs
            )
        adjust_text([t1, t2, t3], arrowprops=dict(arrowstyle="-", lw=0.5), save_steps = False, **kwargs)
    else:
        print([str((a - b).count()), str((a + b).count()), str((b - a).count())])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 8)
    plt.tight_layout()
    fig.savefig(outfn, dpi=dpi) if outfn else None
    plt.show()
    plt.close(fig)

def _bw_peaks_(bw_dict, peak_site1, peak_site2, xlabel=(), d = 2000, colmap = {}, ylim = (), figsize = (6.5, 3), ylabel = None):
    # H2AZ signal at cl up and down genes
    peak1_values = []
    for path in bw_dict:
        bw = pyBigWig.open(bw_dict[path])
        for i, line in peak_site1.iterrows():
            tmp = line.values.tolist()
            c = tmp[1] + int((tmp[2]-tmp[1])/2)
            v = bw.values(tmp[0], c-d, c+d)
            peak1_values.append(tmp[:3]+[v, path])
    
    ddf1 = pd.DataFrame(peak1_values)
    ddf1 = ddf1.groupby(4).apply(lambda d: pd.DataFrame(d[3].tolist()).mean())
    ddf1['cond'] = [x.replace('_rep1', '').replace('_rep2', '').replace('_rep3', '') for x in ddf1.index.tolist()]
    ddf1 = ddf1.groupby('cond').mean()
    if isinstance(peak_site2, pd.DataFrame):
        peak2_values = []
        for path in bw_dict:
            bw = pyBigWig.open(bw_dict[path])
            for i, line in peak_site2.iterrows():
                tmp = line.values.tolist()
                c = tmp[1] + int((tmp[2]-tmp[1])/2)
                v = bw.values(tmp[0], c-d, c+d)
                peak2_values.append(tmp[:3]+[v, path])

        ddf2 = pd.DataFrame(peak2_values)
        ddf2 = ddf2.groupby(4).apply(lambda d: pd.DataFrame(d[3].tolist()).mean())
        ddf2['cond'] = [x.replace('_rep1', '').replace('_rep2', '').replace('_rep3', '') for x in ddf2.index.tolist()]
        ddf2 = ddf2.groupby('cond').mean()
        
        fig, ax = plt.subplots(ncols = 2, sharey = True, figsize = figsize)
        for x in ddf1.index.tolist():
            ax[0].plot(range(0, ddf1.shape[1]), 
                       ddf1.loc[x,:], label = x, linewidth = 1.5, color=colmap.get(x, None) if len(colmap) > 0 else None)
        ax[0].set(xlabel = 'Peak set1' if len(xlabel) != 2 else xlabel[0], 
                  ylabel = 'ChIP-seq signal' if not ylabel else ylabel)
        ax[0].set_xticks([0, .5*ddf1.shape[1], ddf1.shape[1]], ['-'+str(d/1000)+'k', 0, '+'+str(d/1000)+'k'])
        
        for x in ddf2.index.tolist():
            ax[1].plot(range(0, ddf2.shape[1]), 
                       ddf2.loc[x,:], label = x, linewidth = 1.5, color=colmap.get(x, None) if len(colmap) > 0 else None)
        ax[1].set(xlabel = 'Peak set2' if len(xlabel) != 2 else xlabel[1], ylabel = '')
        ax[1].set_xticks([0, .5*ddf2.shape[1], ddf2.shape[1]], ['-'+str(d/1000)+'k', 0, '+'+str(d/1000)+'k'])
        ax[1].legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon = False)
        sns.despine()
        ax[0].set_ylim(ylim[0], ylim[1]) if len(ylim) > 0 else None
        plt.tight_layout()
        plt.show()
        return(peak1_values, peak2_values)
    else:
        fig, ax = plt.subplots(figsize = figsize)
        for x in ddf1.index.tolist():
            ax.plot(range(0, ddf1.shape[1]), 
                       ddf1.loc[x,:], label = x, linewidth = 1.5, color=colmap.get(x, None) if len(colmap) > 0 else None)
        ax.set(xlabel = 'Peak' if not xlabel else xlabel, 
                  ylabel = 'ChIP-seq signal' if not ylabel else ylabel)
        ax.set_xticks([0, .5*ddf1.shape[1], ddf1.shape[1]], ['-'+str(d/1000)+'k', 0, '+'+str(d/1000)+'k'])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon = False)
        sns.despine()
        ax.set_ylim(ylim[0], ylim[1]) if len(ylim) > 0 else None
        plt.tight_layout()
        plt.show()
        return(peak1_values)


# ### looping analysis by Peakahu
# 

# In[ ]:


## load loop probability for each dot in pooled data

shNT_plusCL_all = collections.defaultdict()
with open('./merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/shNT_plusCL_mapq5_merge.pairs_5000.bedpe') as source:
    for line in source:
        p = line.rstrip().split('\t')
        key = (p[0], int(p[1]), int(p[2]), p[3], int(p[4]), int(p[5]))
        scores = list(map(float, p[-2:])) 
        shNT_plusCL_all[key] = scores

shNT_minusCL_all = collections.defaultdict()
with open('./merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/shNT_minusCL_mapq5_merge.pairs_5000.bedpe') as source:
    for line in source:
        p = line.rstrip().split('\t')
        key = (p[0], int(p[1]), int(p[2]), p[3], int(p[4]), int(p[5]))
        scores = list(map(float, p[-2:])) 
        shNT_minusCL_all[key] = scores
        
KD_plusCL_all = collections.defaultdict()
with open('./merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/KD_plusCL_mapq5_merge.pairs_5000.bedpe') as source:
    for line in source:
        p = line.rstrip().split('\t')
        key = (p[0], int(p[1]), int(p[2]), p[3], int(p[4]), int(p[5]))
        scores = list(map(float, p[-2:])) 
        KD_plusCL_all[key] = scores
#### load loop probability for each dot in each replicate
shNT_plusCL_rep1 = collections.defaultdict()
with open('./merge_pairs_mapq5/remove100bp/replicate/loops/shNT_plusCL_rep1.pairs.mapq5.pairs_5000.bedpe') as source:
    for line in source:
        p = line.rstrip().split('\t')
        key = (p[0], int(p[1]), int(p[2]), p[3], int(p[4]), int(p[5]))
        scores = list(map(float, p[-2:])) 
        shNT_plusCL_rep1[key] = scores

shNT_plusCL_rep2 = collections.defaultdict()
with open('./merge_pairs_mapq5/remove100bp/replicate/loops/shNT_plusCL_rep2.pairs.mapq5.pairs_5000.bedpe') as source:
    for line in source:
        p = line.rstrip().split('\t')
        key = (p[0], int(p[1]), int(p[2]), p[3], int(p[4]), int(p[5]))
        scores = list(map(float, p[-2:])) 
        shNT_plusCL_rep2[key] = scores
shNT_plusCL_rep3 = collections.defaultdict()
with open('./merge_pairs_mapq5/remove100bp/replicate/loops/shNT_plusCL_rep3.pairs.mapq5.pairs_5000.bedpe') as source:
    for line in source:
        p = line.rstrip().split('\t')
        key = (p[0], int(p[1]), int(p[2]), p[3], int(p[4]), int(p[5]))
        scores = list(map(float, p[-2:])) 
        shNT_plusCL_rep3[key] = scores
shNT_plusCL_rep4 = collections.defaultdict()
with open('./merge_pairs_mapq5/remove100bp/replicate/loops/shNT_plusCL_rep4.pairs.mapq5.pairs_5000.bedpe') as source:
    for line in source:
        p = line.rstrip().split('\t')
        key = (p[0], int(p[1]), int(p[2]), p[3], int(p[4]), int(p[5]))
        scores = list(map(float, p[-2:])) 
        shNT_plusCL_rep4[key] = scores

shNT_minusCL_rep1 = collections.defaultdict()
with open('./merge_pairs_mapq5/remove100bp/replicate/loops/shNT_minusCL_rep1.pairs.mapq5.pairs_5000.bedpe') as source:
    for line in source:
        p = line.rstrip().split('\t')
        key = (p[0], int(p[1]), int(p[2]), p[3], int(p[4]), int(p[5]))
        scores = list(map(float, p[-2:])) 
        shNT_minusCL_rep1[key] = scores

shNT_minusCL_rep2 = collections.defaultdict()
with open('./merge_pairs_mapq5/remove100bp/replicate/loops/shNT_minusCL_rep2.pairs.mapq5.pairs_5000.bedpe') as source:
    for line in source:
        p = line.rstrip().split('\t')
        key = (p[0], int(p[1]), int(p[2]), p[3], int(p[4]), int(p[5]))
        scores = list(map(float, p[-2:])) 
        shNT_minusCL_rep2[key] = scores
shNT_minusCL_rep3 = collections.defaultdict()
with open('./merge_pairs_mapq5/remove100bp/replicate/loops/shNT_minusCL_rep3.pairs.mapq5.pairs_5000.bedpe') as source:
    for line in source:
        p = line.rstrip().split('\t')
        key = (p[0], int(p[1]), int(p[2]), p[3], int(p[4]), int(p[5]))
        scores = list(map(float, p[-2:])) 
        shNT_minusCL_rep3[key] = scores
shNT_minusCL_rep4 = collections.defaultdict()
with open('./merge_pairs_mapq5/remove100bp/replicate/loops/shNT_minusCL_rep4.pairs.mapq5.pairs_5000.bedpe') as source:
    for line in source:
        p = line.rstrip().split('\t')
        key = (p[0], int(p[1]), int(p[2]), p[3], int(p[4]), int(p[5]))
        scores = list(map(float, p[-2:])) 
        shNT_minusCL_rep4[key] = scores

### identify confident dots by restricting loop probability > 0.95
prob_cut = 0.95
shNT_plusCL_good = [x for x in shNT_plusCL_all if shNT_plusCL_all[x][0] > prob_cut]
shNT_minusCL_good = [x for x in shNT_minusCL_all if shNT_minusCL_all[x][0] > prob_cut]

shNT_union_good = collections.defaultdict()
for x in shNT_plusCL_good + shNT_minusCL_good:
    if x in shNT_union_good:
        continue
    v1 = shNT_plusCL_all[x] if x in shNT_plusCL_all else [0, 0]
    v2 = shNT_minusCL_all[x] if x in shNT_minusCL_all else [0, 0]
    v3 = KD_plusCL_all[x] if x in KD_plusCL_all else [0, 0]
    v11 = shNT_plusCL_rep1[x] if x in shNT_plusCL_rep1 else [0, 0]
    v12 = shNT_plusCL_rep2[x] if x in shNT_plusCL_rep2 else [0, 0]
    v13 = shNT_plusCL_rep3[x] if x in shNT_plusCL_rep3 else [0, 0]
    v14 = shNT_plusCL_rep4[x] if x in shNT_plusCL_rep4 else [0, 0]
    v21 = shNT_minusCL_rep1[x] if x in shNT_minusCL_rep1 else [0, 0]
    v22 = shNT_minusCL_rep2[x] if x in shNT_minusCL_rep2 else [0, 0]
    v23 = shNT_minusCL_rep3[x] if x in shNT_minusCL_rep3 else [0, 0]
    v24 = shNT_minusCL_rep4[x] if x in shNT_minusCL_rep4 else [0, 0]
    shNT_union_good[x] = v1+v2+v3+v11+v12+v13+v14+v21+v22+v23+v24

## generate a dataframe, include all the information, all scores in pooled and replicated data
loop_dots = []
for x in list(shNT_union_good.keys()):
    loop_dots.append(list(x)+shNT_union_good[x])
loop_dots = pd.DataFrame(loop_dots, columns = ['chrom1', 'start1', 'end1',
                                              'chrom2', 'start2', 'end2',
                                              'prob_shNT_plusCL', 'score_shNT_plusCL', 
                                              'prob_shNT_minusCL', 'score_shNT_minusCL', 
                                              'prob_KD_plusCL', 'score_KD_plusCL', 
                                              'prob_shNT_plusCL_rep1', 'score_shNT_plusCL_rep1', 
                                              'prob_shNT_plusCL_rep2', 'score_shNT_plusCL_rep2', 
                                              'prob_shNT_plusCL_rep3', 'score_shNT_plusCL_rep3', 
                                              'prob_shNT_plusCL_rep4', 'score_shNT_plusCL_rep4', 
                                              'prob_shNT_minusCL_rep1', 'score_shNT_minusCL_rep1', 
                                              'prob_shNT_minusCL_rep2', 'score_shNT_minusCL_rep2', 
                                               'prob_shNT_minusCL_rep3', 'score_shNT_minusCL_rep3', 
                                              'prob_shNT_minusCL_rep4', 'score_shNT_minusCL_rep4'])

loop_dots.to_csv('loop_dots.csv', index = None)


# In[11]:


### identify differential dots by 2-fold change thresthold
loop_dots = pd.read_csv('loop_dots.csv')
cord = ['chrom1','start1','end1','chrom2','start2','end2']

loop_dots['log2fc'] = loop_dots.apply(lambda row: np.log2(row['prob_shNT_plusCL'])-np.log2(np.clip(row['prob_shNT_minusCL'], a_min=0.0001, a_max = 1)), axis = 1)
loop_dots['label'] = loop_dots[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)

loop_dots.index = loop_dots['label'].tolist()
## define up and down dots using 2-fold as cutoff
up_dots = loop_dots.query('log2fc > 1')
up_dots['Type'] = 'plusCL specific'
down_dots = loop_dots.query('log2fc < -1')
down_dots['Type'] = 'minusCL specific'
stable_dots = loop_dots.query('log2fc >= -1 and log2fc <= 1')
stable_dots['Type'] = 'No'

# ## output dots to pool them as loops in Peakachu
# cord = ['chrom1','start1','end1','chrom2','start2','end2']
# up_dots[cord+['prob_shNT_plusCL','score_shNT_plusCL']].to_csv('./merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/shNT_plusCL_specific_dots_5000.bedpe',
#                                                              sep = '\t', index = None, header = None)
# down_dots[cord+['prob_shNT_minusCL','score_shNT_minusCL']].to_csv('./merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/shNT_minusCL_specific_dots_5000.bedpe',
#                                                              sep = '\t', index = None, header = None)
# stable_dots[cord+['prob_shNT_minusCL','score_shNT_minusCL']].to_csv('./merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/shNT_plusCL_minusCL_stable_dots_5000.bedpe',
#                                                              sep = '\t', index = None, header = None)


# In[25]:


stable_dots[cord+['prob_shNT_plusCL','score_shNT_plusCL']].to_csv('./merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/shNT_plusCL_minusCL_stable_dots_5000_plusCLprob.bedpe',
                                                             sep = '\t', index = None, header = None)
   


# In[26]:


stable_dots_copy = stable_dots[cord].copy()
stable_dots_copy['prob'] = stable_dots[['prob_shNT_plusCL', 'prob_shNT_minusCL']].T.max().tolist()
stable_dots_copy['score'] = stable_dots[['score_shNT_plusCL', 'score_shNT_minusCL']].T.max().tolist()

stable_dots_copy.to_csv('./merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/shNT_plusCL_minusCL_stable_dots_5000_maxProb.bedpe',
                                                             sep = '\t', index = None, header = None)
                                                                    


# In[12]:


### read in the Peakachu pooled loops
cord = ['chrom1','start1','end1','chrom2','start2','end2']

up_loop = pd.read_csv('merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/shNT_plusCL_specific_loops_5000.pool.bedpe',
                      sep = '\t', header = None)
up_loop.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'prob', 'score']
up_loop['size1'] = up_loop['end1'] - up_loop['start1']
up_loop['size2'] = up_loop['end2'] - up_loop['start2']
up_loop['distance'] = up_loop['start2'] - up_loop['start1']

###
down_loop = pd.read_csv('merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/shNT_minusCL_specific_loops_5000.pool.bedpe',
                       sep = '\t', header = None)
down_loop.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'prob', 'score']
down_loop['size1'] = down_loop['end1'] - down_loop['start1']
down_loop['size2'] = down_loop['end2'] - down_loop['start2']
down_loop['distance'] = down_loop['start2'] - down_loop['start1']
####
stable_loop = pd.read_csv('merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/shNT_plusCL_minusCL_stable_loops_5000.pool.bedpe',
                       sep = '\t', header = None)
stable_loop.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'prob', 'score']
stable_loop['size1'] = stable_loop['end1'] - stable_loop['start1']
stable_loop['size2'] = stable_loop['end2'] - stable_loop['start2']
stable_loop['distance'] = stable_loop['start2'] - stable_loop['start1']

up_loop['label'] = up_loop[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)
down_loop['label'] = down_loop[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)
stable_loop['label'] = stable_loop[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)

## show the interaction frequency of up and down loops
##
df1 = loop_dots[loop_dots['label'].isin(up_loop['label'])]
df1['Type'] = 'Up'
df2 = loop_dots[loop_dots['label'].isin(down_loop['label'])]
df2['Type'] = 'Down'
df3 = loop_dots[loop_dots['label'].isin(stable_loop['label'])]
df3['Type'] = 'No'
union_loops = pd.concat([df1, df2, df3])
union_loops['loop_length'] = np.abs(union_loops['start1'] - union_loops['start2'])


# In[30]:


## compare stable loop for different pool
stable_plusCLprob = pd.read_csv('merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/shNT_plusCL_minusCL_stable_dots_5000_plusCLprob.pool.bedpe',
                               sep = '\t', header = None)
stable_plusCLprob.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'prob', 'score']
stable_plusCLprob['label'] = stable_plusCLprob[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)

stable_maxProb = pd.read_csv('merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/shNT_plusCL_minusCL_stable_dots_5000_maxProb.pool.bedpe',
                               sep = '\t', header = None)
stable_maxProb.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'prob', 'score']
stable_maxProb['label'] = stable_maxProb[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)


# In[31]:


len(set(stable_loop['label'].tolist()) & set(stable_plusCLprob['label'].tolist()))


# In[32]:


len(set(stable_maxProb['label'].tolist()) & set(stable_plusCLprob['label'].tolist()))


# In[8]:


## quantile-normalization and show the loop scores for diff loops
plot_df = pd.DataFrame(quantile_norm(union_loops[['score_shNT_plusCL', 'score_shNT_minusCL']]),
            columns = ['shNT_plusCL', 'shNT_minusCL'])
plot_df['Type']=union_loops['Type'].tolist()
pplot_df = plot_df.melt(id_vars = ['Type'])
pplot_df['variable'] = pd.Categorical(pplot_df['variable'], ['shNT_minusCL', 'shNT_plusCL'])
pplot_df['Type'] = pd.Categorical(pplot_df['Type'], ['Up', 'Down', 'No'])

## plot
fig, ax = plt.subplots(figsize = (5, 3))
sns.violinplot(data = pplot_df, x = 'Type', y = 'value',
               hue = 'variable', ax = ax, scale = 'width', 
               palette={'shNT_minusCL':'grey', 
                        'shNT_plusCL':plt.cm.get_cmap('tab10')(1)})
ax.set(xlabel = '', ylabel = 'Contact Frequency')
ax.legend(title = '', loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.show()
fig.savefig('Figures/diff_loop_violin.pdf')


# In[9]:


plot_df.groupby('Type').apply(lambda df: ttest_ind(df['shNT_plusCL'], df['shNT_minusCL']))


# ### APA plot

# In[10]:


r = 10000
cond1 = 'shNT_plusCL'
cond2 = 'shNT_minusCL'
cond3 = 'KD_plusCL'
## load cool file for APA plots
cool1_path = './merge_pairs_mapq5/remove100bp/cool/'+cond1+'_mapq5_merge.pairs_%s.cool'%r
cool2_path = './merge_pairs_mapq5/remove100bp/cool/'+cond2+'_mapq5_merge.pairs_%s.cool'%r
cool3_path = './merge_pairs_mapq5/remove100bp/cool/'+cond3+'_mapq5_merge.pairs_%s.cool'%r

clr1 = cooler.Cooler(cool1_path)
clr2 = cooler.Cooler(cool2_path)
clr3 = cooler.Cooler(cool3_path)

# Use bioframe to fetch the genomic features from the UCSC.
# mm10_chromsizes = bioframe.fetch_chromsizes('mm10')
# mm10_cens = bioframe.fetch_centromeres('mm10')
# mm10_arms = bioframe.make_chromarms(mm10_chromsizes, mm10_cens)
d = pk.load(open('mm10_bioframe.pk', 'rb'))
mm10_chromsizes=d['mm10_chromsizes']
mm10_cens=d['mm10_cens']
mm10_arms=d['mm10_arms']
# Select only chromosomes that are present in the cooler. 
# This step is typically not required! we call it only because the test data are reduced. 
mm10_arms = mm10_arms.set_index("chrom").loc[clr1.chromnames].reset_index()
# call this to automaticly assign names to chromosomal arms:
mm10_arms = bioframe.make_viewframe(mm10_arms)

### compute expected loop scores as background
# expected1 = expected_cis(
#             clr1,
#             ignore_diags=0,
#             view_df=mm10_arms,
#             chunksize=100000)
# out = open('expected1.pk', 'wb')
# pk.dump(expected1, out)
# out.close()

# expected2 = expected_cis(
#             clr2,
#             ignore_diags=0,
#             view_df=mm10_arms,
#             chunksize=1000000)
# out = open('expected2.pk', 'wb')
# pk.dump(expected2, out)
# out.close()

# expected3 = expected_cis(
#             clr3,
#             ignore_diags=0,
#             view_df=mm10_arms,
#             chunksize=1000000)
# out = open('expected3.pk', 'wb')
# pk.dump(expected3, out)
# out.close()

out = open('expected1_10k.pk', 'rb')
expected1 = pk.load(out)
out.close()

out = open('expected2_10k.pk', 'rb')
expected2 = pk.load(out)
out.close()

out = open('expected3_10k.pk', 'rb')
expected3 = pk.load(out)
out.close()


# In[ ]:


### plot APA for up loop in three conditions
pup1 = coolpup.pileup(clr1, up_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected1,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)


plotpup.plot(pup1,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, 
#              vmax = 4, vmin = 1.3,
             height=1.5)
plt.show()

pup2 = coolpup.pileup(clr2, up_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected2,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup2,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 4, vmin = 1.3,
             height=1.5)
plt.show()

pup3 = coolpup.pileup(clr3, up_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected2,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup3,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 4, vmin = 1.3,
             height=1.5)
plt.show()


# In[ ]:


### plot APA for down loop in three conditions

pup1 = coolpup.pileup(clr1, down_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected1,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)


plotpup.plot(pup1,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, 
             vmax = 4.1, vmin = 1.3,
             height=1.5)
plt.show()

pup2 = coolpup.pileup(clr2, down_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected2,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup2,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
#              vmax = 3, vmin = 1.3,
             height=1.5)
plt.show()

pup3 = coolpup.pileup(clr3, down_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected2,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup3,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 4.1, vmin = 1.3,
             height=1.5)
plt.show()


# In[ ]:


### plot APA for no change loop in three conditions

pup1 = coolpup.pileup(clr1, stable_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected1,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)


plotpup.plot(pup1,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, 
#              vmax = 4, vmin = 1.3,
             height=1.5)
plt.show()

pup2 = coolpup.pileup(clr2, stable_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected2,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup2,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 6.2, vmin = 1.5,
             height=1.5)
plt.show()

pup3 = coolpup.pileup(clr3, stable_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected2,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup3,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 6.2, vmin = 1.5,
             height=1.5)
plt.show()


# ### integrate loop change with gene expression
# 

# In[6]:


## load differential expression result and expression TPM matrix
out = open('20221212_RNA-seq_diffexp_repCorrect.pk', 'rb')
diffexp = pk.load(out)
out.close()

exp_tpm = pd.read_csv('tpm_matrix.csv', index_col = 0)


# In[ ]:


# ### diff loop with gene expression change in boxplot
# up_loop_anchor = up_loop[['chrom1', 'start1', 'end1']].values.tolist()+up_loop[['chrom2', 'start2', 'end2']].values.tolist()
# up_loop_anchor = pd.DataFrame(up_loop_anchor, columns = ['chr', 'start', 'end'])
# down_loop_anchor =down_loop[['chrom1', 'start1', 'end1']].values.tolist()+down_loop[['chrom2', 'start2', 'end2']].values.tolist()
# down_loop_anchor = pd.DataFrame(down_loop_anchor, columns = ['chr', 'start', 'end'])
# stable_loop_anchor =stable_loop[['chrom1', 'start1', 'end1']].values.tolist()+stable_loop[['chrom2', 'start2', 'end2']].values.tolist()
# stable_loop_anchor = pd.DataFrame(stable_loop_anchor, columns = ['chr', 'start', 'end'])

# pdf = PdfPages('Figures/all_loop_mRNA_fc.pdf')
# plot_df = _plot_loop_exp_(mm10_tss_ann, up_loop_anchor, down_loop_anchor, 
#                 stable_loop_anchor, diffexp['shNT_plusCL_vs_minusCL'], col='stat', d = 5000, figsize = (4,4.5),
#                          pdf = pdf)
# pdf.close()


# In[13]:


### classify loops into E-P, P-P, P-O, and check them with gene expression change in boxplot

## EP
q1 = '(r1_promoter == True and r2_ELS == True and r2_K27ac == True) or (r2_promoter == True and r1_ELS == True and r1_K27ac == True)'
## PP
q2 = '(r1_promoter == True and r2_promoter == True)'
## P-other
q3 = '(r1_promoter == True and (r2_ELS != True or r2_K27ac != True) and r2_promoter != True) or (r2_promoter == True and (r1_ELS != True or r1_K27ac != True) and r1_promoter != True)'

d = 5000

## promoter define from mm10 genome
mm10_promoter = mm10_tss_ann.copy()
mm10_promoter['start'] = mm10_tss_ann['start'] - d
mm10_promoter['end'] = mm10_tss_ann['end'] + d
mm10_promoter = mm10_promoter[mm10_promoter['start']>0]

mm10_promoter_bed = pybedtools.BedTool.from_dataframe(mm10_promoter)

## define enhancer
mm10_cRE = pd.read_csv('./mm10-cCREs.bed', sep = '\t', header = None)
mm10_ELS = mm10_cRE[mm10_cRE[5].str.contains('ELS')]
mm10_ELS_bed = pybedtools.BedTool.from_dataframe(mm10_ELS)
mm10_ELS_bed = mm10_ELS_bed.intersect(mm10_promoter_bed, v = True, wa = True)
## our data 
k27ac_plusCL_rep1 = pd.read_csv('./histone/H3K27ac_plusCL.rep1_peaks.narrowPeak', sep = '\t', header = None)
k27ac_plusCL_rep2 = pd.read_csv('./histone/H3K27ac_plusCL.rep2_peaks.narrowPeak', sep = '\t', header = None)
k27ac_minusCL_rep1 = pd.read_csv('./histone/H3K27ac_minusCL.rep1_peaks.narrowPeak', sep = '\t', header = None)
k27ac_minusCL_rep2 = pd.read_csv('./histone/H3K27ac_minusCL.rep2_peaks.narrowPeak', sep = '\t', header = None)
k27ac_union_peaks = pd.concat([k27ac_plusCL_rep1, k27ac_plusCL_rep2, k27ac_minusCL_rep1, k27ac_minusCL_rep2])
k27ac_union_peaks_bed = pybedtools.BedTool.from_dataframe(k27ac_union_peaks.iloc[:,0:3])
k27ac_union_peaks_bed = k27ac_union_peaks_bed.intersect(mm10_promoter_bed, v = True, wa = True)

ELS_k27ac_union_peaks_bed = mm10_ELS_bed.intersect(k27ac_union_peaks_bed, wa = True)


### annotate 
# union_loops_genes = pd.concat([up_dots, down_dots, stable_dots]) 
union_loops_genes = union_loops.copy()
score_qnorm = pd.DataFrame(quantile_norm(union_loops[['score_shNT_plusCL', 'score_shNT_minusCL']]),
            columns = ['shNT_plusCL', 'shNT_minusCL'])
union_loops_genes['score_shNT_plusCL'] = score_qnorm['shNT_plusCL'].tolist()
union_loops_genes['score_shNT_minusCL'] = score_qnorm['shNT_minusCL'].tolist()

union_loops_genes['r1'] = union_loops_genes['chrom1']+':'+union_loops_genes['start1'].astype('str')+'-'+union_loops_genes['end1'].astype('str')
union_loops_genes['r2'] = union_loops_genes['chrom2']+':'+union_loops_genes['start2'].astype('str')+'-'+union_loops_genes['end2'].astype('str')
## overlap with promoter
r1 = pybedtools.BedTool.from_dataframe(union_loops_genes[['chrom1','start1','end1', 'r1']])
r2 = pybedtools.BedTool.from_dataframe(union_loops_genes[['chrom2','start2','end2', 'r2']])
r1_promoter = r1.intersect(mm10_promoter_bed, wao = True).to_dataframe().drop_duplicates()
r2_promoter = r2.intersect(mm10_promoter_bed, wao = True).to_dataframe().drop_duplicates()

## overlap with ELS
r1_ELS = r1.intersect(mm10_ELS_bed, wao = True).to_dataframe().drop_duplicates()
r2_ELS = r2.intersect(mm10_ELS_bed, wao = True).to_dataframe().drop_duplicates()

## overlap with K27ac
r1_K27ac = r1.intersect(k27ac_union_peaks_bed, wao = True).to_dataframe().drop_duplicates()
r2_K27ac = r2.intersect(k27ac_union_peaks_bed, wao = True).to_dataframe().drop_duplicates()

## overlap with ELS_K27ac
r1_ELS_K27ac = r1.intersect(ELS_k27ac_union_peaks_bed, wao = True).to_dataframe().drop_duplicates()
r2_ELS_K27ac = r2.intersect(ELS_k27ac_union_peaks_bed, wao = True).to_dataframe().drop_duplicates()

union_loops_genes['r1_promoter'] = union_loops_genes['r1'].isin(r1_promoter.query('score != "."')['name'])
union_loops_genes['r2_promoter'] = union_loops_genes['r2'].isin(r2_promoter.query('score != "."')['name'])
union_loops_genes['r1_ELS'] = union_loops_genes['r1'].isin(r1_ELS.query('score != "."')['name'])
union_loops_genes['r2_ELS'] = union_loops_genes['r2'].isin(r2_ELS.query('score != "."')['name'])
union_loops_genes['r1_K27ac'] = union_loops_genes['r1'].isin(r1_K27ac.query('score != "."')['name'])
union_loops_genes['r2_K27ac'] = union_loops_genes['r2'].isin(r2_K27ac.query('score != "."')['name'])
union_loops_genes['r1_ELS_K27ac'] = union_loops_genes['r1'].isin(r2_ELS_K27ac.query('score != "."')['name'])
union_loops_genes['r2_ELS_K27ac'] = union_loops_genes['r2'].isin(r2_ELS_K27ac.query('score != "."')['name'])


r1_genes = r1_promoter.query('score != "."').groupby('name').apply(lambda df: ';'.join(set([x.split(':')[-1] for x in df['thickEnd'].unique().tolist()]))).to_dict()
r2_genes = r2_promoter.query('score != "."').groupby('name').apply(lambda df: ';'.join(set([x.split(':')[-1] for x in df['thickEnd'].unique().tolist()]))).to_dict()


union_loops_genes['r1_gene'] = [r1_genes.get(x, np.nan) for x in union_loops_genes['r1'].tolist()]
union_loops_genes['r2_gene'] = [r2_genes.get(x, np.nan) for x in union_loops_genes['r2'].tolist()]

union_loops_genes['r1_mRNA_fc'] = [diffexp['shNT_plusCL_vs_minusCL'].reindex(x.split(';'))['log2FoldChange'].mean() if not pd.isna(x) else np.nan for x in union_loops_genes['r1_gene'].tolist()]
union_loops_genes['r2_mRNA_fc'] = [diffexp['shNT_plusCL_vs_minusCL'].reindex(x.split(';'))['log2FoldChange'].mean() if not pd.isna(x) else np.nan for x in union_loops_genes['r2_gene'].tolist()]

## define EP, PP, and PO
EP_df = union_loops_genes.query(q1)

PP_df = union_loops_genes.query(q2)

PO_df = union_loops_genes.query(q3)
OO_df = union_loops_genes[~union_loops_genes['label'].isin(EP_df['label'].tolist()+PP_df['label'].tolist()+PO_df['label'].tolist())]

## calculate scores
EP_df['score_delta'] = EP_df['score_shNT_plusCL'] - EP_df['score_shNT_minusCL']
EP_df['mRNA_fc'] = EP_df[['r1_mRNA_fc', 'r2_mRNA_fc']].apply(lambda row: row.dropna().mean(), axis = 1)

PP_df['score_delta'] = PP_df['score_shNT_plusCL'] - PP_df['score_shNT_minusCL']
PP_df['mRNA_fc'] = PP_df[['r1_mRNA_fc', 'r2_mRNA_fc']].apply(lambda row: row.dropna().mean(), axis = 1)

PO_df['score_delta'] = PO_df['score_shNT_plusCL'] - PO_df['score_shNT_minusCL']
PO_df['mRNA_fc'] = PO_df[['r1_mRNA_fc', 'r2_mRNA_fc']].apply(lambda row: row.dropna().mean(), axis = 1)

### unique gene and unique loop for each row, convenient for visulization 
EP_ddf = _get_values_(EP_df)
PP_ddf = _get_values_(PP_df)
PO_ddf = _get_values_(PO_df)


# In[35]:


### classify different loop types for up, down, and no change loops 
up_EP_loop, up_PP_loop, up_PO_loop, up_OO_loop = _define_EP_PP_(up_dots[up_dots['label'].isin(up_loop['label'])])
up_EP_loop['group'] = 'E-P'
up_PP_loop['group'] = 'P-P'
up_PO_loop['group'] = 'P-O'
up_OO_loop['group'] = 'O-O'
up_loop_union = pd.concat([up_EP_loop, up_PP_loop, up_PO_loop, up_OO_loop])

down_EP_loop, down_PP_loop, down_PO_loop, down_OO_loop = _define_EP_PP_(down_dots[down_dots['label'].isin(down_loop['label'])])
down_EP_loop['group'] = 'E-P'
down_PP_loop['group'] = 'P-P'
down_PO_loop['group'] = 'P-O'
down_OO_loop['group'] = 'O-O'
down_loop_union = pd.concat([down_EP_loop, down_PP_loop, down_PO_loop, down_OO_loop])


# In[ ]:


## generate gene table for E-P, P-P, and P-O
tab_res = []
for i, line in EP_ddf.iterrows():
    loop = line['label']
    delta = line['score_delta']
    lowess_value = line['lowess']
    for g in line[['r1_gene', 'r2_gene']].dropna().tolist():
        tab_res.append([loop, delta, 'E-P', g, lowess_value])
        
for i, line in PP_ddf.iterrows():
    loop = line['label']
    delta = line['score_delta']
    lowess_value = line['lowess']
    for g in line[['r1_gene', 'r2_gene']].dropna().tolist():
        tab_res.append([loop, delta, 'P-P', g, lowess_value])
        
for i, line in PO_ddf.iterrows():
    loop = line['label']
    delta = line['score_delta']
    lowess_value = line['lowess']
    for g in line[['r1_gene', 'r2_gene']].dropna().tolist():
        tab_res.append([loop, delta, 'P-O', g, lowess_value])
        
tab_res = pd.DataFrame(tab_res, columns = ['loop', 'delta_contact', 'loop_type', 'gene', 'gene_lowess'])
tab_res = pd.merge(tab_res, diffexp['shNT_plusCL_vs_minusCL'], left_on = 'gene', right_index = True)
tab_res = tab_res.sort_values(['delta_contact', 'log2FoldChange'], ascending = False)
## diff loop
loop_de_label = []
for x in tab_res['loop'].tolist():
    if x in up_dots.index.tolist():
        loop_de_label.append('Up')
    elif x in down_dots.index.tolist():
        loop_de_label.append('Down')
    else:
        loop_de_label.append('No')
tab_res['diff_loop'] = loop_de_label
    
tab_res.to_csv('plusCL_vs_minusCL_loop_target_gene_diffexp.tsv', sep = '\t', index = None)


# In[29]:


### plot loop delta contact vs. gene expression change
fig, ax = plt.subplots()

## align different loop type range to the same, so that can be better compared
xmax = min([max(EP_df['score_delta']), max(PP_df['score_delta']), max(PO_df['score_delta'])])
xmin = max([min(EP_ddf['score_delta']), min(PP_ddf['score_delta']), min(PO_ddf['score_delta'])])
q = 'score_delta >= @xmin and score_delta <= @xmax'
ax.plot(EP_ddf.drop_duplicates(['label', 'gene']).query(q)['score_delta'], EP_ddf.drop_duplicates(['label', 'gene']).query(q)['lowess'], label = 'E-P')
ax.plot(PP_ddf.drop_duplicates(['label', 'gene']).query(q)['score_delta'], PP_ddf.drop_duplicates(['label', 'gene']).query(q)['lowess'], label = 'P-P')
ax.plot(PO_ddf.drop_duplicates(['label', 'gene']).query(q)['score_delta'], PO_ddf.drop_duplicates(['label', 'gene']).query(q)['lowess'], label = 'P-Other')
ax.fill_between(EP_ddf.drop_duplicates(['label', 'gene']).query(q)['score_delta'],EP_ddf.drop_duplicates(['label', 'gene']).query(q)['lower_lowess'],EP_ddf.drop_duplicates(['label', 'gene']).query(q)['upper_lowess'],alpha=.33)
ax.fill_between(PP_ddf.drop_duplicates(['label', 'gene']).query(q)['score_delta'],PP_ddf.drop_duplicates(['label', 'gene']).query(q)['lower_lowess'],PP_ddf.drop_duplicates(['label', 'gene']).query(q)['upper_lowess'],alpha=.33)
ax.fill_between(PO_ddf.drop_duplicates(['label', 'gene']).query(q)['score_delta'],PO_ddf.drop_duplicates(['label', 'gene']).query(q)['lower_lowess'],PO_ddf.drop_duplicates(['label', 'gene']).query(q)['upper_lowess'],alpha=.33)
ax.set(xlabel = 'Delta Contact Frequency', ylabel = 'mRNA log2FoldChange',
      title = 'd=%s'%d)
ax.tick_params(axis = 'x', rotation = 45)
ax.legend(title = '', loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()
sns.despine()
# fig.savefig('Figures/diffloop_vs_mRNA_fc.pdf')
plt.show()



pdf = PdfPages('Figures/pie_chart_percentage_diff_category.pdf')
for x,i in [[EP_df, 'E-P'],
          [PP_df, 'P-P'], 
          [PO_df, 'P-Other'], 
          [OO_df, 'Other-Other']]:
    ### up vs down
    n1 = x[x['label'].isin(up_dots['label'])].shape[0]
    n2 = x[x['label'].isin(down_dots['label'])].shape[0]
    n3 = x[x['label'].isin(stable_dots['label'])].shape[0]

    fig, ax = plt.subplots()
    ax.pie([n1, n2, n3], labels=['Up', 'Down', 'No'], autopct='%1.1f%%', 
           colors=['#CD5555', '#1E90FF', 'lightgrey'])
    ax.set(title = i)
    plt.tight_layout()
    pdf.savefig(fig)
    plt.show()
    fig, ax = plt.subplots()
    ax.pie([n1, n2, n3], labels=['', '', ''], 
           colors=['#CD5555', '#1E90FF', 'lightgrey'])
    ax.set(title = i)
    plt.tight_layout()
    pdf.savefig(fig)
#     plt.show()
pdf.close()

# In[32]:


### read loops identified for each condition
shNT_plusCL_allLoop = pd.read_csv('./merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/shNT_plusCL_mapq5_merge.pairs_5000_loops.0.95.bedpe',
                                 sep = '\t', header = None)
shNT_plusCL_allLoop.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'prob', 'score']
shNT_plusCL_allLoop['label'] = shNT_plusCL_allLoop[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)

shNT_minusCL_allLoop = pd.read_csv('./merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/shNT_minusCL_mapq5_merge.pairs_5000_loops.0.95.bedpe',
                                  sep = '\t', header = None)
shNT_minusCL_allLoop.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'prob', 'score']
shNT_minusCL_allLoop['label'] = shNT_minusCL_allLoop[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)

KD_plusCL_allLoop = pd.read_csv('./merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/KD_plusCL_mapq5_merge.pairs_5000_loops.0.95.bedpe',
                                 sep = '\t', header = None)
KD_plusCL_allLoop.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'prob', 'score']
KD_plusCL_allLoop['label'] = KD_plusCL_allLoop[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)


### define E-P, P-P, P-O loops for each condition
shNT_plusCL_EP, shNT_plusCL_PP, shNT_plusCL_PO, shNT_plusCL_other_df = _define_EP_PP_(shNT_plusCL_allLoop)
shNT_plusCL_EP['group'] = 'EP'
shNT_plusCL_PP['group'] = 'PP'
shNT_plusCL_PO['group'] = 'PO'

shNT_minusCL_EP, shNT_minusCL_PP, shNT_minusCL_PO, shNT_minusCL_other_df = _define_EP_PP_(shNT_minusCL_allLoop)
shNT_minusCL_EP['group'] = 'EP'
shNT_minusCL_PP['group'] = 'PP'
shNT_minusCL_PO['group'] = 'PO'

KD_plusCL_EP, KD_plusCL_PP, KD_plusCL_PO, KD_plusCL_other_df = _define_EP_PP_(KD_plusCL_allLoop)
KD_plusCL_EP['group'] = 'EP'
KD_plusCL_PP['group'] = 'PP'
KD_plusCL_PO['group'] = 'PO'


# ### H2AZ KD Micro-C

# In[2]:


#### define H2AZ KD diff loops
cord = ['chrom1','start1','end1','chrom2','start2','end2']
loop_dots_copy = pd.read_csv('loop_dots.csv')

loop_dots_copy['log2fc'] = loop_dots_copy.apply(lambda row: np.log2(row['prob_KD_plusCL'])-np.log2(np.clip(row['prob_shNT_plusCL'], a_min=0.0001, a_max = 1)), axis = 1)
loop_dots_copy['delta'] = loop_dots_copy['prob_KD_plusCL'] - loop_dots_copy['prob_shNT_plusCL']
loop_dots_copy['label'] = loop_dots_copy[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)

## define up and down dots using 2-fold as cutoff
KD_up_dots = loop_dots_copy.query('log2fc > 1')
KD_up_dots['Type'] = 'KD_vs_shNT_plusCL_up'
KD_down_dots = loop_dots_copy.query('log2fc < -1')
KD_down_dots['Type'] = 'KD_vs_shNT_plusCL_down'
KD_stable_dots = loop_dots_copy.query('log2fc >= -1 and log2fc <= 1')
KD_stable_dots['Type'] = 'KD_vs_shNT_plusCL_No'

### read in the Peakachu pooled loops
KD_up_loop = pd.read_csv('merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/KD_vs_shNT_plusCL_up_loops_5000.pool.bedpe',
                      sep = '\t', header = None)
KD_up_loop.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'prob', 'score']
KD_up_loop['size1'] = KD_up_loop['end1'] - KD_up_loop['start1']
KD_up_loop['size2'] = KD_up_loop['end2'] - KD_up_loop['start2']
KD_up_loop['distance'] = KD_up_loop['start2'] - KD_up_loop['start1']

###
KD_down_loop = pd.read_csv('merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/KD_vs_shNT_plusCL_down_loops_5000.pool.bedpe',
                       sep = '\t', header = None)
KD_down_loop.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'prob', 'score']
KD_down_loop['size1'] = KD_down_loop['end1'] - KD_down_loop['start1']
KD_down_loop['size2'] = KD_down_loop['end2'] - KD_down_loop['start2']
KD_down_loop['distance'] =KD_down_loop['start2'] - KD_down_loop['start1']
####
KD_stable_loop = pd.read_csv('merge_pairs_mapq5/remove100bp/loops/Peakahu_analysis/KD_vs_shNT_plusCL_stable_loops_5000.pool.bedpe',
                       sep = '\t', header = None)
KD_stable_loop.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'prob', 'score']
KD_stable_loop['size1'] = KD_stable_loop['end1'] - KD_stable_loop['start1']
KD_stable_loop['size2'] = KD_stable_loop['end2'] - KD_stable_loop['start2']
KD_stable_loop['distance'] = KD_stable_loop['start2'] - KD_stable_loop['start1']

KD_up_loop['label'] = KD_up_loop[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)
KD_down_loop['label'] = KD_down_loop[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)
KD_stable_loop['label'] = KD_stable_loop[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)



# In[38]:


### barplot showing the number of differential loops upon H2AZ KD
print({'Up':KD_up_loop.shape[0], 'Down':KD_down_loop.shape[0], 'Stable':KD_stable_loop.shape[0]})
fig, ax = plt.subplots(figsize = (5,5))
ax.bar(['Up', 'Down', 'No'], 
       [KD_up_loop.shape[0], KD_down_loop.shape[0], KD_stable_loop.shape[0]],
      color = 'grey')
ax.set(xlabel='Loop category\n(KD_plusCL vs shNT_plusCL)',
      ylabel = 'Number of loops')
sns.despine()
plt.tight_layout()
fig.savefig('Figures/Number_diffLoop_upon_H2AZ_KD.pdf')
plt.show()


# In[39]:


### classify E-P, P-P, P-O loops for down-regulated loops upon H2AZ KD
KD_down_EP_df, KD_down_PP_df, KD_down_PO_df, KD_down_other_df = _define_EP_PP_(KD_down_loop, comp = 'plusCL_KD_vs_shNT')

KD_down_EP_df['group'] = 'E-P'
KD_down_PP_df['group'] = 'P-P'
KD_down_PO_df['group'] = 'P-O'
KD_down_other_df['group'] = 'O-O'
KD_down_loop_union = pd.concat([KD_down_EP_df, KD_down_PP_df, KD_down_PO_df, KD_down_other_df])
KD_down_loop_union.index = KD_down_loop_union['label']


# In[ ]:


### APA plot

### up by H2AZ KD
pup1 = coolpup.pileup(clr1, KD_up_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected1,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)


plotpup.plot(pup1,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, 
             vmax = 4.2, vmin = 1.2,
             height=1.5)
plt.show()
plotpup.plot(pup1,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, 
             vmax = 4.2, vmin = 1.2,
             height=1.5)
plt.show()

pup1 = coolpup.pileup(clr3, KD_up_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected3,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)


plotpup.plot(pup1,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, 
#              vmax = 3.4, vmin = 1.1,
             height=1.5)
plt.show()


plotpup.plot(pup1,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, 
#              vmax = 3.4, vmin = 1.1,
             height=1.5)
plt.show()

### down by H2AZ KD

pup2 = coolpup.pileup(clr1, KD_down_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected1,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup2,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
#              vmax = 3.4, vmin = 1.1,
             height=1.5)
plt.show()

plotpup.plot(pup2,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
#              vmax = 3.4, vmin = 1.1,
             height=1.5)
plt.show()

pup2 = coolpup.pileup(clr3, KD_down_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected3,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup2,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 4, vmin = 1.3,
             height=1.5)
plt.show()

plotpup.plot(pup2,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 4, vmin = 1.3,
             height=1.5)
plt.show()

### no change by H2AZ KD
pup3 = coolpup.pileup(clr1, KD_stable_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected1,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup3,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
#              vmax = 3.4, vmin = 1.1,
             height=1.5)
plt.show()
plotpup.plot(pup3,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
#              vmax = 3.4, vmin = 1.1,
             height=1.5)
plt.show()

pup3 = coolpup.pileup(clr3, KD_stable_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected3,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup3,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 6.6, vmin = 1.5,
             height=1.5)
plt.show()
plotpup.plot(pup3,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 6.6, vmin = 1.5,
             height=1.5)
plt.show()




# ### H2AZ-dependant loops - CL induce and H2AZ KD down loops - 2fold change

# In[42]:


## identify H2AZ-dependant loops - CL induce and KD down loops - 2fold change
tmp = pd.DataFrame(quantile_norm(loop_dots[['score_shNT_plusCL', 'score_shNT_minusCL', 'score_KD_plusCL']].copy()), 
                   columns = ['score_shNT_plusCL', 'score_shNT_minusCL', 'score_KD_plusCL'], index = loop_dots['label'].tolist())

CL_induce_KD_down = up_dots[up_dots['label'].isin(KD_down_dots['label']) & up_dots['label'].isin(up_loop['label'])]
CL_induce_KD_down['score_KD_plusCL'] = tmp.loc[CL_induce_KD_down['label'].tolist(), 'score_KD_plusCL'].tolist()
CL_induce_KD_down['score_shNT_plusCL'] = tmp.loc[CL_induce_KD_down['label'].tolist(), 'score_shNT_plusCL'].tolist()
CL_induce_KD_down['score_shNT_minusCL'] = tmp.loc[CL_induce_KD_down['label'].tolist(), 'score_shNT_minusCL'].tolist()


CL_induce_KD_down_EP, CL_induce_KD_down_PP, CL_induce_KD_down_PO, CL_induce_KD_down_OO = _define_EP_PP_(CL_induce_KD_down, comp = 'plusCL_KD_vs_shNT')

CL_induce_KD_down_tmp = pd.concat([CL_induce_KD_down_EP, CL_induce_KD_down_PP, CL_induce_KD_down_PO, CL_induce_KD_down_OO])
CL_induce_KD_down_tmp['group'] = ['E-P']*CL_induce_KD_down_EP.shape[0]+['P-P']*CL_induce_KD_down_PP.shape[0]+['P-O']*CL_induce_KD_down_PO.shape[0]+['O-O']*CL_induce_KD_down_OO.shape[0]

CL_induce_KD_down_tmp['score_delta'] = CL_induce_KD_down_tmp['score_KD_plusCL'] - CL_induce_KD_down_tmp['score_shNT_plusCL']
CL_induce_KD_down_tmp['mRNA_fc'] = CL_induce_KD_down_tmp[['r1_mRNA_fc', 'r2_mRNA_fc']].apply(lambda row: row.dropna().mean(), axis = 1)

CL_induce_KD_down_exp = _get_values_(CL_induce_KD_down_tmp, comp = 'plusCL_KD_vs_shNT')
CL_induce_KD_down_z = _get_values_(CL_induce_KD_down_tmp, comp = 'plusCL_KD_vs_shNT', value = 'stat')


# In[41]:


## APA plot for H2AZ-dependant loops - CL induce and KD down loops - 2fold change
tmp_loop = up_dots[up_dots['label'].isin(KD_down_dots['label']) & up_dots['label'].isin(up_loop['label'])]

pup2 = coolpup.pileup(clr3, tmp_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected3,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup2,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 4, vmin = 1.3,
             height=1.5)
plt.show()

plotpup.plot(pup2,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 4, vmin = 1.3,
             height=1.5)
plt.show()

## CL induce and KD down loops - 2fold change, for EP, PP, PO
##### E-P
tmp_loop = up_dots[up_dots['label'].isin(KD_down_dots['label']) & 
                   up_dots['label'].isin(up_loop_union.query('group == "E-P"')['label'])]

pup2 = coolpup.pileup(clr3, tmp_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected3,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup2,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 3.8, vmin = 1.3,
             height=1.5)
plt.show()

plotpup.plot(pup2,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 3.8, vmin = 1.3,
             height=1.5)
plt.show()
##### P-P
tmp_loop = up_dots[up_dots['label'].isin(KD_down_dots['label']) & 
                   up_dots['label'].isin(up_loop_union.query('group == "P-P"')['label'])]

pup2 = coolpup.pileup(clr3, tmp_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected3,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup2,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 3.5, vmin = 1.3,
             height=1.5)
plt.show()

plotpup.plot(pup2,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 3.5, vmin = 1.3,
             height=1.5)
plt.show()

#### P-O 
tmp_loop = up_dots[up_dots['label'].isin(KD_down_dots['label']) & 
                   up_dots['label'].isin(up_loop_union.query('group == "P-O"')['label'])]

pup2 = coolpup.pileup(clr3, tmp_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected3,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup2,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 3.9, vmin = 1.3,
             height=1.5)
plt.show()

plotpup.plot(pup2,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 3.9, vmin = 1.3,
             height=1.5)
plt.show()



# In[44]:


### analyze the loop change and differental expression upon H2AZ KD, particularly for H2AZ-dependant loops and genes
KD_stable_loop['Type'] = 'No'
KD_up_loop['Type'] = 'Up'
KD_down_loop['Type'] = 'Down'

KD_union_loop = pd.concat([KD_up_loop, KD_down_loop, KD_stable_loop])
KD_union_loop_EP, KD_union_loop_PP, KD_union_loop_PO, KD_union_loop_OO = _define_EP_PP_(KD_union_loop, comp = 'plusCL_KD_vs_shNT')
KD_union_loop_EP['group'] = 'EP'
KD_union_loop_PP['group'] = 'PP'
KD_union_loop_PO['group'] = 'PO'
KD_union_loop_OO['group'] = 'OO'
KD_union_loop = pd.concat([KD_union_loop_EP, KD_union_loop_PP, KD_union_loop_PO, KD_union_loop_OO])

tmp = pd.DataFrame(quantile_norm(loop_dots[['score_shNT_plusCL', 'score_shNT_minusCL', 'score_KD_plusCL']].copy()), 
                   columns = ['score_shNT_plusCL', 'score_shNT_minusCL', 'score_KD_plusCL'], index = loop_dots['label'].tolist())
KD_union_loop['score_KD_plusCL'] = tmp.loc[KD_union_loop['label'].tolist(), 'score_KD_plusCL'].tolist()
KD_union_loop['score_shNT_plusCL'] = tmp.loc[KD_union_loop['label'].tolist(), 'score_shNT_plusCL'].tolist()
KD_union_loop['score_delta'] = KD_union_loop['score_KD_plusCL'] - KD_union_loop['score_shNT_plusCL']

tmp = pd.DataFrame(quantile_norm(loop_dots[['prob_shNT_plusCL', 'prob_KD_plusCL']].copy()), 
                   columns = ['prob_shNT_plusCL', 'prob_KD_plusCL'], index = loop_dots['label'].tolist())
KD_union_loop['prob_KD_plusCL'] = tmp.loc[KD_union_loop['label'].tolist(), 'prob_KD_plusCL'].tolist()
KD_union_loop['prob_shNT_plusCL'] = tmp.loc[KD_union_loop['label'].tolist(), 'prob_shNT_plusCL'].tolist()
KD_union_loop['prob_delta'] = KD_union_loop['prob_KD_plusCL'] - KD_union_loop['prob_shNT_plusCL']

KD_union_loop_df = _get_values_(KD_union_loop, comp = 'plusCL_KD_vs_shNT', value = 'stat')

## CL induce and H2AZ KD impaired genes
g1 = diffexp['plusCL_KD_vs_shNT'].query('log2FoldChange < 0')
cut = np.log2(1.5)
g2 = diffexp['shNT_plusCL_vs_minusCL'].query('log2FoldChange > @cut and padj < 0.05')
g = np.intersect1d(g1.index, g2.index)

## down
index1 = [len(np.intersect1d(x.split(';'), g)) > 0 if not pd.isna(x) else False for x in KD_union_loop_df['r1_gene'].tolist()]
index2 = [len(np.intersect1d(x.split(';'), g)) > 0 if not pd.isna(x) else False for x in KD_union_loop_df['r2_gene'].tolist()]
focus_df = KD_union_loop_df.loc[np.array(index1) | np.array(index2),]

focus_df['r1_gene'] = [';'.join(list(np.intersect1d(x.split(';'), g))) if not pd.isna(x) else x for x in focus_df['r1_gene'].tolist()]
focus_df['r2_gene'] = [';'.join(list(np.intersect1d(x.split(';'), g))) if not pd.isna(x) else x for x in focus_df['r2_gene'].tolist()]
focus_df['r1_gene'] = [np.nan if x == '' else x for x in focus_df['r1_gene'].tolist()]
focus_df['r2_gene'] = [np.nan if x == '' else x for x in focus_df['r2_gene'].tolist()]

focus_df = _get_values_(focus_df.drop(['stat'], axis = 1), comp = 'plusCL_KD_vs_shNT', value = 'stat')


fig, ax = plt.subplots(figsize = (5,5))
sns.scatterplot(data = focus_df.query('group != "OO"'), x = 'stat', y = 'score_delta',
                edgecolor = 'none', s = 5, alpha = .7)
ax.set(xlabel = 'Differential mRNA expression\n(z score of H2AZ KD vs WT)',
      ylabel = 'Delta contact frequency\n(H2AZ KD - WT)')
ax.set_ylim(-0.015, 0.015)
plt.tight_layout()
# fig.savefig('Figures/CL_induce_H2AZKD_down_gene_associated_loop_changes_scatter.pdf')
plt.show()


# ### H2AZ ChIP-seq data integration

# In[4]:


## ChIP-seq path
H2AZ_path = {'ChIP_04082022_H2AZ_plusCL_rep1':'../H2AZ_Danpos/ChIP_04082022_H2AZ_plusCL_vs_IgG_rep1/pooled/ChIP_04082022_H2AZ_plusCL_unique_rep1.sorted.bam.dedup.bgsub.Fnor.peaks.xls',
            'ChIP_04082022_H2AZ_plusCL_rep2':'../H2AZ_Danpos/ChIP_04082022_H2AZ_plusCL_vs_IgG_rep2/pooled/ChIP_04082022_H2AZ_plusCL_unique_rep2.sorted.bam.dedup.bgsub.Fnor.peaks.xls',
            'ChIP_04082022_H2AZ_plusCL_rep3':'../H2AZ_Danpos/ChIP_04082022_H2AZ_plusCL_vs_IgG_rep3/pooled/ChIP_04082022_H2AZ_plusCL_unique_rep3.sorted.bam.dedup.bgsub.Fnor.peaks.xls',
            }

H2AZ_peaks = {}
for x in H2AZ_path:
    tmp = pd.read_csv(H2AZ_path[x], sep = '\t').query('width_above_cutoff > 147')
    H2AZ_peaks[x] = pybedtools.BedTool.from_dataframe(tmp)
    
### highly reproduciable peaks by replicates
H2AZ_plusCL_reprod_peaks = H2AZ_peaks['ChIP_04082022_H2AZ_plusCL_rep2'].intersect(H2AZ_peaks['ChIP_04082022_H2AZ_plusCL_rep1'], wa=True).intersect(H2AZ_peaks['ChIP_04082022_H2AZ_plusCL_rep3'], wa=True).to_dataframe().drop_duplicates()
H2AZ_plusCL_reprod_peaks_bed = pybedtools.BedTool.from_dataframe(H2AZ_plusCL_reprod_peaks)

H2AZ_bw = {'ChIP_04082022_H2AZ_plusCL_rep1':'../H2AZ_Macs/H2AZCLrep1.rep1_treat_pileup.bw',
            'ChIP_04082022_H2AZ_plusCL_rep2':'../H2AZ_Macs/H2AZCLrep2.rep1_treat_pileup.bw',
            'ChIP_04082022_H2AZ_plusCL_rep3':'../H2AZ_Macs/H2AZCLrep3.rep1_treat_pileup.bw',
}


# In[48]:


### identify H2A.Z occupied loops within H2AZ-dependant loops
EP_res = _plot_H2AZ_overlap_v2_(loops=CL_induce_KD_down_EP, peaks={'H2AZ_peak': H2AZ_plusCL_reprod_peaks_bed}, title = 'CL_induce_KD_down_EP_H2AZ_overlap',flip = True)
PP_res = _plot_H2AZ_overlap_v2_(loops=CL_induce_KD_down_PP, peaks={'H2AZ_peak': H2AZ_plusCL_reprod_peaks_bed}, title = 'CL_induce_KD_down_PP_H2AZ_overlap',flip = True)
PO_res = _plot_H2AZ_overlap_v2_(loops=CL_induce_KD_down_PO, peaks={'H2AZ_peak': H2AZ_plusCL_reprod_peaks_bed}, title = 'CL_induce_KD_down_PO_H2AZ_overlap',flip = True)
OO_res = _plot_H2AZ_overlap_v2_(loops=CL_induce_KD_down_OO, peaks={'H2AZ_peak': H2AZ_plusCL_reprod_peaks_bed}, title = 'CL_induce_KD_down_OO_H2AZ_overlap',flip = True)



# In[49]:


### stacked barplot showing H2AZ loop percentage

df = pd.DataFrame({'E-P':EP_res['H2AZ_peak'], 'P-P': PP_res['H2AZ_peak'], 'P-O':PO_res['H2AZ_peak'], 'Other': OO_res['H2AZ_peak']})
df1 = df.loc[['anchor1', 'anchor2', 'both_anchor', 'either_anchor']]*100/df.loc['total_loop']
df1 = df1.reset_index().melt(id_vars = ['index'])
df1 = df1[df1['index'].isin(['both_anchor', 'either_anchor'])]

fig, ax = plt.subplots(figsize = (5, 4))
sns.barplot(data = df1.groupby('variable').sum().reset_index(), 
            x = 'variable', y = 'value', color = 'navy', order = ['E-P', 'P-P', 'P-O', 'Other'],
           label = 'either_anchor')
sns.barplot(data = df1.query('index == "both_anchor"'), 
            x = 'variable', y = 'value', color = 'lightblue', order = ['E-P', 'P-P', 'P-O', 'Other'],
           label = 'both_anchor')
# ax.set_ylim(0, 100)
ax.tick_params(axis = 'x', rotation = 90)
# ax.tick_params(axis = 'x', rotation = 90)
ax.set(xlabel = '', ylabel = '% of loops \noverlapped with H2A.Z peaks')
ax.set_title(label = '', pad = 20)
#     ax.set_xlim(0, 45)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon = False)
plt.tight_layout()
sns.despine()
# fig.savefig('CL_induce_H2AZ_down_percent_H2AZ_stacked.pdf')
plt.show()



# In[63]:


## label loops with H2AZ binding status
both_anchor_loops = pd.DataFrame()
for t, r, l in [['EP', EP_res, CL_induce_KD_down_EP],
         ['PP', PP_res, CL_induce_KD_down_PP],
         ['PO', PO_res, CL_induce_KD_down_PO],
         ['OO', OO_res, CL_induce_KD_down_OO]]:
    both_anchor_loops = pd.concat([both_anchor_loops, 
                                   r['H2AZ_peak']['loops']['both']])
    
either_anchor_loops = pd.DataFrame()
for t, r, l in [['EP', EP_res, CL_induce_KD_down_EP],
         ['PP', PP_res, CL_induce_KD_down_PP],
         ['PO', PO_res, CL_induce_KD_down_PO],
         ['OO', OO_res, CL_induce_KD_down_OO]]:
    either_anchor_loops = pd.concat([either_anchor_loops, 
                                   r['H2AZ_peak']['loops']['either']])
    
    
H2AZ_status = []
for x in CL_induce_KD_down_tmp['label'].tolist():
    if x in both_anchor_loops['label'].tolist():
        H2AZ_status.append('both')
    elif x in either_anchor_loops['label'].tolist():
        H2AZ_status.append('either')
    else:
        H2AZ_status.append('unbound')

CL_induce_KD_down_tmp['H2AZ_status'] = H2AZ_status
CL_induce_KD_down_tmp['score_delta_KD_vs_shNT']=CL_induce_KD_down_tmp['score_KD_plusCL'] - CL_induce_KD_down_tmp['score_shNT_plusCL']
CL_induce_KD_down_tmp['score_delta_plusCL_vs_minusCL']=CL_induce_KD_down_tmp['score_shNT_plusCL'] - CL_induce_KD_down_tmp['score_shNT_minusCL']
CL_induce_KD_down_tmp['loop_length'] = CL_induce_KD_down_tmp['start2']-CL_induce_KD_down_tmp['start1']
CL_induce_KD_down_tmp['H2AZ_status_new'] = ['bound' if x != 'unbound' else x for x in CL_induce_KD_down_tmp['H2AZ_status'].tolist()]


# In[62]:


### boxplot showing the loop length of H2AZ-occupied loops

fig, ax = plt.subplots(figsize = (6, 4))
sns.boxplot(data = CL_induce_KD_down_tmp.query('group != "O-O"'), x = 'group', 
               y = 'loop_length', hue = 'H2AZ_status', showfliers = False,
            palette = {'unbound':'lightgrey', 'either': 'lightgreen', 'both':'green'}
           )
sns.despine()
ax.legend(loc = 'center left', bbox_to_anchor=(1, 0.5), fontsize = 8)
plt.tight_layout()
fig.savefig('./Figures/loop_length_H2AZ_occupied.pdf')
plt.show()


# In[68]:


### APA for H2AZ occupied loops
H2AZ_occupied_loop = CL_induce_KD_down_tmp.query('H2AZ_status_new == "bound"').drop_duplicates('label')
### plusCL
pup1 = coolpup.pileup(clr1, H2AZ_occupied_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected3,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup1,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
#              vmax = 4, vmin = 1.3,
             height=1.5)
plt.show()

plotpup.plot(pup1,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
#              vmax = 4, vmin = 1.3,
             height=1.5)
plt.show()

### minusCL
pup2 = coolpup.pileup(clr2, H2AZ_occupied_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected3,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup2,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 3.4, vmin = 1.2,
             height=1.5)
plt.show()

plotpup.plot(pup2,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 3.4, vmin = 1.2,
             height=1.5)
plt.show()


### H2AZ KD
pup3 = coolpup.pileup(clr3, H2AZ_occupied_loop,
                      features_format='bedpe', view_df=mm10_arms,
                      expected_df = expected3,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup3,
             score=False,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 3.4, vmin = 1.2,
             height=1.5)
plt.show()

plotpup.plot(pup3,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 3.4, vmin = 1.2,
             height=1.5)
plt.show()



# ### Integrate ATAC-seq data with H2AZ and loops

# In[5]:


## ATAC-seq path
ATAC_path = {'shNT_minusCL_rep2':'../ATAC/danpos/e200/shNT_minusCL_rep2/pooled/bam_shNT_minusCL_rep2_unique.sorted.bam.dedup.Fnor.peaks.xls',
            'shNT_minusCL_rep3':'../ATAC/danpos/e200/shNT_minusCL_rep3/pooled/bam_shNT_minusCL_rep3_unique.sorted.bam.dedup.Fnor.peaks.xls',
            'shNT_plusCL_rep2':'../ATAC/danpos/e200/shNT_plusCL_rep2/pooled/bam_shNT_plusCL_rep2_unique.sorted.bam.dedup.Fnor.peaks.xls',
            'shNT_plusCL_rep3':'../ATAC/danpos/e200/shNT_plusCL_rep3/pooled/bam_shNT_plusCL_rep3_unique.sorted.bam.dedup.Fnor.peaks.xls',
            'KD_plusCL_rep2':'../ATAC/danpos/e200/KD_plusCL_rep2/pooled/bam_KD_plusCL_rep2_unique.sorted.bam.dedup.Fnor.peaks.xls',
            'KD_plusCL_rep3':'../ATAC/danpos/e200/KD_plusCL_rep3/pooled/bam_KD_plusCL_rep3_unique.sorted.bam.dedup.Fnor.peaks.xls',
            }

ATAC_peaks = {}
for x in ATAC_path:
    tmp = pd.read_csv(ATAC_path[x], sep = '\t')
    ATAC_peaks[x] = pybedtools.BedTool.from_dataframe(tmp)
## reproduciable peaks in two replicates
ATAC_peaks_reprod =  {'shNT_minusCL':ATAC_peaks['shNT_minusCL_rep2'].intersect(ATAC_peaks['shNT_minusCL_rep3'], wa=True).to_dataframe().drop_duplicates(),
                      'shNT_plusCL':ATAC_peaks['shNT_plusCL_rep2'].intersect(ATAC_peaks['shNT_plusCL_rep3'], wa=True).to_dataframe().drop_duplicates(),
                      'KD_plusCL':ATAC_peaks['KD_plusCL_rep2'].intersect(ATAC_peaks['KD_plusCL_rep3'], wa=True).to_dataframe().drop_duplicates()}

ATAC_peaks_reprod_bed =  {x:pybedtools.BedTool.from_dataframe(ATAC_peaks_reprod[x]) for x in ATAC_peaks_reprod}


bw_path_list = {#'shNT_minusCL_rep1':'../ATAC/bw/shNT_minusCL_rep1.rep1_treat_pileup.bw',
            'shNT_minusCL_rep2':'../ATAC/bw/shNT_minusCL_rep2.rep1_treat_pileup.bw',
            'shNT_minusCL_rep3':'../ATAC/bw/shNT_minusCL_rep3.rep1_treat_pileup.bw',
            #'shNT_plusCL_rep1':'../ATAC/bw/shNT_plusCL_rep1.rep1_treat_pileup.bw',
            'shNT_plusCL_rep2':'../ATAC/bw/shNT_plusCL_rep2.rep1_treat_pileup.bw',
            'shNT_plusCL_rep3':'../ATAC/bw/shNT_plusCL_rep3.rep1_treat_pileup.bw',
            #'KD_plusCL_rep1':'../ATAC/bw/KD_plusCL_rep1.rep1_treat_pileup.bw',
            'KD_plusCL_rep2':'../ATAC/bw/KD_plusCL_rep2.rep1_treat_pileup.bw',
            'KD_plusCL_rep3':'../ATAC/bw/KD_plusCL_rep3.rep1_treat_pileup.bw',
            }


# In[6]:


out = open('loop_set.pk', 'rb')
d = pk.load(out)
out.close()

union_loops= d['union_loops']
union_loops_genes = d['union_loops_genes']
up_loop = d['up_loop']
down_loop = d['down_loop']
up_loop_union = d['up_loop_union']
CL_induce_KD_down = d['CL_induce_KD_down']
CL_induce_KD_down_tmp = d['CL_induce_KD_down_tmp']
KD_down_loop = d['KD_down_loop']
shNT_plusCL_allLoop = d['shNT_plusCL_allLoop']
EP_res=d['EP_res']
PP_res=d['PP_res']
PO_res=d['PO_res']
OO_res=d['OO_res']
CL_induce_KD_down_EP=d['CL_induce_KD_down_EP']
CL_induce_KD_down_PP=d['CL_induce_KD_down_PP']
CL_induce_KD_down_PO=d['CL_induce_KD_down_PO']
CL_induce_KD_down_OO=d['CL_induce_KD_down_OO']
        


# In[71]:


### venn diagram for ATAC total peaks and H2AZ total peaks
venn_mpl(a=H2AZ_plusCL_reprod_peaks_bed, b=ATAC_peaks_reprod_bed['shNT_plusCL'], colors=['green', 'purple'], 
             outfn='Figures/H2AZ_ATAC_peak_venn_shNT_plusCL.pdf', labels=['H2A.Z', 'ATAC'],text=False,
             dpi=300, figsize = (5, 4))

### venn diagram for ATAC peaks and H2AZ peaks in plucCL loops
anchors = shNT_plusCL_allLoop[['chrom1', 'start1', 'end1']].values.tolist()+shNT_plusCL_allLoop[['chrom2', 'start2', 'end2']].values.tolist()
anchors = pd.DataFrame(anchors, columns = ['chrom', 'start', 'end']).drop_duplicates()
anchor_bed = pybedtools.BedTool.from_dataframe(anchors)

H2AZ_plusCL_in_loop = H2AZ_plusCL_reprod_peaks_bed.intersect(anchor_bed, wa = True)
ATAC_plusCL_in_loop = ATAC_peaks_reprod_bed['shNT_plusCL'].intersect(anchor_bed, wa = True)

venn_mpl(a=H2AZ_plusCL_in_loop, b=ATAC_plusCL_in_loop, colors=['green', 'purple'], 
             outfn='Figures/H2AZ_ATAC_peak_venn_shNT_plusCL_loop.pdf', labels=['H2A.Z', 'ATAC'],text=False,
             dpi=300, figsize = (5, 4))


# In[7]:


### signal of ATAC-seq in H2AZ peak overlapped peaks and non-overlapped peaks - all peaks
ATAC_withH2AZ =ATAC_peaks_reprod_bed['shNT_plusCL'].intersect(H2AZ_plusCL_reprod_peaks_bed, wa = True).to_dataframe().drop_duplicates()
ATAC_noH2AZ = ATAC_peaks_reprod_bed['shNT_plusCL'].intersect(H2AZ_plusCL_reprod_peaks_bed, v=True, wa = True).to_dataframe().drop_duplicates()

bw_value2 = {'bind':[],'unbind':[]}
d = 2000
bw = pyBigWig.open(bw_path_list['shNT_plusCL_rep2'])
for i, line in ATAC_withH2AZ.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    bw_value2['bind'].append(bw.values(c, center-d, center+d))

for i, line in ATAC_noH2AZ.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    bw_value2['unbind'].append(bw.values(c, center-d, center+d))

bw_value3 = {'bind':[],'unbind':[]}
bw = pyBigWig.open(bw_path_list['shNT_plusCL_rep3'])
for i, line in ATAC_withH2AZ.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    bw_value3['bind'].append(bw.values(c, center-d, center+d))

for i, line in ATAC_noH2AZ.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    bw_value3['unbind'].append(bw.values(c, center-d, center+d))
shNT_bind = np.array(bw_value2['bind']+bw_value3['bind']).mean(axis = 0)
shNT_unbind = np.array(bw_value2['unbind']+bw_value3['unbind']).mean(axis = 0)

    
bw_value2 = {'bind':[],'unbind':[]}
d = 2000
bw = pyBigWig.open(bw_path_list['KD_plusCL_rep2'])
for i, line in ATAC_withH2AZ.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    bw_value2['bind'].append(bw.values(c, center-d, center+d))

for i, line in ATAC_noH2AZ.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    bw_value2['unbind'].append(bw.values(c, center-d, center+d))

bw_value3 = {'bind':[],'unbind':[]}
bw = pyBigWig.open(bw_path_list['KD_plusCL_rep3'])
for i, line in ATAC_withH2AZ.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    bw_value3['bind'].append(bw.values(c, center-d, center+d))

for i, line in ATAC_noH2AZ.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    bw_value3['unbind'].append(bw.values(c, center-d, center+d))
KD_bind = np.array(bw_value2['bind']+bw_value3['bind']).mean(axis = 0)
KD_unbind = np.array(bw_value2['unbind']+bw_value3['unbind']).mean(axis = 0)

fig,ax = plt.subplots(figsize = (5.5, 4))
ax.plot(range(1,4001), shNT_bind, label = 'H2A.Z bind', color = 'purple')
ax.plot(range(1,4001), shNT_unbind, label = 'H2A.Z unbind', color = 'purple', linestyle = 'dashed')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set(xlabel = 'position from peak start', ylabel = 'ATAC-seq signal',
      title = 'All ATAC peaks')
plt.tight_layout()
sns.despine()
fig.savefig('Figures/ATAC_signal_H2AZ_bind_vs_unbind_plusCL.pdf')
plt.show()
plt.close()


### delta ATAC-seq signal
df_KD_shNT_bind = KD_bind - shNT_bind
df_KD_shNT_unbind = KD_unbind - shNT_unbind

fig,ax = plt.subplots(figsize = (5.5, 4))
ax.plot(range(1,4001), df_KD_shNT_bind, label = 'bind', color = 'purple')
ax.plot(range(1,4001), df_KD_shNT_unbind, label = 'unbind', color = 'purple', linestyle = 'dashed')
ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
# ax.set_ylim(0, 5.5)
ax.set(xlabel = 'position from peak start', ylabel = 'ATAC-seq signal (KD - shNT)',
      title = 'All ATAC peaks')
plt.tight_layout()
fig.savefig('Figures/ATAC_all_peak_KD_delta_signal_curve_H2AZ_bind_vs_unbind.pdf')
plt.show()
plt.close()


# In[19]:


## H2AZ signal of ATAC-seq peak aligned - all peaks
ATAC_withH2AZ =ATAC_peaks_reprod_bed['shNT_plusCL'].intersect(H2AZ_plusCL_reprod_peaks_bed, wa = True).to_dataframe().drop_duplicates()
ATAC_noH2AZ = ATAC_peaks_reprod_bed['shNT_plusCL'].intersect(H2AZ_plusCL_reprod_peaks_bed, v=True, wa = True).to_dataframe().drop_duplicates()

H2AZ_bw_value1 = {'bind':[]}
d = 2000
bw = pyBigWig.open(H2AZ_bw['ChIP_04082022_H2AZ_plusCL_rep1'])
for i, line in ATAC_withH2AZ.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    H2AZ_bw_value1['bind'].append(bw.values(c, center-d, center+d))

H2AZ_bw_value2 = {'bind':[]}
bw = pyBigWig.open(H2AZ_bw['ChIP_04082022_H2AZ_plusCL_rep2'])
for i, line in ATAC_withH2AZ.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    H2AZ_bw_value2['bind'].append(bw.values(c, center-d, center+d))

H2AZ_bw_value3 = {'bind':[]}
bw = pyBigWig.open(H2AZ_bw['ChIP_04082022_H2AZ_plusCL_rep3'])
for i, line in ATAC_withH2AZ.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    H2AZ_bw_value3['bind'].append(bw.values(c, center-d, center+d))

hbv_all = np.array(H2AZ_bw_value1['bind']+H2AZ_bw_value2['bind']+H2AZ_bw_value3['bind']).mean(axis = 0)


fig,ax = plt.subplots(figsize = (5.5, 4))
ax.plot(range(1,4001), hbv_all, label = 'H2A.Z bind', color = 'green')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set(xlabel = 'position from peak start\n(ATAC-seq peak)', ylabel = 'H2AZ ChIP-seq signal',
      title = 'All ATAC peaks')
plt.tight_layout()
sns.despine()
fig.savefig('Figures/H2AZ_signal_H2AZ_bind_ATAC_peaks_plusCL.pdf')
plt.show()
plt.close()


# In[10]:


### signal of ATAC-seq in H2AZ peak overlapped peaks and non-overlapped peaks -- within plusCL loop anchors

anchors = shNT_plusCL_allLoop[['chrom1', 'start1', 'end1']].values.tolist()+shNT_plusCL_allLoop[['chrom2', 'start2', 'end2']].values.tolist()
anchors = pd.DataFrame(anchors, columns = ['chrom', 'start', 'end']).drop_duplicates()
anchor_bed = pybedtools.BedTool.from_dataframe(anchors)

H2AZ_plusCL_in_loop = H2AZ_plusCL_reprod_peaks_bed.intersect(anchor_bed, wa = True)
ATAC_plusCL_in_loop = ATAC_peaks_reprod_bed['shNT_plusCL'].intersect(anchor_bed, wa = True)

ATAC_withH2AZ_inloop = ATAC_plusCL_in_loop.intersect(H2AZ_plusCL_reprod_peaks_bed, wa = True).to_dataframe().drop_duplicates()
ATAC_noH2AZ_inloop = ATAC_plusCL_in_loop.intersect(H2AZ_plusCL_reprod_peaks_bed, v=True, wa = True).to_dataframe().drop_duplicates()

bw_value2 = {'bind':[],'unbind':[]}
d = 2000
bw = pyBigWig.open(bw_path_list['shNT_plusCL_rep2'])
for i, line in ATAC_withH2AZ_inloop.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    bw_value2['bind'].append(bw.values(c, center-d, center+d))

for i, line in ATAC_noH2AZ_inloop.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    bw_value2['unbind'].append(bw.values(c, center-d, center+d))

bw_value3 = {'bind':[],'unbind':[]}
bw = pyBigWig.open(bw_path_list['shNT_plusCL_rep3'])
for i, line in ATAC_withH2AZ_inloop.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    bw_value3['bind'].append(bw.values(c, center-d, center+d))

for i, line in ATAC_noH2AZ_inloop.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    bw_value3['unbind'].append(bw.values(c, center-d, center+d))
shNT_bind_inloop = np.array(bw_value2['bind']+bw_value3['bind']).mean(axis = 0)
shNT_unbind_inloop = np.array(bw_value2['unbind']+bw_value3['unbind']).mean(axis = 0)

    
bw_value2 = {'bind':[],'unbind':[]}
d = 2000
bw = pyBigWig.open(bw_path_list['KD_plusCL_rep2'])
for i, line in ATAC_withH2AZ_inloop.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    bw_value2['bind'].append(bw.values(c, center-d, center+d))

for i, line in ATAC_noH2AZ_inloop.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    bw_value2['unbind'].append(bw.values(c, center-d, center+d))

bw_value3 = {'bind':[],'unbind':[]}
bw = pyBigWig.open(bw_path_list['KD_plusCL_rep3'])
for i, line in ATAC_withH2AZ_inloop.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    bw_value3['bind'].append(bw.values(c, center-d, center+d))

for i, line in ATAC_noH2AZ_inloop.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    bw_value3['unbind'].append(bw.values(c, center-d, center+d))
KD_bind_inloop = np.array(bw_value2['bind']+bw_value3['bind']).mean(axis = 0)
KD_unbind_inloop = np.array(bw_value2['unbind']+bw_value3['unbind']).mean(axis = 0)

fig,ax = plt.subplots(figsize = (5.5, 4))
ax.plot(range(1,4001), shNT_bind_inloop, label = 'H2A.Z bind', color = 'purple')
ax.plot(range(1,4001), shNT_unbind_inloop, label = 'H2A.Z unbind', color = 'purple', linestyle = 'dashed')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set(xlabel = 'position from peak start', ylabel = 'ATAC-seq signal',
      title = 'All ATAC peaks in shNT plusCL loop')
plt.tight_layout()
sns.despine()
fig.savefig('Figures/ATAC_signal_H2AZ_bind_vs_unbind_plusCL_in_shNT_plusCL_loops.v2.pdf')
plt.show()
plt.close()


### delta ATAC-seq signal
df_KD_shNT_bind_inloop = KD_bind_inloop - shNT_bind_inloop
df_KD_shNT_unbind_inloop = KD_unbind_inloop - shNT_unbind_inloop

fig,ax = plt.subplots(figsize = (5.5, 4))
ax.plot(range(1,4001), df_KD_shNT_bind_inloop, label = 'bind', color = 'purple')
ax.plot(range(1,4001), df_KD_shNT_unbind_inloop, label = 'unbind', color = 'purple', linestyle = 'dashed')
ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
# ax.set_ylim(0, 5.5)
ax.set(xlabel = 'position from peak start', ylabel = 'ATAC-seq signal (KD - shNT)',
      title = 'All ATAC peaks in plusCL loop')
plt.tight_layout()
fig.savefig('Figures/ATAC_in_plusCL_loop_KD_delta_signal_curve_H2AZ_bind_vs_unbind.pdf')
plt.show()
plt.close()


# In[18]:


## H2AZ signal of ATAC-seq peak aligned -- within plusCL loop anchors

H2AZ_bw_value1 = {'bind':[]}
d = 2000
bw = pyBigWig.open(H2AZ_bw['ChIP_04082022_H2AZ_plusCL_rep1'])
for i, line in ATAC_withH2AZ.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    H2AZ_bw_value1['bind'].append(bw.values(c, center-d, center+d))

H2AZ_bw_value2 = {'bind':[]}
bw = pyBigWig.open(H2AZ_bw['ChIP_04082022_H2AZ_plusCL_rep2'])
for i, line in ATAC_withH2AZ.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    H2AZ_bw_value2['bind'].append(bw.values(c, center-d, center+d))

H2AZ_bw_value3 = {'bind':[]}
bw = pyBigWig.open(H2AZ_bw['ChIP_04082022_H2AZ_plusCL_rep3'])
for i, line in ATAC_withH2AZ.iterrows():
    c, s, e = line['chrom'], line['start'], line['end']
    center = s+int((np.abs(s-e)/2))
    H2AZ_bw_value3['bind'].append(bw.values(c, center-d, center+d))

hbv = np.array(H2AZ_bw_value1['bind']+H2AZ_bw_value2['bind']+H2AZ_bw_value3['bind']).mean(axis = 0)


fig,ax = plt.subplots(figsize = (5.5, 4))
ax.plot(range(1,4001), hbv, label = 'H2A.Z bind', color = 'green')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set(xlabel = 'position from peak start\n(ATAC-seq peak)', ylabel = 'H2AZ ChIP-seq signal',
      title = 'All ATAC peaks in shNT plusCL loop')
plt.tight_layout()
sns.despine()
fig.savefig('Figures/H2AZ_signal_H2AZ_bind_ATAC_peaks_plusCL_in_shNT_plusCL_loops.pdf')
plt.show()
plt.close()


# In[13]:


### ATAC-seq signal at H2AZ site in H2AZ occupied loops in three conditions
bw_H2AZ_site_values={}
for t, r, l in [['EP', EP_res, CL_induce_KD_down_EP],
         ['PP', PP_res, CL_induce_KD_down_PP],
         ['PO', PO_res, CL_induce_KD_down_PO],
         ['OO', OO_res, CL_induce_KD_down_OO]]:
    H2AZ_bind_loops = pd.concat([r['H2AZ_peak']['loops']['either'], r['H2AZ_peak']['loops']['both']])
    H2AZ_bind_loops['anchor1'] = H2AZ_bind_loops['chrom1']+':'+H2AZ_bind_loops['start1'].astype('str')+'-'+H2AZ_bind_loops['end1'].astype('str')
    H2AZ_bind_loops['anchor2'] = H2AZ_bind_loops['chrom2']+':'+H2AZ_bind_loops['start2'].astype('str')+'-'+H2AZ_bind_loops['end2'].astype('str')
   
    bind_anchors = pd.DataFrame(H2AZ_bind_loops[['chrom1', 'start1', 'end1']].values.tolist()+H2AZ_bind_loops[['chrom2', 'start2', 'end2']].values.tolist())
    bind_anchors_bed = pybedtools.BedTool.from_dataframe(bind_anchors)
    H2AZ_site = H2AZ_plusCL_reprod_peaks_bed.intersect(bind_anchors_bed, wa = True).to_dataframe().drop_duplicates()
    H2AZ_site['start'] = H2AZ_site['start'] + ((H2AZ_site['start']-H2AZ_site['end'])/2).abs().astype('int')
    
    d = 2000
    bw_H2AZ_site_values[t] = []
    for path in bw_path_list:
        bw = pyBigWig.open(bw_path_list[path])
        for i, line in H2AZ_site.iterrows():
            line['bw'] = bw.values(line['chrom'], line['start']-d, line['start']+d)
            bw_H2AZ_site_values[t].append(line.values.tolist()+[path])


df = pd.concat([pd.DataFrame(bw_H2AZ_site_values['EP']),
               pd.DataFrame(bw_H2AZ_site_values['PP']),
               pd.DataFrame(bw_H2AZ_site_values['PO']),
               pd.DataFrame(bw_H2AZ_site_values['OO'])]).groupby(9).apply(lambda d: pd.DataFrame(d[8].tolist()).mean())

df['cond'] = [x.replace('_rep2', '').replace('_rep3', '') for x in df.index.tolist()]
df = df.groupby('cond').mean()

df22 = df.loc[:,range(1000, 3000)]
fig,ax = plt.subplots(figsize = (6, 4))
ax.plot(range(0, df22.shape[1]), df22.loc['KD_plusCL',:], label = 'KD_plusCL', linewidth = 1.5)
ax.plot(range(0, df22.shape[1]), df22.loc['shNT_plusCL',:], label = 'shNT_plusCL', linewidth = 1.5)
ax.plot(range(0, df22.shape[1]), df22.loc['shNT_minusCL',:], label = 'shNT_minusCL', linewidth = 1.5, color = 'grey')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon = False)
ax.set(xlabel = 'H2AZ peak center', ylabel = 'ATAC-seq signal')
sns.despine()
ax.set_ylim(0.6, 2.3)
plt.tight_layout()
fig.savefig('Figures/ATAC_signal_at_H2AZ_site_in_H2AZ_dependant_loops_three_condition.pdf')
plt.show()


# ## output supplementary tables
# 

# In[ ]:


# Supplementary Table 1. Chromatin loops in vehicle-treated and b3-AR stimulated brown adipocytes
with pd.ExcelWriter('SuppTable/SuppTable1.xlsx') as writer:
    shNT_minusCL_allLoop.iloc[:,0:8].drop_duplicates().to_excel(writer, sheet_name='vehicle-treated', index = None)
    shNT_plusCL_allLoop.iloc[:,0:8].drop_duplicates().to_excel(writer, sheet_name='b3-AR stimulated', index = None)

# Supplementary Table 2. Differential loops between vehicle-treated and b3-AR stimulated brown adipocytes
with pd.ExcelWriter('SuppTable/SuppTable2.xlsx') as writer:
    up_loop_union.loc[up_loop['label'],:].iloc[:,:10].drop_duplicates().to_excel(writer, sheet_name='UP', index = None)
    down_loop_union.loc[down_loop['label'],:].iloc[:,:10].drop_duplicates().to_excel(writer, sheet_name='DOWN', index = None)

# Supplementary Table 3. E-P, P-P, and P-other loops in vehicle-treated and b3-AR stimulated brown adipocytes
with pd.ExcelWriter('SuppTable/SuppTable3.xlsx') as writer:
    shNT_minusCL_EP.iloc[:,0:8].drop_duplicates().to_excel(writer, sheet_name='vehicle-treated_EP', index = None)
    shNT_minusCL_PP.iloc[:,0:8].drop_duplicates().to_excel(writer, sheet_name='vehicle-treated_PP', index = None)
    shNT_minusCL_PO.iloc[:,0:8].drop_duplicates().to_excel(writer, sheet_name='vehicle-treated_PO', index = None)
    shNT_plusCL_EP.iloc[:,0:8].drop_duplicates().to_excel(writer, sheet_name='b3-AR stimulated_EP', index = None)
    shNT_plusCL_PP.iloc[:,0:8].drop_duplicates().to_excel(writer, sheet_name='b3-AR stimulated_PP', index = None)
    shNT_plusCL_PO.iloc[:,0:8].drop_duplicates().to_excel(writer, sheet_name='b3-AR stimulated_PO', index = None)

# Supplementary Table 4. Differential E-P, P-P, and P-other loops between vehicle-treated and b3-AR stimulated brown adipocytes
with pd.ExcelWriter('SuppTable/SuppTable4.xlsx') as writer:
    up_loop_union.query('group == "E-P"').iloc[:,:10].drop_duplicates().to_excel(writer, sheet_name='UP_EP', index = None)
    up_loop_union.query('group == "P-P"').iloc[:,:10].drop_duplicates().to_excel(writer, sheet_name='UP_PP', index = None)
    up_loop_union.query('group == "P-O"').iloc[:,:10].drop_duplicates().to_excel(writer, sheet_name='UP_PO', index = None)
    down_loop_union.query('group == "E-P"').iloc[:,:10].drop_duplicates().to_excel(writer, sheet_name='DOWN_EP', index = None)
    down_loop_union.query('group == "P-P"').iloc[:,:10].drop_duplicates().to_excel(writer, sheet_name='DOWN_PP', index = None)
    down_loop_union.query('group == "P-O"').iloc[:,:10].drop_duplicates().to_excel(writer, sheet_name='DOWN_PO', index = None)

# Supplementary Table 5. Chromatin loops in b3-AR stimulated H2A.Z KD brown adipocytes
with pd.ExcelWriter('SuppTable/SuppTable5.xlsx') as writer:
    KD_plusCL_allLoop.iloc[:,0:8].drop_duplicates().to_excel(writer, sheet_name='b3-AR stimulated H2A.Z KD', index = None)

# Supplementary Table 6. Differential loops by H2A.Z KD
with pd.ExcelWriter('SuppTable/SuppTable6.xlsx') as writer:
    tmp = loop_dots.loc[KD_down_loop['label'],:].iloc[:,0:12].drop(['prob_shNT_minusCL','score_shNT_minusCL'], axis = 1)
    tmp.drop_duplicates().to_excel(writer, sheet_name='DOWN', index = None)
    tmp = loop_dots.loc[KD_up_loop['label'],:].iloc[:,0:12].drop(['prob_shNT_minusCL','score_shNT_minusCL'], axis = 1)
    tmp.drop_duplicates().to_excel(writer, sheet_name='UP', index = None)

# Supplementary Table 7. Dynamic loops affected by H2A.Z KD
with pd.ExcelWriter('SuppTable/SuppTable7.xlsx') as writer:
    CL_induce_KD_down.iloc[:,0:12].drop_duplicates().to_excel(writer, sheet_name='b3-AR_UP_H2AZ_KD_DOWN', index = None)



# ### Micro-C integrates with PRO-seq for nasent transcription

# In[ ]:


## load PRO-seq differential analysis data and save into a dict
plus_KD_vs_NT_de = pd.read_csv('/temp_work/ch228298/H2AZ_GEO_2025/GEO/geo_submission_PROseq/b3AR_stim_KD_vs_Control_nasent_transcript_PROseq_DEG.txt.gz', compression = 'gzip',
                              sep = '\t')
tmp = plus_KD_vs_NT_de[['symbol', 'log2FoldChange_KDplus.vs.NTplus', 'padj_KDplus.vs.NTplus', 'stat_KDplus.vs.NTplus']]
tmp.columns = ['gene', 'log2FoldChange', 'padj', 'stat']
tmp['padj'] = tmp['padj'].fillna(1)
tmp['log2FoldChange'] = tmp['log2FoldChange'].fillna(0)
tmp['log2FoldChange_abs'] = tmp['log2FoldChange'].abs()
## deduplicate by taking largest abs fold change
tmp = tmp.sort_values('log2FoldChange_abs', ascending = False)
tmp = tmp[~tmp['gene'].duplicated()]
tmp.index = tmp['gene'].tolist()
pro_diffexp['plusCL_KD_vs_shNT'] = tmp.copy()
del tmp

plus_vs_minus_de = pd.read_csv('/temp_work/ch228298/H2AZ_GEO_2025/GEO/geo_submission_PROseq/Control_b3AR_stim_vs_vehicle_nasent_transcript_PROseq_DEG.txt.gz', compression = 'gzip',
                               sep = '\t')

tmp = plus_vs_minus_de[['symbol', 'log2FoldChange_NTplus.vs.NTminus', 'padj_NTplus.vs.NTminus', 'stat_NTplus.vs.NTminus']]
tmp.columns = ['gene', 'log2FoldChange', 'padj', 'stat']
tmp['padj'] = tmp['padj'].fillna(1)
tmp['log2FoldChange'] = tmp['log2FoldChange'].fillna(0)
tmp['log2FoldChange_abs'] = tmp['log2FoldChange'].abs()
## deduplicate by taking largest abs fold change
tmp = tmp.sort_values('log2FoldChange_abs', ascending = False)
tmp = tmp[~tmp['gene'].duplicated()]
tmp.index = tmp['gene'].tolist()
pro_diffexp['shNT_plusCL_vs_minusCL'] = tmp.copy()
del tmp


# In[ ]:


### function to get values for loop change and expression change, lowess smoothing applied
def _get_values_pro_(df, comp = 'shNT_plusCL_vs_minusCL', value = 'log2FoldChange'):
    dx, dy = [], []
    ddf = []
    for i, line in df.iterrows():
        g = line[['r1_gene', 'r2_gene']].dropna().values.tolist()
        for j in g:
            j = j.split(';')
            for jj in j:
                ytmp = pro_diffexp[comp].loc[jj, value] if jj in pro_diffexp[comp].index.tolist() else np.nan
                dy.append(ytmp)
                dx.append(line['score_delta'])
                ddf.append(line.tolist()+[ytmp, jj])

    ddf = pd.DataFrame(ddf, columns = df.columns.tolist()+[value, 'gene'])
    # return(ddf)
    ddf = ddf[~pd.isna(ddf[value]) & ~pd.isna(ddf['score_delta'])]
    ddf = ddf.sort_values('score_delta')
    x, y = ddf['score_delta'], ddf[value]
    l = loess(x,y)
    l.fit()
    pred = l.predict(x, stderror=True)
    conf = pred.confidence()
    lowess = pred.values
    ll = conf.lower
    ul = conf.upper
    ddf['lowess'] = lowess
    ddf['lower_lowess'] = ll
    ddf['upper_lowess'] = ul
    return(ddf)
## match loop with PRO-seq detected genes, and their differential values
EP_ddf_pro = _get_values_pro_(EP_df)
PP_ddf_pro = _get_values_pro_(PP_df)
PO_ddf_pro = _get_values_pro_(PO_df)

### plot loop delta contact vs. gene expression change
fig, ax = plt.subplots()

## align different loop type range to the same, so that can be better compared
xmax = min([max(EP_ddf_pro['score_delta']), max(PP_ddf_pro['score_delta']), max(PO_ddf_pro['score_delta'])])
xmin = max([min(EP_ddf_pro['score_delta']), min(PP_ddf_pro['score_delta']), min(PO_ddf_pro['score_delta'])])

q = 'score_delta >= @xmin and score_delta <= @xmax'
ax.plot(EP_ddf_pro.drop_duplicates(['label', 'gene']).query(q)['score_delta'], EP_ddf_pro.drop_duplicates(['label', 'gene']).query(q)['lowess'], label = 'E-P')
ax.plot(PP_ddf_pro.drop_duplicates(['label', 'gene']).query(q)['score_delta'], PP_ddf_pro.drop_duplicates(['label', 'gene']).query(q)['lowess'], label = 'P-P')
ax.plot(PO_ddf_pro.drop_duplicates(['label', 'gene']).query(q)['score_delta'], PO_ddf_pro.drop_duplicates(['label', 'gene']).query(q)['lowess'], label = 'P-Other')
ax.fill_between(EP_ddf_pro.drop_duplicates(['label', 'gene']).query(q)['score_delta'],EP_ddf_pro.drop_duplicates(['label', 'gene']).query(q)['lower_lowess'],EP_ddf_pro.drop_duplicates(['label', 'gene']).query(q)['upper_lowess'],alpha=.33)
ax.fill_between(PP_ddf_pro.drop_duplicates(['label', 'gene']).query(q)['score_delta'],PP_ddf_pro.drop_duplicates(['label', 'gene']).query(q)['lower_lowess'],PP_ddf_pro.drop_duplicates(['label', 'gene']).query(q)['upper_lowess'],alpha=.33)
ax.fill_between(PO_ddf_pro.drop_duplicates(['label', 'gene']).query(q)['score_delta'],PO_ddf_pro.drop_duplicates(['label', 'gene']).query(q)['lower_lowess'],PO_ddf_pro.drop_duplicates(['label', 'gene']).query(q)['upper_lowess'],alpha=.33)

ax.set(xlabel = 'Delta Contact Frequency', ylabel = 'mRNA log2FoldChange',
      title = 'd=%s'%d)
ax.tick_params(axis = 'x', rotation = 45)
ax.legend(title = '', loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()
sns.despine()
plt.show()
plt.close()



# In[ ]:


### analyze the loop change and differental expression upon H2AZ KD, particularly for H2AZ-dependant loops and genes
KD_union_loop_df = _get_values_pro_(KD_union_loop, comp = 'plusCL_KD_vs_shNT', value = 'stat')

## CL induce and H2AZ KD impaired genes
g1 = pro_diffexp['plusCL_KD_vs_shNT'].query('log2FoldChange < 0')
cut = np.log2(1.5)
g2 = pro_diffexp['shNT_plusCL_vs_minusCL'].query('log2FoldChange > @cut and padj < 0.05')
g = np.intersect1d(g1.index, g2.index)

index1 = [len(np.intersect1d(x.split(';'), g)) > 0 if not pd.isna(x) else False for x in KD_union_loop_df['r1_gene'].tolist()]
index2 = [len(np.intersect1d(x.split(';'), g)) > 0 if not pd.isna(x) else False for x in KD_union_loop_df['r2_gene'].tolist()]
focus_df = KD_union_loop_df.loc[np.array(index1) | np.array(index2),]

focus_df['r1_gene'] = [';'.join(list(np.intersect1d(x.split(';'), g))) if not pd.isna(x) else x for x in focus_df['r1_gene'].tolist()]
focus_df['r2_gene'] = [';'.join(list(np.intersect1d(x.split(';'), g))) if not pd.isna(x) else x for x in focus_df['r2_gene'].tolist()]
focus_df['r1_gene'] = [np.nan if x == '' else x for x in focus_df['r1_gene'].tolist()]
focus_df['r2_gene'] = [np.nan if x == '' else x for x in focus_df['r2_gene'].tolist()]

focus_ddf2 = _get_values_pro_(focus_df.drop(['stat'], axis = 1), comp = 'plusCL_KD_vs_shNT', value = 'stat')

## plot
fig, ax = plt.subplots(figsize = (5,5))
sns.scatterplot(data = focus_ddf2.query('group != "OO"'), x = 'stat', y = 'score_delta',
                edgecolor = 'none', s = 5, alpha = .7)
ax.set(xlabel = 'Differential mRNA expression\n(z score of H2AZ KD vs WT)',
      ylabel = 'Delta contact frequency\n(H2AZ KD - WT)')
ax.set_ylim(-0.015, 0.015)
plt.tight_layout()
plt.show()
plt.close()

print('# of delta > 0 and < 0 dots:', focus_ddf2.query('score_delta > 0').shape[0], focus_ddf2.query('score_delta < 0').shape[0])



# In[ ]:


## compare differential PRO-seq expression for different loop types upon H2AZ KD

tmp = pd.DataFrame(quantile_norm(loop_dots[['score_shNT_plusCL', 'score_shNT_minusCL', 'score_KD_plusCL']].copy()), 
                   columns = ['score_shNT_plusCL', 'score_shNT_minusCL', 'score_KD_plusCL'], index = loop_dots['label'].tolist())
loop_dots['log2fc'] = loop_dots.apply(lambda row: np.log2(row['prob_shNT_plusCL'])-np.log2(np.clip(row['prob_shNT_minusCL'], a_min=0.0001, a_max = 1)), axis = 1)
up_dots = loop_dots.query('log2fc > 1')
up_dots['Type'] = 'plusCL specific'
down_dots = loop_dots.query('log2fc < -1')
down_dots['Type'] = 'minusCL specific'
stable_dots = loop_dots.query('log2fc >= -1 and log2fc <= 1')
stable_dots['Type'] = 'No'
CL_induce_KD_down = up_dots[up_dots['label'].isin(KD_down_dots['label']) & up_dots['label'].isin(up_loop['label'])]

CL_induce_KD_down['score_KD_plusCL'] = tmp.loc[CL_induce_KD_down['label'].tolist(), 'score_KD_plusCL'].tolist()
CL_induce_KD_down['score_shNT_plusCL'] = tmp.loc[CL_induce_KD_down['label'].tolist(), 'score_shNT_plusCL'].tolist()
CL_induce_KD_down['score_shNT_minusCL'] = tmp.loc[CL_induce_KD_down['label'].tolist(), 'score_shNT_minusCL'].tolist()


CL_induce_KD_down_EP, CL_induce_KD_down_PP, CL_induce_KD_down_PO, CL_induce_KD_down_OO = _define_EP_PP_(CL_induce_KD_down, comp = 'plusCL_KD_vs_shNT')
CL_induce_KD_down_EP['group'] = 'EP'
CL_induce_KD_down_PP['group'] = 'PP'
CL_induce_KD_down_PO['group'] = 'PO'
CL_induce_KD_down_OO['group'] = 'OO'
CL_induce_KD_down_union = pd.concat([CL_induce_KD_down_EP, CL_induce_KD_down_PP, 
                                       CL_induce_KD_down_PO, CL_induce_KD_down_OO])

CL_induce_KD_down_tmp = pd.concat([CL_induce_KD_down_EP, CL_induce_KD_down_PP, CL_induce_KD_down_PO, CL_induce_KD_down_OO])
CL_induce_KD_down_tmp['group'] = ['E-P']*CL_induce_KD_down_EP.shape[0]+['P-P']*CL_induce_KD_down_PP.shape[0]+['P-O']*CL_induce_KD_down_PO.shape[0]+['O-O']*CL_induce_KD_down_OO.shape[0]

CL_induce_KD_down_tmp['score_delta'] = CL_induce_KD_down_tmp['score_KD_plusCL'] - CL_induce_KD_down_tmp['score_shNT_plusCL']

CL_induce_KD_down_exp = _get_values_pro_(CL_induce_KD_down_tmp, comp = 'plusCL_KD_vs_shNT')
CL_induce_KD_down_z = _get_values_pro_(CL_induce_KD_down_tmp, comp = 'plusCL_KD_vs_shNT', value = 'stat')

df = []
for i, line in CL_induce_KD_down_z.iterrows():
    if not pd.isna(line['r1_gene']):
        for g in line['r1_gene'].split(';'):
            df.append([g, line['group'], line['r1_gene'], line['r2_gene'], line['label']])
    if not pd.isna(line['r2_gene']):
        for g in line['r2_gene'].split(';'):
            df.append([g, line['group'], line['r1_gene'], line['r2_gene'], line['label']])
df = pd.DataFrame(df, columns = ['gene', 'group', 'r1_gene', 'r2_gene', 'loop']).drop_duplicates()
df = df[df['gene'].isin(pro_diffexp['plusCL_KD_vs_shNT'].index)]
df['stat'] = [pro_diffexp['plusCL_KD_vs_shNT'].loc[x, 'stat'] for x in df['gene'].tolist()]

l1 = list(set(df.query('group == "E-P"')['gene'].tolist()))
l2 = list(set(df.query('group == "P-P"')['gene'].tolist()))
l3 = list(set(df.query('group == "P-O"')['gene'].tolist()))

l = list(np.intersect1d(l1, l2))+list(np.intersect1d(l1, l3))+list(np.intersect1d(l2, l3))

fig, ax = plt.subplots(figsize = (4,4))
sns.boxplot(data = df[~df['r1_gene'].isin(l) & ~df['r2_gene'].isin(l)], x = 'group', y = 'stat', 
            order = ['E-P', 'P-P', 'P-O'], whis = 1.5, showfliers = False, palette = 'tab10')
ax.hlines([0], *ax.get_xlim(), linestyle = 'dashed', color = 'black', linewidth = .5)
ax.set(xlabel = '', ylabel = 'Differential nascent expression\n(z score of H2AZ KD vs WT)')
plt.tight_layout()
plt.show()
plt.close()



# ### loop connectivity analysis

# In[ ]:


## perform loop connectivity analysis on plusCL loops
## the connected loops were defined as if any anchors of two loops were the same or next to each other
shNT_plusCL_allLoop['anchor1'] = ['-'.join(x.split('-')[:3]) for x in shNT_plusCL_allLoop['label'].tolist()]
shNT_plusCL_allLoop['anchor2'] = ['-'.join(x.split('-')[3:]) for x in shNT_plusCL_allLoop['label'].tolist()]
shNT_plusCL_allLoop.index = shNT_plusCL_allLoop['label'].tolist()
shNT_plusCL_loop_connect_res = []
for i1, x1 in shNT_plusCL_allLoop.iterrows():
    ## exactly match
    tmp = shNT_plusCL_allLoop[((shNT_plusCL_allLoop['anchor1'] == x1['anchor1']) | (shNT_plusCL_allLoop['anchor2'] == x1['anchor1']) | 
                 (shNT_plusCL_allLoop['anchor1'] == x1['anchor2']) | (shNT_plusCL_allLoop['anchor2'] == x1['anchor2'])) &
                       (shNT_plusCL_allLoop.index != i1)]
    for i2 in tmp.index.tolist():
        shNT_plusCL_loop_connect_res.append([i1, i2, 1])
    ## allow the loop anchor just next to each other
    x1_anchor1_chr, x1_anchor1_start, x1_anchor1_end = x1['anchor1'].split('-')
    x1_anchor1_left = '-'.join([x1_anchor1_chr, str(int(x1_anchor1_start) - 5000), str(int(x1_anchor1_end) - 5000)])
    x1_anchor2_chr, x1_anchor2_start, x1_anchor2_end = x1['anchor2'].split('-')
    x1_anchor2_left = '-'.join([x1_anchor2_chr, str(int(x1_anchor2_start) - 5000), str(int(x1_anchor2_end) - 5000)])
    tmp = shNT_plusCL_allLoop[((shNT_plusCL_allLoop['anchor1'] == x1_anchor1_left) | (shNT_plusCL_allLoop['anchor2'] == x1_anchor1_left) | 
                 (shNT_plusCL_allLoop['anchor1'] == x1_anchor2_left) | (shNT_plusCL_allLoop['anchor2'] == x1_anchor2_left)) &
                       (shNT_plusCL_allLoop.index != i1)]
    for i2 in tmp.index.tolist():
        shNT_plusCL_loop_connect_res.append([i1, i2, 1])
    ## right shit one
    x1_anchor1_right = '-'.join([x1_anchor1_chr, str(int(x1_anchor1_start) + 5000), str(int(x1_anchor1_end) + 5000)])
    x1_anchor2_right = '-'.join([x1_anchor2_chr, str(int(x1_anchor2_start) + 5000), str(int(x1_anchor2_end) + 5000)])
    tmp = shNT_plusCL_allLoop[((shNT_plusCL_allLoop['anchor1'] == x1_anchor1_right) | (shNT_plusCL_allLoop['anchor2'] == x1_anchor1_right) | 
                 (shNT_plusCL_allLoop['anchor1'] == x1_anchor2_right) | (shNT_plusCL_allLoop['anchor2'] == x1_anchor2_right)) &
                       (shNT_plusCL_allLoop.index != i1)]
    for i2 in tmp.index.tolist():
        shNT_plusCL_loop_connect_res.append([i1, i2, 1])
## same the connectivity into a dataframe format
shNT_plusCL_loop_connect_res = pd.DataFrame(shNT_plusCL_loop_connect_res)
shNT_plusCL_loop_connect_res.columns = ['loop1', 'loop2', 'connect']
shNT_plusCL_loop_connect_res['chr'] = [x.split('-')[0] for x in shNT_plusCL_loop_connect_res['loop1'].tolist()]

## identify CL up and H2AZ KD down loops using shNT plusCL anchors
shNT_plusCL_allLoop_df = loop_dots[loop_dots.index.isin(shNT_plusCL_allLoop['label'])]
## up in CL
shNT_plusCL_allLoop_df['log2fc'] = shNT_plusCL_allLoop_df.apply(lambda row: np.log2(row['prob_shNT_plusCL'])-np.log2(np.clip(row['prob_shNT_minusCL'], a_min=0.0001, a_max = 1)),
                                                                axis = 1)
shNT_plusCL_allLoop_up_df = shNT_plusCL_allLoop_df.query('log2fc > 2')
## down in KD
shNT_plusCL_allLoop_up_df['log2fc'] = shNT_plusCL_allLoop_up_df.apply(lambda row: np.log2(row['prob_KD_plusCL'])-np.log2(np.clip(row['prob_shNT_plusCL'], a_min=0.0001, a_max = 1)), axis = 1)
shNT_plusCL_allLoop_up_down_df = shNT_plusCL_allLoop_up_df.query('log2fc < -2')

## further label if the loop anchor occupied by H2AZ
anchor1 = shNT_plusCL_allLoop_up_down_df.iloc[:,0:3]
anchor1['name'] = anchor1.apply(lambda row: '-'.join(row.astype('str')), axis = 1)

anchor2 = shNT_plusCL_allLoop_up_down_df.iloc[:,3:6]
anchor2['name'] = anchor2.apply(lambda row: '-'.join(row.astype('str')), axis = 1)

shNT_plusCL_allLoop_up_down_df['anchor1'] = anchor1['name'].tolist()
shNT_plusCL_allLoop_up_down_df['anchor2'] = anchor2['name'].tolist()

anchor1_h2az = pybedtools.BedTool.from_dataframe(anchor1).intersect(H2AZ_plusCL_reprod_peaks_bed, wa = True).to_dataframe().drop_duplicates()
anchor2_h2az = pybedtools.BedTool.from_dataframe(anchor2).intersect(H2AZ_plusCL_reprod_peaks_bed, wa = True).to_dataframe().drop_duplicates()

shNT_plusCL_allLoop_up_down_df['anchor1_h2az'] = shNT_plusCL_allLoop_up_down_df['anchor1'].isin(anchor1_h2az['name'])
shNT_plusCL_allLoop_up_down_df['anchor2_h2az'] = shNT_plusCL_allLoop_up_down_df['anchor2'].isin(anchor2_h2az['name'])

## separate loop set into H2AZ occupied and H2AZ non-occupied, because latter analysis will compare them
shNT_plusCL_allLoop_up_down_H2AZ_occ = shNT_plusCL_allLoop_up_down_df.query('anchor1_h2az == True or anchor2_h2az == True')
shNT_plusCL_allLoop_up_down_H2AZ_nonocc = shNT_plusCL_allLoop_up_down_df.query('anchor1_h2az == False and anchor2_h2az == False')


## ---- extract loop connectivity result between H2AZ occupied and H2AZ non-occupied loop sets -----
## H2AZ non-occupied loops but connected with H2AZ occupied and depandant loops
connected_loops = []
connected_loops.extend(shNT_plusCL_loop_connect_res[shNT_plusCL_loop_connect_res['loop1'].isin(shNT_plusCL_allLoop_up_down_H2AZ_occ.index)]['loop2'].tolist())
connected_loops.extend(shNT_plusCL_loop_connect_res[shNT_plusCL_loop_connect_res['loop2'].isin(shNT_plusCL_allLoop_up_down_H2AZ_occ.index)]['loop1'].tolist())
### sometime, connectivity happens among many loops, so iterate to find all connected loops
stop = False
while stop != True:
    new_connect = shNT_plusCL_loop_connect_res[shNT_plusCL_loop_connect_res['loop1'].isin(connected_loops)]['loop2'].tolist()
    new_connect.extend(shNT_plusCL_loop_connect_res[shNT_plusCL_loop_connect_res['loop2'].isin(connected_loops)]['loop1'].tolist())
    new_connect = np.setdiff1d(new_connect, connected_loops)
    print(len(new_connect))
    if len(new_connect)>0:
        connected_loops.extend(new_connect)
    else:
        stop = True ## stop when no new loop connected

connected_loops = list(set(connected_loops) - set(shNT_plusCL_allLoop_up_down_H2AZ_occ.index)) ## connected but not H2AZ occu and depand
print('Total connected loops: %s'%len(connected_loops))
print('fraction of H2AZ non-occ loop:', len(np.intersect1d(connected_loops, shNT_plusCL_allLoop_up_down_H2AZ_nonocc.index)) / len(shNT_plusCL_allLoop_up_down_H2AZ_nonocc))


# In[ ]:


## visualize contact frequecy of H2AZ non-occupied but connected to H2AZ occupied loops
dot_df = loop_dots[['score_shNT_plusCL', 'score_shNT_minusCL', 'score_KD_plusCL']]
dot_df = pd.DataFrame(quantile_norm(dot_df), 
                   columns = ['score_shNT_plusCL', 'score_shNT_minusCL', 'score_KD_plusCL'], index = dot_df.index.tolist())
dot_df['KD_vs_shNT'] = dot_df['score_KD_plusCL'] - dot_df['score_shNT_plusCL']
dot_df['plusCL_vs_minusCL'] = dot_df['score_shNT_plusCL'] - dot_df['score_shNT_minusCL']

plot_df = dot_df[dot_df.index.isin(connected_loops)][['score_shNT_plusCL', 'score_shNT_minusCL', 'score_KD_plusCL']]
fig, ax = plt.subplots(figsize = (4, 5))
sns.boxplot(plot_df.melt(), x = 'variable', y = 'value', hue = 'variable', 
            order = ['score_shNT_minusCL', 'score_shNT_plusCL', 'score_KD_plusCL'],
            palette={'score_shNT_minusCL': 'grey', 'score_shNT_plusCL':plt.cm.get_cmap("tab10")(1),
                     'score_KD_plusCL':plt.cm.get_cmap("tab10")(0)},
            ax = ax, showfliers = False, whis = 0.8)
ax.tick_params(axis = 'x', rotation = 75)
ax.set(xlabel='', ylabel='loop score',
       title = 'H2AZ-occupied connected loop')
sns.despine()
plt.tight_layout()
plt.show()
plt.close()


## random sampling for background control
np.random.seed(1234)
random_loop = np.random.choice(dot_df.index, len(connected_loops))

plot_df = dot_df[dot_df.index.isin(random_loop)][['score_shNT_plusCL', 'score_shNT_minusCL', 'score_KD_plusCL']]
fig, ax = plt.subplots(figsize = (4, 5))
sns.boxplot(plot_df.melt(), x = 'variable', y = 'value', hue = 'variable', 
            order = ['score_shNT_minusCL', 'score_shNT_plusCL', 'score_KD_plusCL'],
            palette={'score_shNT_minusCL': 'grey', 'score_shNT_plusCL':plt.cm.get_cmap("tab10")(1),
                     'score_KD_plusCL':plt.cm.get_cmap("tab10")(0)},
            ax = ax, showfliers = False, whis = .8)
ax.tick_params(axis = 'x', rotation = 75)
ax.set(xlabel='', ylabel='loop score',
       title = 'Random')
sns.despine()
plt.tight_layout()
plt.show()
plt.close()



# ### H3K27ac for enhancer activity
# 

# In[ ]:


### check H3K27ac signal around H2AZ sites in loops and compare between plusCL and minusCL group
## H3K27ac signal from ChIP-seq bigwig
k27_cl_bw = {'plusCL': '../DataProcess/ChIP_H3K27ac/20230428_ChIPseq_Day11/H3K27ac_plusCL.rep1/H3K27ac_plusCL.rep1_treat_pileup.bw',
         'minusCL': '../DataProcess/ChIP_H3K27ac/20230428_ChIPseq_Day11/H3K27ac_minusCL.rep1/H3K27ac_minusCL.rep1_treat_pileup.bw'}
k27_kd_bw = {'shNT_plusCL_rep1': '../DataProcess/ChIP_H3K27ac/20230616_ChIP_H2AZ_Histone/shNT_plusCL_H3K27ac.rep1/shNT_plusCL_H3K27ac.rep1_treat_pileup.bw',
         'shNT_plusCL_rep2': '../DataProcess/ChIP_H3K27ac/20230616_ChIP_H2AZ_Histone/shNT_plusCL_H3K27ac.rep2/shNT_plusCL_H3K27ac.rep2_treat_pileup.bw',
            'KD_plusCL_rep1': '../DataProcess/ChIP_H3K27ac/20230616_ChIP_H2AZ_Histone/KD_plusCL_H3K27ac.rep1/KD_plusCL_H3K27ac.rep1_treat_pileup.bw',
         'KD_plusCL_rep2': '../DataProcess/ChIP_H3K27ac/20230616_ChIP_H2AZ_Histone/KD_plusCL_H3K27ac.rep2/KD_plusCL_H3K27ac.rep2_treat_pileup.bw'}

## extract H3K27ac signal in plusCL and minusCL samples
k27ac_cl_bw_h2az_site_values={}
for t, r, l in [['EP', EP_res, CL_induce_KD_down_EP],
         ['PP', PP_res, CL_induce_KD_down_PP],
         ['PO', PO_res, CL_induce_KD_down_PO],
         ['OO', OO_res, CL_induce_KD_down_OO]]:
    H2AZ_bind_loops = pd.concat([r['H2AZ_peak']['loops']['either'], r['H2AZ_peak']['loops']['both']])
    H2AZ_bind_loops['anchor1'] = H2AZ_bind_loops['chrom1']+':'+H2AZ_bind_loops['start1'].astype('str')+'-'+H2AZ_bind_loops['end1'].astype('str')
    H2AZ_bind_loops['anchor2'] = H2AZ_bind_loops['chrom2']+':'+H2AZ_bind_loops['start2'].astype('str')+'-'+H2AZ_bind_loops['end2'].astype('str')
   
    bind_anchors = pd.DataFrame(H2AZ_bind_loops[['chrom1', 'start1', 'end1']].values.tolist()+H2AZ_bind_loops[['chrom2', 'start2', 'end2']].values.tolist())
    bind_anchors_bed = pybedtools.BedTool.from_dataframe(bind_anchors)
    H2AZ_site = H2AZ_plusCL_reprod_peaks_bed.intersect(bind_anchors_bed, wa = True).to_dataframe().drop_duplicates()
    H2AZ_site['start'] = H2AZ_site['start'] + ((H2AZ_site['start']-H2AZ_site['end'])/2).abs().astype('int')
    
    d = 5000
    k27ac_cl_bw_h2az_site_values[t] = []
    for path in k27_cl_bw:
        bw = pyBigWig.open(k27_cl_bw[path])
        for i, line in H2AZ_site.iterrows():
            line['bw'] = bw.values(line['chrom'], line['start']-d, line['start']+d)
            k27ac_cl_bw_h2az_site_values[t].append(line.values.tolist()+[path])

## E-P loops were focused because for checking enhancer activity
ddf = pd.DataFrame(k27ac_cl_bw_h2az_site_values['EP'])
ddf['mean'] = [np.mean(x[4500:5500]) for x in ddf[8].tolist()]

fig, ax = plt.subplots(figsize = (3,3.5))
sns.boxplot(data = ddf, x = 9, y = 'mean', 
            showfliers = False, order = ['minusCL', 'plusCL'], width = .5, 
            palette = {'plusCL': 'darkorange', 'minusCL': 'grey'})
ax.set(xlabel = '', ylabel = 'H3K27ac ChIP-seq', title = 'All K27ac in occu loop (1000bp)')
sns.despine()
plt.show()
plt.close()
print('plusCL vs minusCL p value: ', wilcoxon(ddf[ddf[9] == "plusCL"]['mean'], ddf[ddf[9] == "minusCL"]['mean']))


# In[ ]:


## specifically focus on up-regulated H3K27ac peaks 
def _volcano_(dat, pcol = 'height_diff_FDR', fcol = 'height_log2FC', pcut = .05, fcut = 0.58,
             xlabel = 'Log2 Fold Change', ylabel = '-Log10(FDR)', title = '', pdf = None):
    pdat = dat[[pcol, fcol]]
    label = []
    for p, f in pdat[[pcol, fcol]].values.tolist():
        if f > fcut and p < pcut:
            label.append('Up')
        elif f < -fcut and p < pcut:
            label.append('Down')
        else:
            label.append('No')
    pdat['label'] = label
    pdat = pdat.sort_values('label')
    pdat['-log10(significane)'] = -np.log10(pdat[pcol])
    
    fig, ax = plt.subplots()
    sns.scatterplot(data = pdat, x = fcol, y = '-log10(significane)', hue = 'label', 
                    palette={'Up':'#FF4500', 'Down':'#1E90FF', 'No': 'lightgrey'},
                    edgecolor = 'none', s = 3, alpha = .8)
    ax.legend(title = '', markerscale=3, frameon = False, 
              bbox_to_anchor=[1.1, 0.5], loc='center')
    up_set = dat.loc[pdat.query('label == "Up"').index,:]
    down_set = dat.loc[pdat.query('label == "Down"').index,:]
    title = title+'\n'+'Up: %s, Down: %s'%(up_set.shape[0], down_set.shape[0])
    ax.set(xlabel = xlabel, ylabel = ylabel, title = title)
    sns.despine()
    plt.tight_layout()
    plt.show()
    pdf.savefig(fig) if pdf else None
    plt.close()
    return(up_set, down_set, label)
    
### DANPOS was used to run differential H3K27ac peaks
k27ac_diff_cl = pd.read_csv('../DataProcess/ChIP_H3K27ac/plusCL_vs_minus_CL/plusCL_K27_vs_IgG_pooled_20230428_ChIPseq_Day11_H3K27ac_plusCL_unique.sorted.bgsub.Fnor-minusCL_K27_vs_IgG_pooled_20230428_ChIPseq_Day11_H3K27ac_minusCL_unique.sorted.bgsub.Fnor.peaks.integrative.xls', sep = '\t')

## focus on differential H3K27ac peaks in EP loop anchors
k27ac_diff_cl_bed1 = pybedtools.BedTool.from_dataframe(k27ac_diff_cl.iloc[:,:3].reset_index()[['chr', 'start', 'end', 'index']])
H2AZ_bind_loops = pd.concat([EP_res['H2AZ_peak']['loops']['either'], EP_res['H2AZ_peak']['loops']['both']])
H2AZ_bind_loops['anchor1'] = H2AZ_bind_loops['chrom1']+':'+H2AZ_bind_loops['start1'].astype('str')+'-'+H2AZ_bind_loops['end1'].astype('str')
H2AZ_bind_loops['anchor2'] = H2AZ_bind_loops['chrom2']+':'+H2AZ_bind_loops['start2'].astype('str')+'-'+H2AZ_bind_loops['end2'].astype('str')
bind_anchors = pd.DataFrame(H2AZ_bind_loops[['chrom1', 'start1', 'end1']].values.tolist()+H2AZ_bind_loops[['chrom2', 'start2', 'end2']].values.tolist())
bind_anchors_bed = pybedtools.BedTool.from_dataframe(bind_anchors)
## differential H3K27ac peaks overlapped with H2AZ sites were focused
tmp = k27ac_diff_cl_bed1.intersect(H2AZ_plusCL_reprod_peaks_bed.intersect(bind_anchors_bed, wa = True), wa = True).to_dataframe().drop_duplicates()
## significanly up k27ac peaks
tmp_up = k27ac_diff_cl.loc[tmp['name'],].query('height_log2FC > 0')
tmp_up_k27ac = k27ac_union_peaks_bed.intersect(pybedtools.BedTool.from_dataframe(tmp_up.iloc[:,:3]),
                                              wa = True).to_dataframe().drop_duplicates()


# In[ ]:


## extract bigwig signal for plusCL and minusCL samples
tmp_up_k27_plusCL_minusCL = _bw_peaks_(bw_dict=k27_cl_bw, peak_site1=tmp_up_k27ac, peak_site2=None,
                               colmap = {'plusCL': 'darkorange', 'minusCL': 'grey'}, figsize = (4.5,3))
tmp_up_k27_plusCL_minusCL_df = pd.DataFrame(tmp_up_k27_plusCL_minusCL)
tmp_up_k27_plusCL_minusCL_df['peak'] = tmp_up_k27_plusCL_minusCL_df.iloc[:,:3].apply(lambda row: '-'.join(row.astype('str').tolist()), axis = 1)
tmp_up_k27_plusCL_minusCL_df['mean'] = [np.median(x[1500:2500]) for x in tmp_up_k27_plusCL_minusCL_df[3].tolist()]
## plot
fig, ax = plt.subplots(figsize = (3,3.5))
sns.boxplot(data = tmp_up_k27_plusCL_minusCL_df, x = 4, y = 'mean', order = ['minusCL', 'plusCL'],
            showfliers = False, width = .7, palette = {'plusCL': 'darkorange', 'minusCL': 'grey'})
ax.set(xlabel = '', ylabel = 'H3K27ac ChIP-seq', title = 'up K27ac in occu loop (1000bp)')
sns.despine()
plt.show()
plt.close()
print('pvalue: ', wilcoxon(tmp_up_k27_plusCL_minusCL_df[tmp_up_k27_plusCL_minusCL_df[4] == "plusCL"]['mean'],
         tmp_up_k27_plusCL_minusCL_df[tmp_up_k27_plusCL_minusCL_df[4] == "minusCL"]['mean'])
)


# In[ ]:


## extract bigwig signal for H2AZ KD plusCL and control plusCL samples
tmp_up_k27_KD_shNT = _bw_peaks_(bw_dict=k27_kd_bw, peak_site1=tmp_up_k27ac, peak_site2=None,
                               colmap = {'shNT_plusCL': 'darkorange', 'KD_plusCL': plt.cm.get_cmap('tab10')(0)}, figsize = (4.5,3))
tmp_up_k27_KD_shNT_df = pd.DataFrame(tmp_up_k27_KD_shNT)
tmp_up_k27_KD_shNT_df['peak'] = tmp_up_k27_KD_shNT_df.iloc[:,:3].apply(lambda row: '-'.join(row.astype('str').tolist()), axis = 1)
tmp_up_k27_KD_shNT_df['cond'] = tmp_up_k27_KD_shNT_df[4].str.replace('_rep1', '').str.replace('_rep2', '')
tmp_up_k27_KD_shNT_df = tmp_up_k27_KD_shNT_df.groupby(['peak', 'cond']).apply(lambda d: pd.DataFrame(d[3].tolist()).mean())#.reset_index()
tmp_up_k27_KD_shNT_df_plot = tmp_up_k27_KD_shNT_df.iloc[:,1500:2500].T.mean().reset_index()
tmp_up_k27_KD_shNT_df_plot.columns = ['peak', 'cond', 'value']
## plot
fig, ax = plt.subplots(figsize = (3,3.5))
sns.boxplot(data = tmp_up_k27_KD_shNT_df_plot, x = 'cond', y = 'value', showfliers = False, width = .7,
           order = ['shNT_plusCL', 'KD_plusCL'], palette={'shNT_plusCL': 'darkorange', 'KD_plusCL': plt.cm.get_cmap('tab10')(0)})
ax.set(xlabel = '', ylabel = 'H3K27ac ChIP-seq', title = 'up K27ac in occu loop (1000bp)')
sns.despine()
plt.show()
plt.close()
print('pvalue:', wilcoxon(tmp_up_k27_KD_shNT_df_plot[tmp_up_k27_KD_shNT_df_plot['cond'] == "shNT_plusCL"]['value'],
         tmp_up_k27_KD_shNT_df_plot[tmp_up_k27_KD_shNT_df_plot['cond'] == "KD_plusCL"]['value'])
)


# ### CTCF, SMC1, MED1, and ATAC-seq integration
# 

# In[ ]:


ctcf_bw = {'shNT_plusCL_rep1': '../DataProcess/ChIP_CTCF/08212024/NT_plusCL_CTCF_run.rep1/NT_plusCL_CTCF_run.rep1_treat_pileup.bw',
          'shNT_plusCL_rep2': '../DataProcess/ChIP_CTCF/11042024/NT_plusCL_CTCF.rep1/NT_plusCL_CTCF.rep1_treat_pileup.bw',
          'shNT_minusCL_rep1': '../DataProcess/ChIP_CTCF/08212024/NT_minusCL_CTCF_run.rep1/NT_minusCL_CTCF_run.rep1_treat_pileup.bw',
          'shNT_minusCL_rep2': '../DataProcess/ChIP_CTCF/11042024/NT_minusCL_CTCF.rep1/NT_minusCL_CTCF.rep1_treat_pileup.bw'}

smc1_bw = {'shNT_plusCL_rep1': '../DataProcess/ChIP_SMC1/shNT_plusCL_SMC1_r1.rep1_treat_pileup.bw',
          'shNT_plusCL_rep2': '../DataProcess/ChIP_SMC1/shNT_plusCL_SMC1_r2.rep1_treat_pileup.bw',
          'shNT_minusCL_rep1': '../DataProcess/ChIP_SMC1/shNT_minusCL_SMC1_r1.rep1_treat_pileup.bw',
          'shNT_minusCL_rep2': '../DataProcess/ChIP_SMC1/shNT_minusCL_SMC1_r2.rep1_treat_pileup.bw'}

med1_bw = {'shNT_plusCL_rep1': '../DataProcess/ChIP_MED1/shNT_plusCL_MED1_r1.rep1_treat_pileup.bw',
          'shNT_plusCL_rep2': '../DataProcess/ChIP_MED1/shNT_plusCL_MED1_r2.rep1_treat_pileup.bw',
          'shNT_minusCL_rep1': '../DataProcess/ChIP_MED1/shNT_minusCL_MED1_r1.rep1_treat_pileup.bw',
          'shNT_minusCL_rep2': '../DataProcess/ChIP_MED1/shNT_minusCL_MED1_r2.rep1_treat_pileup.bw',
          'KD_plusCL_rep1': '../DataProcess/ChIP_MED1/KD_plusCL_MED1_r1.rep1_treat_pileup.bw',
          'KD_plusCL_rep2': '../DataProcess/ChIP_MED1/KD_plusCL_MED1_r2.rep1_treat_pileup.bw'}

atac_bw = {'shNT_plusCL_rep1': '../DataProcess/ATAC/bw/shNT_plusCL_rep2.rep1_treat_pileup.bw',
          'shNT_plusCL_rep2': '../DataProcess/ATAC/bw/shNT_plusCL_rep3.rep1_treat_pileup.bw',
          'shNT_minusCL_rep1': '../DataProcess/ATAC/bw/shNT_minusCL_rep2.rep1_treat_pileup.bw',
          'shNT_minusCL_rep2': '../DataProcess/ATAC/bw/shNT_minusCL_rep3.rep1_treat_pileup.bw',
          'KD_plusCL_rep1': '../DataProcess/ATAC/bw/KD_plusCL_rep2.rep1_treat_pileup.bw',
          'KD_plusCL_rep2': '../DataProcess/ATAC/bw/KD_plusCL_rep3.rep1_treat_pileup.bw'}


# #### CTCF and SMC1

# In[ ]:


## CTCF_SMC1 peaks in different loop sets
CTCF_SMC1_peaks_H2AZ_occu_loop = CTCF_SMC1_peaks_plusCL_shNT_bed.intersect(
    pybedtools.BedTool.from_dataframe(H2AZ_bind_loops_anchors), wa = True
).to_dataframe().drop_duplicates()

CTCF_SMC1_peaks_up_loop = CTCF_SMC1_peaks_plusCL_shNT_bed.intersect(
    pybedtools.BedTool.from_dataframe(pd.DataFrame(up_loop.iloc[:,0:3].values.tolist()+up_loop.iloc[:,3:6].values.tolist())), 
    wa = True
).to_dataframe().drop_duplicates()

CTCF_SMC1_peaks_stable_loop = CTCF_SMC1_peaks_plusCL_shNT_bed.intersect(
    pybedtools.BedTool.from_dataframe(pd.DataFrame(stable_loop_union.iloc[:,0:3].values.tolist()+stable_loop_union.iloc[:,3:6].values.tolist())), 
    wa = True
).to_dataframe().drop_duplicates()

## CTCF SMC1 peaks in both occupied loop
CTCF_SMC1_peaks_H2AZ_occu_loop_bothanchor = CTCF_SMC1_peaks_plusCL_shNT_bed.intersect(
    pybedtools.BedTool.from_dataframe(pd.DataFrame(CTCF_SMC1_H2AZ_bind_loop['both'].iloc[:,0:3].values.tolist()+
                                     CTCF_SMC1_H2AZ_bind_loop['both'].iloc[:,3:6].values.tolist())), wa = True
).to_dataframe().drop_duplicates()

CTCF_SMC1_peaks_up_loop_bothanchor = CTCF_SMC1_peaks_plusCL_shNT_bed.intersect(
    pybedtools.BedTool.from_dataframe(pd.DataFrame(CTCF_SMC1_up_loop['both'].iloc[:,0:3].values.tolist()+
                                     CTCF_SMC1_up_loop['both'].iloc[:,3:6].values.tolist())), 
    wa = True
).to_dataframe().drop_duplicates()

CTCF_SMC1_peaks_stable_loop_bothanchor = CTCF_SMC1_peaks_plusCL_shNT_bed.intersect(
    pybedtools.BedTool.from_dataframe(pd.DataFrame(CTCF_SMC1_stable_loop['both'].iloc[:,0:3].values.tolist()+
                                     CTCF_SMC1_stable_loop['both'].iloc[:,3:6].values.tolist())), 
    wa = True
).to_dataframe().drop_duplicates()



# In[ ]:


## CTCF signal was not change at CTCF site in H2AZ occu loop in overall between plusCL and minusCL
bw_CTCF_site_H2AZ_loop_values=[]
d = 2000
for path in ctcf_bw:
    bw = pyBigWig.open(ctcf_bw[path])
    for i, line in CTCF_SMC1_peaks_H2AZ_occu_loop.iterrows():
        c = line['start']+int((line['end']-line['start'])/2)
        line['bw'] = bw.values(line['chrom'], c-d, c+d)
        bw_CTCF_site_H2AZ_loop_values.append(line.values.tolist()+[path])

## plot
df = pd.DataFrame(bw_CTCF_site_H2AZ_loop_values)
df['value'] = [np.mean(x[1500:2500]) for x in df[10].tolist()]
df['cond'] = [x.replace('_rep1', '').replace('_rep2', '') for x in df[11].tolist()]

fig, ax = plt.subplots(figsize = (4, 4))
sns.violinplot(data = df.query('cond != "KD_plusCL"'), x = 'cond', y = 'value',
               cut = 0, width = .5, palette={'shNT_minusCL':'grey', 'shNT_plusCL':'darkorange'}, 
               order = ['shNT_minusCL', 'shNT_plusCL'])
ax.set(xlabel = '', ylabel = 'CTCF ChIP-seq signal', title = 'CTCF_SMC1 peak (500bp)')
sns.despine()
plt.tight_layout()
plt.show()
plt.close()
print('p values: ', ranksums(df.query('cond == "shNT_plusCL"')['value'], df.query('cond == "shNT_minusCL"')['value']))


## SMC1 signal was not change at CTCF site in H2AZ occu loop in overall between plusCL and minusCL
smc1_bw_CTCF_site_H2AZ_loop_values=[]
d = 2000
for path in smc1_bw:
    bw = pyBigWig.open(smc1_bw[path])
    for i, line in CTCF_SMC1_peaks_H2AZ_occu_loop.iterrows():
        c = line['start']+int((line['end']-line['start'])/2)
        line['bw'] = bw.values(line['chrom'], c-d, c+d)
        smc1_bw_CTCF_site_H2AZ_loop_values.append(line.values.tolist()+[path])
print('p values: ', ranksums(df.query('cond == "shNT_plusCL"')['value'], df.query('cond == "shNT_minusCL"')['value']))



# In[ ]:


## CTCF ChIP-seq signal in CTCF_SMC1 copeak sites were lower in up loop anchor than stable loops
ctcf_bw_CTCF_site_loop_both_anchor_values=[]
d = 250
for path in ctcf_bw:
    bw = pyBigWig.open(ctcf_bw[path])

    for i, line in CTCF_SMC1_peaks_up_loop_bothanchor.iterrows():
        c = line['start']+int((line['end']-line['start'])/2)
        line['bw'] = bw.stats(line['chrom'], c-d, c+d, exact = True)
        ctcf_bw_CTCF_site_loop_both_anchor_values.append(line.values.tolist()+[path, 'Up_loop'])

    for i, line in CTCF_SMC1_peaks_stable_loop_bothanchor.iterrows():
        c = line['start']+int((line['end']-line['start'])/2)
        line['bw'] = bw.stats(line['chrom'], c-d, c+d, exact = True)
        ctcf_bw_CTCF_site_loop_both_anchor_values.append(line.values.tolist()+[path, 'Stable'])

ctcf_bw_CTCF_site_loop_both_anchor_values = pd.DataFrame(ctcf_bw_CTCF_site_loop_both_anchor_values)
ctcf_bw_CTCF_site_loop_both_anchor_values[10] = [x[0] for x in ctcf_bw_CTCF_site_loop_both_anchor_values[10].tolist()]
ctcf_bw_CTCF_site_loop_both_anchor_values['cond'] = ctcf_bw_CTCF_site_loop_both_anchor_values[11].str.replace('_rep1', '').str.replace('_rep2', '')
ctcf_bw_CTCF_site_loop_both_anchor_values['rep'] = [x.split('_')[-1] for x in ctcf_bw_CTCF_site_loop_both_anchor_values[11].tolist()]
## make a dataframe for plot
plot_df = ctcf_bw_CTCF_site_loop_both_anchor_values.groupby([3, 12, 'cond'])[10].mean().reset_index()
plot_df.columns = ['peak', 'Type', 'cond', 'value']

pplot_df = plot_df.query('cond == "shNT_plusCL"')
fig, ax = plt.subplots(figsize = (3,3.5))
sns.boxplot(data = pplot_df,
            x = 'Type', y = 'value', showfliers = False, color = 'lightblue')
ax.set(ylabel = 'CTCF ChIP-seq signal', xlabel = 'loop type', title = 'CTCF_SMC1 peak (500bp)')
plt.tight_layout()
plt.show()
plt.close()
print('p values: ', ranksums(pplot_df.query('Type == "Up_loop"')['value'], pplot_df.query('Type == "Stable"')['value']))



# In[ ]:


## SMC1 ChIP-seq signal in CTCF_SMC1 copeak sites were lower in up loop anchor than stable loops
smc1_bw_CTCF_site_loop_both_anchor_values=[]
d = 250
for path in smc1_bw:
    bw = pyBigWig.open(smc1_bw[path])
    for i, line in CTCF_SMC1_peaks_H2AZ_occu_loop_bothanchor.iterrows():
        c = line['start']+int((line['end']-line['start'])/2)
        line['bw'] = bw.stats(line['chrom'], c-d, c+d, exact = True)
        smc1_bw_CTCF_site_loop_both_anchor_values.append(line.values.tolist()+[path, 'H2AZ_occu'])
        
    for i, line in CTCF_SMC1_peaks_up_loop_bothanchor.iterrows():
        c = line['start']+int((line['end']-line['start'])/2)
        line['bw'] = bw.stats(line['chrom'], c-d, c+d, exact = True)
        smc1_bw_CTCF_site_loop_both_anchor_values.append(line.values.tolist()+[path, 'Up_loop'])

    for i, line in CTCF_SMC1_peaks_stable_loop_bothanchor.iterrows():
        c = line['start']+int((line['end']-line['start'])/2)
        line['bw'] = bw.stats(line['chrom'], c-d, c+d, exact = True)
        smc1_bw_CTCF_site_loop_both_anchor_values.append(line.values.tolist()+[path, 'Stable'])

smc1_bw_CTCF_site_loop_both_anchor_values = pd.DataFrame(smc1_bw_CTCF_site_loop_both_anchor_values)
smc1_bw_CTCF_site_loop_both_anchor_values[10] = [x[0] for x in smc1_bw_CTCF_site_loop_both_anchor_values[10].tolist()]
smc1_bw_CTCF_site_loop_both_anchor_values['cond'] = smc1_bw_CTCF_site_loop_both_anchor_values[11].str.replace('_rep1', '').str.replace('_rep2', '')
smc1_bw_CTCF_site_loop_both_anchor_values['rep'] = [x.split('_')[-1] for x in smc1_bw_CTCF_site_loop_both_anchor_values[11].tolist()]
## make a dataframe for plot
plot_df = smc1_bw_CTCF_site_loop_both_anchor_values.groupby([3, 12, 'cond'])[10].mean().reset_index()
plot_df.columns = ['peak', 'Type', 'cond', 'value']

pplot_df = plot_df.query('cond == "shNT_plusCL"')
fig, ax = plt.subplots(figsize = (3,3.5))
sns.boxplot(data = pplot_df,
            x = 'Type', y = 'value', showfliers = False, color = 'lightblue')
ax.set(ylabel = 'SMC1 ChIP-seq signal', xlabel = 'loop type', title = 'CTCF_SMC1 peak (500bp)')
plt.tight_layout()
plt.show()
plt.close()
print('p values: ', ranksums(pplot_df.query('Type == "Up_loop"')['value'], pplot_df.query('Type == "Stable"')['value']))


# #### MED1 and ATAC

# In[ ]:


### this analysis focus on MED1 signal change vs ATAC-seq signal change at ATAC site in loop anchors

## get union ATAC-seq peaks across conditions for comparison analysis between plusCL vs minusCL and KD vs Control plusCL
ATAC_union_peaks = pd.concat([ATAC_peaks_reprod_bed['shNT_plusCL'].to_dataframe(),
                              ATAC_peaks_reprod_bed['shNT_minusCL'].to_dataframe(),
                              ATAC_peaks_reprod_bed['KD_plusCL'].to_dataframe()]).drop_duplicates()
ATAC_union_peaks_bed = pybedtools.BedTool.from_dataframe(ATAC_union_peaks)


### ATAC-seq signal at H2AZ site in H2AZ occupied loops in three conditions
ATAC_bw_ATAC_site_values={}
for t, r, l in [['EP', EP_res, CL_induce_KD_down_EP],
         ['PP', PP_res, CL_induce_KD_down_PP],
         ['PO', PO_res, CL_induce_KD_down_PO],
         ['OO', OO_res, CL_induce_KD_down_OO]]:
    H2AZ_bind_loops = pd.concat([r['H2AZ_peak']['loops']['either'], r['H2AZ_peak']['loops']['both']])
    H2AZ_bind_loops['anchor1'] = H2AZ_bind_loops['chrom1']+':'+H2AZ_bind_loops['start1'].astype('str')+'-'+H2AZ_bind_loops['end1'].astype('str')
    H2AZ_bind_loops['anchor2'] = H2AZ_bind_loops['chrom2']+':'+H2AZ_bind_loops['start2'].astype('str')+'-'+H2AZ_bind_loops['end2'].astype('str')
   
    bind_anchors = pd.DataFrame(H2AZ_bind_loops[['chrom1', 'start1', 'end1']].values.tolist()+H2AZ_bind_loops[['chrom2', 'start2', 'end2']].values.tolist())
    bind_anchors_bed = pybedtools.BedTool.from_dataframe(bind_anchors)
    ATAC_site = ATAC_union_peaks_bed.intersect(bind_anchors_bed, wa = True).to_dataframe().drop_duplicates()
    ATAC_site['start'] = ATAC_site['start'] + ((ATAC_site['start']-ATAC_site['end'])/2).abs().astype('int')
    
    d = 3000
    ATAC_bw_ATAC_site_values[t] = []
    for path in atac_bw:
        bw = pyBigWig.open(atac_bw[path])
        for i, line in ATAC_site.iterrows():
            line['bw'] = bw.values(line['chrom'], line['start']-d, line['start']+d)
            ATAC_bw_ATAC_site_values[t].append(line.values.tolist()+[path])

### MED1 signal at H2AZ site in H2AZ occupied loops in three conditions
MED1_bw_ATAC_site_values={}
for t, r, l in [['EP', EP_res, CL_induce_KD_down_EP],
         ['PP', PP_res, CL_induce_KD_down_PP],
         ['PO', PO_res, CL_induce_KD_down_PO],
         ['OO', OO_res, CL_induce_KD_down_OO]]:
    H2AZ_bind_loops = pd.concat([r['H2AZ_peak']['loops']['either'], r['H2AZ_peak']['loops']['both']])
    H2AZ_bind_loops['anchor1'] = H2AZ_bind_loops['chrom1']+':'+H2AZ_bind_loops['start1'].astype('str')+'-'+H2AZ_bind_loops['end1'].astype('str')
    H2AZ_bind_loops['anchor2'] = H2AZ_bind_loops['chrom2']+':'+H2AZ_bind_loops['start2'].astype('str')+'-'+H2AZ_bind_loops['end2'].astype('str')
   
    bind_anchors = pd.DataFrame(H2AZ_bind_loops[['chrom1', 'start1', 'end1']].values.tolist()+H2AZ_bind_loops[['chrom2', 'start2', 'end2']].values.tolist())
    bind_anchors_bed = pybedtools.BedTool.from_dataframe(bind_anchors)
    ATAC_site = ATAC_union_peaks_bed.intersect(bind_anchors_bed, wa = True).to_dataframe().drop_duplicates()
    ATAC_site['start'] = ATAC_site['start'] + ((ATAC_site['start']-ATAC_site['end'])/2).abs().astype('int')
    
    d = 3000
    MED1_bw_ATAC_site_values[t] = []
    for path in med1_bw:
        bw = pyBigWig.open(med1_bw[path])
        for i, line in ATAC_site.iterrows():
            line['bw'] = bw.values(line['chrom'], line['start']-d, line['start']+d)
            MED1_bw_ATAC_site_values[t].append(line.values.tolist()+[path])

### make bw signal into a dataframe for plotting
med1_df = pd.concat([pd.DataFrame(MED1_bw_ATAC_site_values['EP']),
               pd.DataFrame(MED1_bw_ATAC_site_values['PP']),
               pd.DataFrame(MED1_bw_ATAC_site_values['PO']),
               pd.DataFrame(MED1_bw_ATAC_site_values['OO'])])
med1_df['peaks'] = med1_df.iloc[:,0:3].apply(lambda row: '-'.join(row.astype('str').tolist()), axis = 1).tolist()
med1_df['signal'] = [np.mean(x) for x in med1_df[8].tolist()]
med1_df['cond'] = [x.replace('_rep1', '').replace('_rep2', '') for x in med1_df[9].tolist()]

med1_mean_df = med1_df.groupby(['peaks', 'cond'])['signal'].mean().reset_index()
med1_mean_df = med1_mean_df.pivot_table(index = 'peaks', columns = 'cond', values = 'signal')


atac_df = pd.concat([pd.DataFrame(ATAC_bw_ATAC_site_values['EP']),
               pd.DataFrame(ATAC_bw_ATAC_site_values['PP']),
               pd.DataFrame(ATAC_bw_ATAC_site_values['PO']),
               pd.DataFrame(ATAC_bw_ATAC_site_values['OO'])])
atac_df['peaks'] = atac_df.iloc[:,0:3].apply(lambda row: '-'.join(row.astype('str').tolist()), axis = 1).tolist()
atac_df['signal'] = [np.mean(x) for x in atac_df[8].tolist()]
atac_df['cond'] = [x.replace('_rep1', '').replace('_rep2', '') for x in atac_df[9].tolist()]

atac_mean_df = atac_df.groupby(['peaks', 'cond'])['signal'].mean().reset_index()
atac_mean_df = atac_mean_df.pivot_table(index = 'peaks', columns = 'cond', values = 'signal')


### correlation between delta ATAC signal and delta MED1 signal upon CL treatment
print(pearsonr(atac_mean_df['shNT_plusCL']-atac_mean_df['shNT_minusCL'],
         med1_mean_df['shNT_plusCL']-med1_mean_df['shNT_minusCL']))

fig, ax = plt.subplots(figsize = (4,4))
sns.regplot(x = atac_mean_df['shNT_plusCL']-atac_mean_df['shNT_minusCL'],
            y = med1_mean_df['shNT_plusCL']-med1_mean_df['shNT_minusCL'],
            scatter_kws = {'s':10, 'color': 'black'}, line_kws={'color': 'darkorange', 'linestyle': 'dashed'})
ax.set(xlabel = 'Delta ATAC-seq signal', ylabel = 'Delta MED1 ChIP-seq signal', 
       title = 'ATAC peak, plusCL vs minusCL')
plt.show()
plt.close()

### correlation between delta ATAC signal and delta MED1 signal upon H2AZ KD in CL treatment condition
print(pearsonr(atac_mean_df['KD_plusCL']-atac_mean_df['shNT_plusCL'],
         med1_mean_df['KD_plusCL']-med1_mean_df['shNT_plusCL']))

fig, ax = plt.subplots(figsize = (4,4))
sns.regplot(x = atac_mean_df['KD_plusCL']-atac_mean_df['shNT_plusCL'],
         y = med1_mean_df['KD_plusCL']-med1_mean_df['shNT_plusCL'],
            scatter_kws = {'s':10, 'color': 'black'}, line_kws={'color': plt.cm.get_cmap('tab10')(0), 'linestyle': 'dashed'})
ax.set(xlabel = 'Delta ATAC-seq signal', ylabel = 'Delta MED1 ChIP-seq signal', title = 'ATAC peak, KD vs NT')
ax.set_ylim(-0.8, 0.12) ## set range to focus on majority dots in visualization
plt.tight_layout()
plt.show()
plt.close()


## MED1 ChIP signal in H2AZ site in loops
H2AZ_bind_loop_anchor = pybedtools.BedTool.from_dataframe(
    pd.DataFrame(H2AZ_bind_loops.iloc[:,:3].values.tolist()+H2AZ_bind_loops.iloc[:,3:6].values.tolist())
)

med1_peaks_h2az_occu_loop = med1_peaks_plusCL_shNT_bed.intersect(H2AZ_bind_loop_anchor, wa = True).to_dataframe().drop_duplicates()
med1_peaks_h2az_occu_loop_center = med1_peaks_h2az_occu_loop.copy()
c = med1_peaks_h2az_occu_loop_center['start']+(med1_peaks_h2az_occu_loop_center['end']-med1_peaks_h2az_occu_loop_center['start']).astype('int')
med1_peaks_h2az_occu_loop_center['start'] = c
med1_peaks_h2az_occu_loop_center['end'] = c+1

med1_peaks_h2az_occu_loop.iloc[:,:3].to_csv('med1_peaks_h2az_occu_loop.bed', sep = '\t', header = None, index = None)
## to generate heatmap using deepTools +- 600bp around peak center, the result in med1_peaks_h2az_occu_loop_MED1_bw.mat file
bw_mat = pd.read_csv('./med1_peaks_h2az_occu_loop_MED1_bw.mat', sep = '\t', header = None)
minusCL_r1 = bw_mat.iloc[:,6:(6+600)]
minusCL_r2 = bw_mat.iloc[:,(6+600):(6+1200)]
plusCL_r1 = bw_mat.iloc[:,(6+1200):(6+1800)]
plusCL_r2 = bw_mat.iloc[:,(6+1800):(6+2400)]
KD_r1 = bw_mat.iloc[:,(6+2400):(6+3000)]
KD_r2 = bw_mat.iloc[:,(6+3000):]

minusCL_avg = pd.DataFrame((np.array(minusCL_r1) + np.array(minusCL_r2))/2)
plusCL_avg = pd.DataFrame((np.array(plusCL_r1) + np.array(plusCL_r2))/2)
## violinplot, 1kb window of MED1 peak center
plot_df = pd.DataFrame([
    minusCL_avg.iloc[:, 250:350].T.mean().tolist()+plusCL_avg.iloc[:, 250:350].T.mean().tolist()+KD_avg.iloc[:, 250:350].T.mean().tolist(),
    ['minusCL']*minusCL_avg.shape[0]+['plusCL']*plusCL_avg.shape[0]+['KD']*KD_avg.shape[0]
], index = ['value', 'cond']).T

fig, ax = plt.subplots(figsize = (2,2.5))
sns. boxplot(data = plot_df.query('cond != "KD"'), x = 'cond', y = 'value', showfliers = False, 
               palette={'minusCL': 'grey', 'plusCL': 'darkorange', 'KD': plt.cm.get_cmap('tab10')(0)}, whis = 0.8)
# sns.swarmplot(data = plot_df.query('cond != "KD"'), x = 'cond', y = 'value', size = 1,
#                palette={'minusCL': 'grey', 'plusCL': 'darkorange', 'KD': plt.cm.get_cmap('tab10')(0)})
ax.set_ylabel('MED1 ChIP-seq Signal', fontsize = 8)
ax.set_xlabel('')
sns.despine()
plt.tight_layout()
plt.show()
plt.close()

## plot correlation between replicate computed from hicrep software
pdf = PdfPages('mouse_microc_hicrep_replicate_heatmap.pdf')

## shNT plusCL
folder = './MicroC/replicate_reproduce/shNT_plusCL'
files = [x for x in os.listdir(folder) if x.startswith('rep') and x.endswith('.txt')]
corr_res = []
for path in files:
    if path == 'replicates_hicrep_corr_mat.txt':
        continue
    with open(os.path.join(folder, path)) as f:
        v = [float(x.rstrip().split(' ')[-1]) for x in f.readlines()[2:]]
        corr_res.append([path.split('_')[0], path.split('_')[1], np.median(v)])
df = pd.DataFrame(corr_res).pivot_table(index = 0, columns = 1, values = 2).reindex(index = ['rep1', 'rep2', 'rep3', 'rep4'],
                                                                              columns = ['rep1', 'rep2', 'rep3', 'rep4'])
## make symatrix
for x1 in ['rep1', 'rep2', 'rep3', 'rep4']:
    col = df.loc[:,x1]
    for x2 in col.index:
        if pd.isna(col[x2]):
            if x1 == x2:
                df.loc[x1, x2] = 1
            else:
                df.loc[x2, x1] = df.loc[x1, x2]

## plot heatmap
fig, ax = plt.subplots(figsize = (4.5,3.5))
sns.heatmap(df, cmap = 'Reds', vmin = 0, annot = df, fmt = '.4f')
ax.set(xlabel = '', ylabel = '', title = 'shNT_plusCL')
plt.tight_layout()
pdf.savefig(fig)
plt.close()
fig, ax = plt.subplots(figsize = (4.5,3.5))
sns.heatmap(df, cmap = 'Reds', vmin = 0)
ax.set(xlabel = '', ylabel = '', title = 'shNT_plusCL')
plt.tight_layout()
pdf.savefig(fig)
plt.close()

## shNT plusCL
folder = './MicroC/replicate_reproduce/shNT_minusCL'
files = [x for x in os.listdir(folder) if x.startswith('rep') and x.endswith('.txt')]
corr_res = []
for path in files:
    if path == 'replicates_hicrep_corr_mat.txt':
        continue
    with open(os.path.join(folder, path)) as f:
        v = [float(x.rstrip().split(' ')[-1]) for x in f.readlines()[2:]]
        corr_res.append([path.split('_')[0], path.split('_')[1], np.median(v)])
df = pd.DataFrame(corr_res).pivot_table(index = 0, columns = 1, values = 2).reindex(index = ['rep1', 'rep2', 'rep3', 'rep4'],
                                                                              columns = ['rep1', 'rep2', 'rep3', 'rep4'])
## make symatrix
for x1 in ['rep1', 'rep2', 'rep3', 'rep4']:
    col = df.loc[:,x1]
    for x2 in col.index:
        if pd.isna(col[x2]):
            if x1 == x2:
                df.loc[x1, x2] = 1
            else:
                df.loc[x2, x1] = df.loc[x1, x2]

## plot heatmap
fig, ax = plt.subplots(figsize = (4.5,3.5))
sns.heatmap(df, cmap = 'Reds', vmin = 0, annot = df, fmt = '.4f')
ax.set(xlabel = '', ylabel = '', title = 'shNT_minusCL')
plt.tight_layout()
pdf.savefig(fig)
plt.close()
fig, ax = plt.subplots(figsize = (4.5,3.5))
sns.heatmap(df, cmap = 'Reds', vmin = 0)
ax.set(xlabel = '', ylabel = '', title = 'shNT_minusCL')
plt.tight_layout()
pdf.savefig(fig)
plt.close()

## shNT plusCL
folder = './MicroC/replicate_reproduce/KD_plusCL'
files = [x for x in os.listdir(folder) if x.startswith('rep') and x.endswith('.txt')]
corr_res = []
for path in files:
    if path == 'replicates_hicrep_corr_mat.txt':
        continue
    with open(os.path.join(folder, path)) as f:
        v = [float(x.rstrip().split(' ')[-1]) for x in f.readlines()[2:]]
        corr_res.append([path.split('_')[0], path.split('_')[1], np.median(v)])
df = pd.DataFrame(corr_res).pivot_table(index = 0, columns = 1, values = 2).reindex(index = ['rep1', 'rep2', 'rep3', 'rep4'],
                                                                              columns = ['rep1', 'rep2', 'rep3', 'rep4'])
## make symatrix
for x1 in ['rep1', 'rep2', 'rep3', 'rep4']:
    col = df.loc[:,x1]
    for x2 in col.index:
        if pd.isna(col[x2]):
            if x1 == x2:
                df.loc[x1, x2] = 1
            else:
                df.loc[x2, x1] = df.loc[x1, x2]

## plot heatmap
fig, ax = plt.subplots(figsize = (4.5,3.5))
sns.heatmap(df, cmap = 'Reds', vmin = 0, annot = df, fmt = '.4f')
ax.set(xlabel = '', ylabel = '', title = 'KD_plusCL')
plt.tight_layout()
pdf.savefig(fig)
plt.close()
fig, ax = plt.subplots(figsize = (4.5,3.5))
sns.heatmap(df, cmap = 'Reds', vmin = 0)
ax.set(xlabel = '', ylabel = '', title = 'KD_plusCL')
plt.tight_layout()
pdf.savefig(fig)
plt.close()
pdf.close()


## diff loop with gene expression change in RNA-seq

## identify H2AZ-dependant loops - CL induce and KD down loops - 2fold change
tmp = pd.DataFrame(quantile_norm(loop_dots[['score_shNT_plusCL', 'score_shNT_minusCL', 'score_KD_plusCL']].copy()), 
                   columns = ['score_shNT_plusCL', 'score_shNT_minusCL', 'score_KD_plusCL'], index = loop_dots['label'].tolist())

CL_induce_KD_down = up_dots[up_dots['label'].isin(KD_down_dots['label']) & up_dots['label'].isin(up_loop['label'])]
CL_induce_KD_down['score_KD_plusCL'] = tmp.loc[CL_induce_KD_down['label'].tolist(), 'score_KD_plusCL'].tolist()
CL_induce_KD_down['score_shNT_plusCL'] = tmp.loc[CL_induce_KD_down['label'].tolist(), 'score_shNT_plusCL'].tolist()
CL_induce_KD_down['score_shNT_minusCL'] = tmp.loc[CL_induce_KD_down['label'].tolist(), 'score_shNT_minusCL'].tolist()


CL_induce_KD_down_EP, CL_induce_KD_down_PP, CL_induce_KD_down_PO, CL_induce_KD_down_OO = _define_EP_PP_(CL_induce_KD_down, comp = 'plusCL_KD_vs_shNT')

CL_induce_KD_down_tmp = pd.concat([CL_induce_KD_down_EP, CL_induce_KD_down_PP, CL_induce_KD_down_PO, CL_induce_KD_down_OO])
CL_induce_KD_down_tmp['group'] = ['E-P']*CL_induce_KD_down_EP.shape[0]+['P-P']*CL_induce_KD_down_PP.shape[0]+['P-O']*CL_induce_KD_down_PO.shape[0]+['O-O']*CL_induce_KD_down_OO.shape[0]

CL_induce_KD_down_tmp['score_delta'] = CL_induce_KD_down_tmp['score_KD_plusCL'] - CL_induce_KD_down_tmp['score_shNT_plusCL']
CL_induce_KD_down_tmp['mRNA_fc'] = CL_induce_KD_down_tmp[['r1_mRNA_fc', 'r2_mRNA_fc']].apply(lambda row: row.dropna().mean(), axis = 1)

CL_induce_KD_down_exp = _get_values_(CL_induce_KD_down_tmp, comp = 'plusCL_KD_vs_shNT')
CL_induce_KD_down_z = _get_values_(CL_induce_KD_down_tmp, comp = 'plusCL_KD_vs_shNT', value = 'stat')

df = []
for i, line in CL_induce_KD_down_z.iterrows():
    if not pd.isna(line['r1_gene']):
        for g in line['r1_gene'].split(';'):
            df.append([g, line['group'], line['r1_gene'], line['r2_gene'], line['label']])
    if not pd.isna(line['r2_gene']):
        for g in line['r2_gene'].split(';'):
            df.append([g, line['group'], line['r1_gene'], line['r2_gene'], line['label']])
df = pd.DataFrame(df, columns = ['gene', 'group', 'r1_gene', 'r2_gene', 'loop']).drop_duplicates()
df = df[df['gene'].isin(diffexp['plusCL_KD_vs_shNT'].index)]
df['stat'] = [diffexp['plusCL_KD_vs_shNT'].loc[x, 'stat'] for x in df['gene'].tolist()]

l1 = list(set(df.query('group == "E-P"')['gene'].tolist()))
l2 = list(set(df.query('group == "P-P"')['gene'].tolist()))
l3 = list(set(df.query('group == "P-O"')['gene'].tolist()))

l = list(np.intersect1d(l1, l2))+list(np.intersect1d(l1, l3))+list(np.intersect1d(l2, l3))

fig, ax = plt.subplots(figsize = (4.5,4.5))
sns.boxplot(data = df[~df['r1_gene'].isin(l) & ~df['r2_gene'].isin(l)], x = 'group', y = 'stat', 
            order = ['E-P', 'P-P', 'P-O'], whis = 1.5, showfliers = False)
ax.hlines([0], *ax.get_xlim(), linestyle = 'dashed', color = 'black', linewidth = .5)
ax.set(xlabel = '', ylabel = 'Differential mRNA expression\n(z score of H2AZ KD vs WT)')
ax.set_ylim(-3.3, 3.3)
plt.tight_layout()
# fig.savefig('Figures/CL_induce_H2AZKD_down_loop_target_gene_diff_box.pdf')
plt.show()


## H2AZ signal in CUT&RUN of in vivo adipocytes
H2AZ_invivo_bw_list2={
    "cold_rep1": './H2AZ_bw/T_H_1_mm10.uniq.bam.bw', ## T for cold treatment
    "cold_rep2": './H2AZ_bw/T_H_2_mm10.uniq.bam.bw',
    'cold_rep3': './H2AZ_bw/T_H_2_mm10.uniq.bam.bw',
    'warm_rep1': './H2AZ_bw/C_H_2_mm10.uniq.bam.bw', # C for control 
    'warm_rep2': './H2AZ_bw/C_H_2_mm10.uniq.bam.bw',
    'warm_rep3': './H2AZ_bw/C_H_3_mm10.uniq.bam.bw',
}

H2AZ_invivo_plusCL_peaks = pd.read_csv('./H2AZ_bw/T_H_dedup.bgsub.Fnor.peaks.xls', sep = '\t').query('width_above_cutoff > 147')
H2AZ_invivo_plusCL_peaks_bed = pybedtools.BedTool.from_dataframe(H2AZ_invivo_plusCL_peaks)
H2AZ_invivo_minusCL_peaks = pd.read_csv('./H2AZ_bw/C_H_dedup.bgsub.Fnor.peaks.xls', sep = '\t').query('width_above_cutoff > 147')
H2AZ_invivo_minusCL_peaks_bed = pybedtools.BedTool.from_dataframe(H2AZ_invivo_minusCL_peaks)
H2AZ_invivo_plus_minusCL_union_bed = pybedtools.BedTool.from_dataframe(pd.concat([H2AZ_invivo_plusCL_peaks, H2AZ_invivo_minusCL_peaks]).query('abs(start - end) > 147').drop_duplicates())

H2AZocc_loop_anchor1_bed = pybedtools.BedTool.from_dataframe(CL_induce_KD_down[['chrom1', 'start1', 'end1']])
H2AZocc_loop_anchor2_bed = pybedtools.BedTool.from_dataframe(CL_induce_KD_down[['chrom2', 'start2', 'end2']])

H2AZ_invivo_overlapped_H2AZocc_anchor1 = H2AZ_invivo_plus_minusCL_union_bed.intersect(H2AZocc_loop_anchor1_bed, wa = True).to_dataframe().drop_duplicates()
H2AZ_invivo_overlapped_H2AZocc_anchor2 = H2AZ_invivo_plus_minusCL_union_bed.intersect(H2AZocc_loop_anchor2_bed, wa = True).to_dataframe().drop_duplicates()
H2AZ_invivo_H2AZocc_overlapped = pd.concat([H2AZ_invivo_overlapped_H2AZocc_anchor1, H2AZ_invivo_overlapped_H2AZocc_anchor2]).iloc[:,:3].drop_duplicates()


## extract bw to plot
d = 2000
bw_H2AZ_invivo_site_values = []
for path in H2AZ_invivo_bw_list2:
    bw = pyBigWig.open(H2AZ_invivo_bw_list2[path])
    for i, line in H2AZ_invivo_H2AZocc_overlapped.iterrows():
        c = line['start'] + int(abs((line['start']-line['end'])/2))
        line['bw'] = bw.values(line['chrom'], c-d, c+d)
        bw_H2AZ_invivo_site_values.append(line.values.tolist()+[path])

df = pd.DataFrame(bw_H2AZ_invivo_site_values).groupby(4).apply(lambda d: pd.DataFrame(d[3].tolist()).mean())
df['cond'] = [x.replace('_rep1', '').replace('_rep2', '').replace('_rep3', '') for x in df.index.tolist()]
df = df.groupby('cond').mean()

df22 = df.loc[:,range(1000, 3000)]

fig,ax = plt.subplots(figsize = (6, 4))
ax.plot(range(0, df.shape[1]), df.loc['cold',:], label = 'Cold', linewidth = 1.5, color = plt.cm.get_cmap('tab10')(1))
ax.plot(range(0, df.shape[1]), df.loc['warm',:], label = 'TH', linewidth = 1.5, color = 'grey')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon = False)
ax.set(xlabel = 'H2AZ peaks in H2AZ-occu loop', ylabel = 'H2AZ signal')
sns.despine()
# ax.set_ylim(0.6, 2.3)
plt.tight_layout()
fig.savefig('H2AZ_invivo_signal_H2AZoccu_loop.pdf')
plt.show()
plt.close()
