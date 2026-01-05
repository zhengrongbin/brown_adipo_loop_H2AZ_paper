#!/usr/bin/env python
# coding: utf-8

# In[3]:


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



# In[2]:


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
def _plot_loop_exp_(mm10_tss_ann, up_loops, down_loops, stable_loops, de, d = 10000, figsize = (5,5), pdf = False):
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
    col = 'log2FoldChange'
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
             dpi=300, figsize = (5, 5)):
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
        t1 = ax.text(center - aradius, center, str((a - b).count()), **kwargs)

        # Unique to B
        t2 = ax.text(center + bradius, center, str((b - a).count()), **kwargs)

        t3 = ax.text(
                center, center, str((a + b).count()), **kwargs
            )
        adjust_text([t1, t2, t3], #arrowprops=dict(arrowstyle="-", lw=0.5), 
                    save_steps = False, **kwargs)
    else:
        print([str((a - b).count()), str((a + b).count()), str((b - a).count())])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 8)
    plt.tight_layout()
    fig.savefig(outfn, dpi=dpi) if outfn else None
    plt.show()
    plt.close(fig)



# ### Micro-C reproducibility

# In[ ]:


### Use hicrep software to meansure the correlation between replication
## script for running hicrep could be found at replicateCorrelation.sh

pdf = PdfPages('human_microc_hicrep_replicate_heatmap.pdf')

## shNT plusCL
folder = './replicate_correlation/'
files = [x for x in os.listdir(folder) if x.startswith('shNT_plusFSK_rep') and x.endswith('.txt')]
corr_res = []
for path in files:
    if path == 'replicates_hicrep_corr_mat.txt':
        continue
    with open(os.path.join(folder, path)) as f:
        v = [float(x.rstrip().split(' ')[-1]) for x in f.readlines()[2:]]
        corr_res.append([path.split('_')[2], path.split('_')[3], np.median(v)])
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
ax.set(xlabel = '', ylabel = '', title = 'shNT_plusFSK')
plt.tight_layout()
pdf.savefig(fig)
plt.close()
fig, ax = plt.subplots(figsize = (4.5,3.5))
sns.heatmap(df, cmap = 'Reds', vmin = 0)
ax.set(xlabel = '', ylabel = '', title = 'shNT_plusFSK')
plt.tight_layout()
pdf.savefig(fig)
plt.close()

## shNT plusCL
folder = './replicate_correlation/'
files = [x for x in os.listdir(folder) if x.startswith('shNT_minusFSK_rep') and x.endswith('.txt')]
corr_res = []
for path in files:
    if path == 'replicates_hicrep_corr_mat.txt':
        continue
    with open(os.path.join(folder, path)) as f:
        v = [float(x.rstrip().split(' ')[-1]) for x in f.readlines()[2:]]
        corr_res.append([path.split('_')[2], path.split('_')[3], np.median(v)])
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
ax.set(xlabel = '', ylabel = '', title = 'shNT_minusFSK')
plt.tight_layout()
pdf.savefig(fig)
plt.close()
fig, ax = plt.subplots(figsize = (4.5,3.5))
sns.heatmap(df, cmap = 'Reds', vmin = 0)
ax.set(xlabel = '', ylabel = '', title = 'shNT_minusFSK')
plt.tight_layout()
pdf.savefig(fig)
plt.close()

## shNT plusCL
folder = './replicate_correlation/'
files = [x for x in os.listdir(folder) if x.startswith('KD_plusFSK_rep') and x.endswith('.txt')]
corr_res = []
for path in files:
    if path == 'replicates_hicrep_corr_mat.txt':
        continue
    with open(os.path.join(folder, path)) as f:
        v = [float(x.rstrip().split(' ')[-1]) for x in f.readlines()[2:]]
        corr_res.append([path.split('_')[2], path.split('_')[3], np.median(v)])
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
ax.set(xlabel = '', ylabel = '', title = 'KD_plusFSK')
plt.tight_layout()
pdf.savefig(fig)
plt.close()
fig, ax = plt.subplots(figsize = (4.5,3.5))
sns.heatmap(df, cmap = 'Reds', vmin = 0)
ax.set(xlabel = '', ylabel = '', title = 'KD_plusFSK')
plt.tight_layout()
pdf.savefig(fig)
plt.close()
pdf.close()


# #### loop analysis

# In[7]:


## TSS annotation
hg38_tss_ann = pd.read_csv('/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/hg38/gencode.v38.annotation.protein_coding.tss.csv',
                          sep = '\t', header = None)
hg38_tss_ann.columns = ['chr', 'start', 'end', 'name', 'label', 'strand']



# In[ ]:


## load loop probability for each dot in pooled data
## the loop probability was computed using Peakachu geomoe_score function followed by pool function
shNT_plusFSK_all = collections.defaultdict()
with open('../../DataProcess/Human_MicroC/merge_pairs_mapq5/combine/loops/NT_plusFSK_merge.pairs_5000.bedpe') as source:
    for line in source:
        p = line.rstrip().split('\t')
        key = (p[0], int(p[1]), int(p[2]), p[3], int(p[4]), int(p[5]))
        scores = list(map(float, p[-2:])) 
        shNT_plusFSK_all[key] = scores

shNT_minusFSK_all = collections.defaultdict()
with open('../../DataProcess/Human_MicroC/merge_pairs_mapq5/combine/loops/NT_minusFSK_merge.pairs_5000.bedpe') as source:
    for line in source:
        p = line.rstrip().split('\t')
        key = (p[0], int(p[1]), int(p[2]), p[3], int(p[4]), int(p[5]))
        scores = list(map(float, p[-2:])) 
        shNT_minusFSK_all[key] = scores
        
KD_plusFSK_all = collections.defaultdict()
with open('../../DataProcess/Human_MicroC/merge_pairs_mapq5/combine/loops/KD_plusFSK_merge.pairs_5000.bedpe') as source:
    for line in source:
        p = line.rstrip().split('\t')
        key = (p[0], int(p[1]), int(p[2]), p[3], int(p[4]), int(p[5]))
        scores = list(map(float, p[-2:])) 
        KD_plusFSK_all[key] = scores
        
### identify confident dots by restricting loop probability > 0.95
prob_cut = 0.95
shNT_plusFSK_good = [x for x in shNT_plusFSK_all if shNT_plusFSK_all[x][0] > prob_cut]
shNT_minusFSK_good = [x for x in shNT_minusFSK_all if shNT_minusFSK_all[x][0] > prob_cut]

shNT_union_good = collections.defaultdict()
for x in shNT_plusFSK_good + shNT_minusFSK_good:
    if x in shNT_union_good:
        continue
    v1 = shNT_plusFSK_all[x] if x in shNT_plusFSK_all else [0, 0]
    v2 = shNT_minusFSK_all[x] if x in shNT_minusFSK_all else [0, 0]
    v3 = KD_plusFSK_all[x] if x in KD_plusFSK_all else [0, 0]
    shNT_union_good[x] = v1+v2+v3

## generate a dataframe, include all the information, all scores in pooled and replicated data
loop_dots = []
for x in list(shNT_union_good.keys()):
    loop_dots.append(list(x)+shNT_union_good[x])
loop_dots = pd.DataFrame(loop_dots, columns = ['chrom1', 'start1', 'end1',
                                              'chrom2', 'start2', 'end2',
                                              'prob_shNT_plusFSK', 'score_shNT_plusFSK', 
                                              'prob_shNT_minusFSK', 'score_shNT_minusFSK', 
                                              'prob_KD_plusFSK', 'score_KD_plusFSK'])

loop_dots.to_csv('human_loop_dots.csv', index = None)


# In[8]:


### ----------- after union all dots, differential dot analysis performed between plusFSK vs minusFSK----------------
### identify differential dots by 2-fold change thresthold
loop_dots = pd.read_csv('human_loop_dots.csv')
cord = ['chrom1','start1','end1','chrom2','start2','end2']

loop_dots['log2fc'] = loop_dots.apply(lambda row: np.log2(row['prob_shNT_plusFSK'])-np.log2(np.clip(row['prob_shNT_minusFSK'], a_min=0.0001, a_max = 1)), axis = 1)
loop_dots['label'] = loop_dots[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)

loop_dots.index = loop_dots['label'].tolist()
## define up and down dots using 2-fold as cutoff
up_dots = loop_dots.query('log2fc > 1')
up_dots['Type'] = 'plusCL specific'
down_dots = loop_dots.query('log2fc < -1')
down_dots['Type'] = 'minusCL specific'
stable_dots = loop_dots.query('log2fc >= -1 and log2fc <= 1')
stable_dots['Type'] = 'No'

## output dots to pool them as loops in Peakachu
cord = ['chrom1','start1','end1','chrom2','start2','end2']
up_dots[cord+['prob_shNT_plusFSK','score_shNT_plusFSK']].to_csv('shNT_plusFSK_specific_dots_5000.bedpe',
                                                             sep = '\t', index = None, header = None)
down_dots[cord+['prob_shNT_minusFSK','score_shNT_minusFSK']].to_csv('shNT_minusFSK_specific_dots_5000.bedpe',
                                                             sep = '\t', index = None, header = None)
stable_dots[cord+['prob_shNT_minusFSK','score_shNT_minusFSK']].to_csv('shNT_plusFSK_minusFSK_stable_dots_5000.bedpe',
                                                             sep = '\t', index = None, header = None)

## differential loop identification from differential dots,
## run peakachu pool to pool dots into loops, e.g.
## peakachu pool -r 5000 -i shNT_plusFSK_minusFSK_stable_dots_5000.bedpe -o shNT_plusFSK_minusFSK_stable_dots_5000.pool.bedpe -t 0.95



# #### downstream analysis for diff loops

# In[ ]:


loop_dots = pd.read_csv('human_loop_dots.csv')
cord = ['chrom1','start1','end1','chrom2','start2','end2']
loop_dots['label'] = loop_dots[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)

### read in the Peakachu pooled loops
up_loop = pd.read_csv('shNT_plusFSK_specific_loops_5000.pool.bedpe',
                      sep = '\t', header = None)
up_loop.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'prob', 'score']
up_loop['size1'] = up_loop['end1'] - up_loop['start1']
up_loop['size2'] = up_loop['end2'] - up_loop['start2']
up_loop['distance'] = up_loop['start2'] - up_loop['start1']

###
down_loop = pd.read_csv('shNT_minusFSK_specific_loops_5000.pool.bedpe',
                       sep = '\t', header = None)
down_loop.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'prob', 'score']
down_loop['size1'] = down_loop['end1'] - down_loop['start1']
down_loop['size2'] = down_loop['end2'] - down_loop['start2']
down_loop['distance'] = down_loop['start2'] - down_loop['start1']
####
stable_loop = pd.read_csv('shNT_plusFSK_minusFSK_stable_loops_5000.pool.bedpe',
                       sep = '\t', header = None)
stable_loop.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'prob', 'score']
stable_loop['size1'] = stable_loop['end1'] - stable_loop['start1']
stable_loop['size2'] = stable_loop['end2'] - stable_loop['start2']
stable_loop['distance'] = stable_loop['start2'] - stable_loop['start1']

up_loop['label'] = up_loop[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)
down_loop['label'] = down_loop[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)
stable_loop['label'] = stable_loop[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)

## union loop set
df1 = loop_dots[loop_dots['label'].isin(up_loop['label'])]
df1['Type'] = 'Up'
df2 = loop_dots[loop_dots['label'].isin(down_loop['label'])]
df2['Type'] = 'Down'
df3 = loop_dots[loop_dots['label'].isin(stable_loop['label'])]
df3['Type'] = 'No'
union_loops = pd.concat([df1, df2, df3])
union_loops['loop_length'] = np.abs(union_loops['start1'] - union_loops['start2'])


# In[ ]:


## quantile-normalization and show the loop scores for diff loops
plot_df = pd.DataFrame(quantile_norm(union_loops[['score_shNT_plusFSK', 'score_shNT_minusFSK']]),
            columns = ['shNT_plusFSK', 'shNT_minusFSK'])
plot_df['Type']=union_loops['Type'].tolist()
pplot_df = plot_df.melt(id_vars = ['Type'])
pplot_df['variable'] = pd.Categorical(pplot_df['variable'], ['shNT_minusFSK', 'shNT_plusFSK'])
pplot_df['Type'] = pd.Categorical(pplot_df['Type'], ['Up', 'Down', 'No'])

## plot
fig, ax = plt.subplots(figsize = (5, 3))
sns.violinplot(data = pplot_df, x = 'Type', y = 'value',
               hue = 'variable', ax = ax, scale = 'width', 
               palette={'shNT_minusFSK':'grey', 
                        'shNT_plusFSK':plt.cm.get_cmap('tab10')(1)})
ax.set(xlabel = '', ylabel = 'Contact Frequency')
ax.legend(title = 'Human MicroC', loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.show()
plt.close()


# In[ ]:


## same analysis for KD vs NT
### KD vs NT
cord = ['chrom1','start1','end1','chrom2','start2','end2']

loop_dots = pd.read_csv('human_loop_dots.csv')
loop_dots_copy = loop_dots.copy()
loop_dots_copy['log2fc'] = loop_dots_copy.apply(lambda row: np.log2(row['prob_KD_plusFSK'])-np.log2(np.clip(row['prob_shNT_plusFSK'], a_min=0.0001, a_max = 1)), axis = 1)
loop_dots_copy['delta'] = loop_dots_copy['prob_KD_plusFSK'] - loop_dots_copy['prob_shNT_plusFSK']
loop_dots_copy['label'] = loop_dots_copy[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)

## define up and down dots using 2-fold as cutoff
KD_up_dots = loop_dots_copy.query('log2fc > 1')
KD_up_dots['Type'] = 'KD_vs_shNT_plusFSK_up'
KD_down_dots = loop_dots_copy.query('log2fc < -1')
KD_down_dots['Type'] = 'KD_vs_shNT_plusFSK_down'
KD_stable_dots = loop_dots_copy.query('log2fc >= -1 and log2fc <= 1')
KD_stable_dots['Type'] = 'KD_vs_shNT_plusFSK_No'

# ## output dots to pool them as loops in Peakachu, prob cut as 0.95
cord = ['chrom1','start1','end1','chrom2','start2','end2']
KD_up_dots[cord+['prob_KD_plusFSK','score_KD_plusFSK']].to_csv('KD_vs_shNT_plusFSK_up_dots_5000.bedpe',
                                                             sep = '\t', index = None, header = None)
KD_down_dots[cord+['prob_shNT_plusFSK','score_shNT_plusFSK']].to_csv('KD_vs_shNT_plusFSK_down_dots_5000.bedpe',
                                                             sep = '\t', index = None, header = None)
KD_stable_dots[cord+['prob_shNT_plusFSK','score_shNT_plusFSK']].to_csv('KD_vs_shNT_plusFSK_stable_dots_5000.bedpe',
                                                             sep = '\t', index = None, header = None)

del loop_dots_copy

## H2AZ KD
### read in the Peakachu pooled loops
KD_up_loop = pd.read_csv('KD_vs_shNT_plusFSK_up_loops_5000.pool.bedpe',
                      sep = '\t', header = None)
KD_up_loop.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'prob', 'score']
KD_up_loop['size1'] = KD_up_loop['end1'] - KD_up_loop['start1']
KD_up_loop['size2'] = KD_up_loop['end2'] - KD_up_loop['start2']
KD_up_loop['distance'] = KD_up_loop['start2'] - KD_up_loop['start1']

###
KD_down_loop = pd.read_csv('KD_vs_shNT_plusFSK_down_loops_5000.pool.bedpe',
                       sep = '\t', header = None)
KD_down_loop.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'prob', 'score']
KD_down_loop['size1'] = KD_down_loop['end1'] - KD_down_loop['start1']
KD_down_loop['size2'] = KD_down_loop['end2'] - KD_down_loop['start2']
KD_down_loop['distance'] =KD_down_loop['start2'] - KD_down_loop['start1']
####
KD_stable_loop = pd.read_csv('KD_vs_shNT_plusFSK_stable_loops_5000.pool.bedpe',
                       sep = '\t', header = None)
KD_stable_loop.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'prob', 'score']
KD_stable_loop['size1'] = KD_stable_loop['end1'] - KD_stable_loop['start1']
KD_stable_loop['size2'] = KD_stable_loop['end2'] - KD_stable_loop['start2']
KD_stable_loop['distance'] = KD_stable_loop['start2'] - KD_stable_loop['start1']

KD_up_loop['label'] = KD_up_loop[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)
KD_down_loop['label'] = KD_down_loop[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)
KD_stable_loop['label'] = KD_stable_loop[cord].apply(lambda x: '-'.join(x.astype('str').tolist()), axis = 1)


### barplot showing the number of differential loops upon H2AZ KD
print({'Up':KD_up_loop.shape[0], 'Down':KD_down_loop.shape[0], 'Stable':KD_stable_loop.shape[0]})
fig, ax = plt.subplots(figsize = (5,5))
ax.bar(['Up', 'Down', 'No'], 
       [KD_up_loop.shape[0], KD_down_loop.shape[0], KD_stable_loop.shape[0]],
      color = 'grey')
ax.set(xlabel='Loop category\n(KD_plusFSK vs shNT_plusFSK)',
      ylabel = 'Number of loops', title = 'Human MicroC')
sns.despine()
plt.tight_layout()
plt.show()
plt.close()

### violin plot to show the loop score change for FSK up loops in control, FSK treatment, and H2AZKD
## quantile-normalization
plot_df = pd.DataFrame(quantile_norm(union_loops[['score_shNT_plusFSK', 'score_shNT_minusFSK', 'score_KD_plusFSK']]),
            columns = ['shNT_plusFSK', 'shNT_minusFSK', 'KD_plusFSK'])
plot_df['Type'] = union_loops['Type'].tolist()
pplot_df = plot_df.melt(id_vars = ['Type'])
pplot_df['variable'] = pd.Categorical(pplot_df['variable'], ['shNT_minusFSK', 'shNT_plusFSK', 'KD_plusFSK'])
pplot_df['Type'] = pd.Categorical(pplot_df['Type'], ['Up', 'Down', 'No'])

fig, ax = plt.subplots(figsize = (4, 4))
sns.violinplot(data = pplot_df.query('Type == "Up"'), x = 'variable', y = 'value',
            ax = ax, scale = 'width', cut = 1,
               palette={'shNT_minusFSK':'grey',  
                        'shNT_plusFSK':plt.cm.get_cmap('tab10')(1),
                       'KD_plusFSK':"#1E90FF"})
ax.set(xlabel = '', ylabel = 'Contact Frequency')
ax.tick_params(axis = 'x', rotation = 90)
plt.tight_layout()
plt.show()
plt.close()


# ### APA plot to showing MicroC signal for up, down, and stable loop

# In[ ]:


r = 10000
cond1 = 'NT_plusFSK'
cond2 = 'NT_minusFSK'
cond3 = 'KD_plusFSK'
## load cool file for APA plots
cool1_path = '../../DataProcess/Human_MicroC/merge_pairs_mapq5/combine/cool/'+cond1+'_merge.pairs_%s.cool'%r
cool2_path = '../../DataProcess/Human_MicroC/merge_pairs_mapq5/combine/cool/'+cond2+'_merge.pairs_%s.cool'%r
cool3_path = '../../DataProcess/Human_MicroC/merge_pairs_mapq5/combine/cool/'+cond3+'_merge.pairs_%s.cool'%r

clr1 = cooler.Cooler(cool1_path)
clr2 = cooler.Cooler(cool2_path)
clr3 = cooler.Cooler(cool3_path)

## Use bioframe to fetch the genomic features from the UCSC.
hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
hg38_cens = bioframe.fetch_centromeres('hg38')
hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)

# Select only chromosomes that are present in the cooler. 
# This step is typically not required! we call it only because the test data are reduced. 
hg38_arms = hg38_arms.set_index("chrom").loc[clr1.chromnames].reset_index()
# call this to automaticly assign names to chromosomal arms:
hg38_arms = bioframe.make_viewframe(hg38_arms)

## compute expected loop scores as background for hg38
expected1 = expected_cis(
            clr1,
            ignore_diags=0,
            view_df=hg38_arms,
            chunksize=100000)
out = open('expected1.pk', 'wb')
pk.dump(expected1, out)
out.close()

expected2 = expected_cis(
            clr2,
            ignore_diags=0,
            view_df=hg38_arms,
            chunksize=1000000)
out = open('expected2.pk', 'wb')
pk.dump(expected2, out)
out.close()

expected3 = expected_cis(
            clr3,
            ignore_diags=0,
            view_df=hg38_arms,
            chunksize=1000000)

out = open('expected3.pk', 'wb')
pk.dump(expected3, out)
out.close()

out = open('expected1.pk', 'rb')
expected1 = pk.load(out)
out.close()

out = open('expected2.pk', 'rb')
expected2 = pk.load(out)
out.close()

out = open('expected3.pk', 'rb')
expected3 = pk.load(out)
out.close()

### plot APA for up loop in three conditions
pup1 = coolpup.pileup(clr1, up_loop,
                      features_format='bedpe', view_df=hg38_arms,
                      expected_df = expected1,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)

plotpup.plot(pup1,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, 
             vmax = 3.3, vmin = 1.2,
             height=1.5)
plt.show()

pup2 = coolpup.pileup(clr2, up_loop,
                      features_format='bedpe', view_df=hg38_arms,
                      expected_df = expected2,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup2,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 3.3, vmin = 1.2,
             height=1.5)
plt.show()

pup3 = coolpup.pileup(clr3, up_loop,
                      features_format='bedpe', view_df=hg38_arms,
                      expected_df = expected2,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup3,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 3.3, vmin = 1.2,
             height=1.5)
plt.show()


### plot APA for down loop in three conditions
pup1 = coolpup.pileup(clr1, down_loop,
                      features_format='bedpe', view_df=hg38_arms,
                      expected_df = expected1,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)

plotpup.plot(pup1,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, 
             vmax = 3.3, vmin = 1.2,
             height=1.5)
plt.show()

pup2 = coolpup.pileup(clr2, down_loop,
                      features_format='bedpe', view_df=hg38_arms,
                      expected_df = expected2,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup2,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 3.3, vmin = 1.2,
             height=1.5)
plt.show()

pup3 = coolpup.pileup(clr3, down_loop,
                      features_format='bedpe', view_df=hg38_arms,
                      expected_df = expected2,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup3,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 3.3, vmin = 1.2,
             height=1.5)
plt.show()

### plot APA for stable loop in three conditions
pup1 = coolpup.pileup(clr1, stable_loop,
                      features_format='bedpe', view_df=hg38_arms,
                      expected_df = expected1,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)

plotpup.plot(pup1,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, 
             vmax = 4.4, vmin = 1.4,
             height=1.5)
plt.show()

pup2 = coolpup.pileup(clr2, stable_loop,
                      features_format='bedpe', view_df=hg38_arms,
                      expected_df = expected2,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup2,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 4.4, vmin = 1.4,
             height=1.5)
plt.show()

pup3 = coolpup.pileup(clr3, stable_loop,
                      features_format='bedpe', view_df=hg38_arms,
                      expected_df = expected2,
#                       rescale=True, rescale_flank=1, 
                        flank=50000, nproc=2, min_diag = 0)
plotpup.plot(pup3,
             score=True,
             center = 1,
             cmap='Reds',
             sym=False, #scale='log', 
             vmax = 4.4, vmin = 1.4,
             height=1.5)
plt.show()




# ### loop integration with gene expression

# In[ ]:


def _get_values_sig_(df, diff_exp_table, value = 'log2FoldChange', x_range = [], gene_sig = False):
    dx, dy = [], []
    ddf = []
    if gene_sig == True:
        diff_exp_table = diff_exp_table.query('padj < 0.05').copy()
    for i, line in df.iterrows():
        g = line[['r1_gene', 'r2_gene']].dropna().values.tolist()
        for j in g:
            j = j.split(';')
            for jj in j:
                ytmp = diff_exp_table.loc[jj, value] if jj in diff_exp_table.index.tolist() else np.nan
                dy.append(ytmp)
                dx.append(line['score_delta'])
                ddf.append(line.tolist()+[ytmp, jj])

    ddf = pd.DataFrame(ddf, columns = df.columns.tolist()+[value, 'gene'])
#     return(ddf)
    ddf = ddf[~pd.isna(ddf[value]) & ~pd.isna(ddf['score_delta'])]
    ddf = ddf.sort_values('score_delta')
    # both_up_down = np.intersect1d(ddf.query('score_delta > 0')['gene'], ddf.query('score_delta < 0')['gene'])
    # ddf = ddf[~ddf['gene'].isin(both_up_down)]
    if len(x_range) == 2:
        ddf = ddf[(ddf['score_delta'] < x_range[1]) & (ddf['score_delta'] > x_range[0])]
        x, y = ddf['score_delta'], ddf[value]
    else:
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
    


# In[ ]:


## differential gene expression identified using DESeq2
diffexp = {'NT_plusFSK_vs_minusFSK': pd.read_csv('./_plusFSK_vs_minusFSK.protein_gene.csv', index_col = 0),
          'plusFSK_KD_vs_NT': pd.read_csv('./_KD_vs_NT.protein_gene.csv', index_col = 0)}

## set gene promoters for hg38
d = 5000
## promoter define from hg38 genome
hg38_promoter = hg38_tss_ann.copy()
hg38_promoter['start'] = hg38_tss_ann['start'] - d
hg38_promoter['end'] = hg38_tss_ann['end'] + d
hg38_promoter = hg38_promoter[hg38_promoter['start']>0]

hg38_promoter_bed = pybedtools.BedTool.from_dataframe(hg38_promoter)

## define enhancer from human cCRE in ENCODE SCREEN project
hg38_cRE = pd.read_csv('./GRCh38-cCREs.bed', sep = '\t', header = None)
hg38_ELS = hg38_cRE[hg38_cRE[5].str.contains('ELS')]
hg38_ELS_bed = pybedtools.BedTool.from_dataframe(hg38_ELS)
hg38_ELS_bed = hg38_ELS_bed.intersect(hg38_promoter_bed, v = True, wa = True)


# In[ ]:


## after inspect the human H3K27ac ChIP-seq quality, peak fold enrichment greater than 5 was used
k27ac_plusFSK_rep1 = pd.read_csv('/temp_work/ch228298/H2AZ_GEO_2025/ChIPseq_04032025_human_H2AZ_H3K27ac/peaks/peaks/NT_plusFSK_H3K27ac_rep1.rep1_5fold_peaks.narrowPeak', sep = '\t', header = None)
k27ac_plusFSK_rep2 = pd.read_csv('/temp_work/ch228298/H2AZ_GEO_2025/ChIPseq_04032025_human_H2AZ_H3K27ac/peaks/peaks/NT_plusFSK_H3K27ac_rep2.rep1_5fold_peaks.narrowPeak', sep = '\t', header = None)
k27ac_minusFSK_rep1 = pd.read_csv('/temp_work/ch228298/H2AZ_GEO_2025/ChIPseq_04032025_human_H2AZ_H3K27ac/peaks/peaks/NT_minusFSK_H3K27ac_rep1.rep1_5fold_peaks.narrowPeak', sep = '\t', header = None)
k27ac_minusFSK_rep2 = pd.read_csv('/temp_work/ch228298/H2AZ_GEO_2025/ChIPseq_04032025_human_H2AZ_H3K27ac/peaks/peaks/NT_minusFSK_H3K27ac_rep2.rep1_5fold_peaks.narrowPeak', sep = '\t', header = None)

## overlapped peaks by two replicates were used
k27ac_plusFSK = pybedtools.BedTool.from_dataframe(k27ac_plusFSK_rep1.iloc[:,0:3]).intersect(pybedtools.BedTool.from_dataframe(k27ac_plusFSK_rep2.iloc[:,0:3]),
                                                                           wa = True).to_dataframe().drop_duplicates()
k27ac_minusFSK = pybedtools.BedTool.from_dataframe(k27ac_minusFSK_rep1.iloc[:,0:3]).intersect(pybedtools.BedTool.from_dataframe(k27ac_minusFSK_rep2.iloc[:,0:3]),
                                                                           wa = True).to_dataframe().drop_duplicates()
## union peaks between plus and minus FSK condition
k27ac_union_peaks = pd.concat([k27ac_plusFSK, k27ac_minusFSK])


# In[ ]:


### classify loops into E-P, P-P, P-O, and integrate with gene expression

## EP
q1 = '(r1_promoter == True and r2_ELS == True and r2_K27ac == True) or (r2_promoter == True and r1_ELS == True and r1_K27ac == True)'
## PP
q2 = '(r1_promoter == True and r2_promoter == True)'
## P-other
q3 = '(r1_promoter == True and (r2_ELS != True or r2_K27ac != True) and r2_promoter != True) or (r2_promoter == True and (r1_ELS != True or r1_K27ac != True) and r1_promoter != True)'

## ELS enhancer overlap with our H3K27ac peaks for system specific ELS
ELS_k27ac_union_peaks_bed = hg38_ELS_bed.intersect(k27ac_union_peaks_bed, wa = True)

### annotate 
union_loops_genes = union_loops.copy()

score_qnorm = pd.DataFrame(quantile_norm(union_loops[['score_shNT_plusFSK', 'score_shNT_minusFSK']]),
            columns = ['NT_plusFSK', 'NT_minusFSK'])
union_loops_genes['score_shNT_plusFSK'] = score_qnorm['NT_plusFSK'].tolist()
union_loops_genes['score_shNT_minusFSK'] = score_qnorm['NT_minusFSK'].tolist()

union_loops_genes['r1'] = union_loops_genes['chrom1']+':'+union_loops_genes['start1'].astype('str')+'-'+union_loops_genes['end1'].astype('str')
union_loops_genes['r2'] = union_loops_genes['chrom2']+':'+union_loops_genes['start2'].astype('str')+'-'+union_loops_genes['end2'].astype('str')
## genomic location for two anchors
r1 = pybedtools.BedTool.from_dataframe(union_loops_genes[['chrom1','start1','end1', 'r1']])
r2 = pybedtools.BedTool.from_dataframe(union_loops_genes[['chrom2','start2','end2', 'r2']])
## overlap with promoter using two anchors
r1_promoter = r1.intersect(hg38_promoter_bed, wao = True).to_dataframe().drop_duplicates()
r2_promoter = r2.intersect(hg38_promoter_bed, wao = True).to_dataframe().drop_duplicates()

## overlap with ELS using two anchors
r1_ELS = r1.intersect(hg38_ELS_bed, wao = True).to_dataframe().drop_duplicates()
r2_ELS = r2.intersect(hg38_ELS_bed, wao = True).to_dataframe().drop_duplicates()

## overlap with K27ac using two anchors
r1_K27ac = r1.intersect(k27ac_union_peaks_bed, wao = True).to_dataframe().drop_duplicates()
r2_K27ac = r2.intersect(k27ac_union_peaks_bed, wao = True).to_dataframe().drop_duplicates()

## overlap with ELS_K27ac using two anchors
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

union_loops_genes['r1_mRNA_fc'] = [diffexp['NT_plusFSK_vs_minusFSK'].reindex(x.split(';'))['log2FoldChange'].mean() if not pd.isna(x) else np.nan for x in union_loops_genes['r1_gene'].tolist()]
union_loops_genes['r2_mRNA_fc'] = [diffexp['NT_plusFSK_vs_minusFSK'].reindex(x.split(';'))['log2FoldChange'].mean() if not pd.isna(x) else np.nan for x in union_loops_genes['r2_gene'].tolist()]

## define EP, PP, and PO

EP_df = union_loops_genes.query(q1)

PP_df = union_loops_genes.query(q2)

PO_df = union_loops_genes.query(q3)
OO_df = union_loops_genes[~union_loops_genes['label'].isin(EP_df['label'].tolist()+PP_df['label'].tolist()+PO_df['label'].tolist())]

## calculate diff scores
EP_df['score_delta'] = EP_df['score_shNT_plusFSK'] - EP_df['score_shNT_minusFSK']
EP_df['mRNA_fc'] = EP_df[['r1_mRNA_fc', 'r2_mRNA_fc']].apply(lambda row: row.dropna().mean(), axis = 1)

PP_df['score_delta'] = PP_df['score_shNT_plusFSK'] - PP_df['score_shNT_minusFSK']
PP_df['mRNA_fc'] = PP_df[['r1_mRNA_fc', 'r2_mRNA_fc']].apply(lambda row: row.dropna().mean(), axis = 1)

PO_df['score_delta'] = PO_df['score_shNT_plusFSK'] - PO_df['score_shNT_minusFSK']
PO_df['mRNA_fc'] = PO_df[['r1_mRNA_fc', 'r2_mRNA_fc']].apply(lambda row: row.dropna().mean(), axis = 1)



# In[ ]:


# ### unique gene and unique loop for each row for a table
EP_ddf = _get_values_sig_(EP_df, diff_exp_table=diffexp['NT_plusFSK_vs_minusFSK'], value = 'log2FoldChange', x_range = [], gene_sig = False) #.query('padj < 0.05 and abs(log2FoldChange) > 0.5'))
EP_ddf['loop_type'] = 'EP'

PP_ddf = _get_values_sig_(PP_df, diff_exp_table = diffexp['NT_plusFSK_vs_minusFSK'], value = 'log2FoldChange', x_range = [], gene_sig = False) #.query('padj < 0.05 and abs(log2FoldChange) > 0.5'))
PP_ddf['loop_type'] = 'PP'

PO_ddf = _get_values_sig_(PO_df, diff_exp_table = diffexp['NT_plusFSK_vs_minusFSK'], value = 'log2FoldChange', x_range = [], gene_sig = False) #.query('padj < 0.05 and abs(log2FoldChange) > 0.5'))
PO_ddf['loop_type'] = 'PO'

## pool into a dataframe for supp table and perform GO analysis
dat = pd.concat([EP_ddf, PP_ddf, PO_ddf])
dat['score_delta_plusFSK_vs_minusFSK'] = dat['score_shNT_plusFSK'] - dat['score_shNT_minusFSK']
dat['score_delta_KD_vs_NT_plusFSK'] = dat['score_KD_plusFSK'] - dat['score_shNT_plusFSK']
dat['log2FoldChange_plusFSK_vs_minusFSK'] = diffexp['NT_plusFSK_vs_minusFSK'].loc[dat['gene'].tolist(), 'log2FoldChange'].tolist()
dat['padj_plusFSK_vs_minusFSK'] = diffexp['NT_plusFSK_vs_minusFSK'].loc[dat['gene'].tolist(), 'padj'].tolist()
dat['log2FoldChange_KD_vs_NT_plusFSK'] = [diffexp['plusFSK_KD_vs_NT'].loc[x,]['log2FoldChange'] if x in diffexp['plusFSK_KD_vs_NT'].index.tolist() else None for x in dat['gene'].tolist()]
dat['padj_KD_vs_NT_plusFSK'] = [diffexp['plusFSK_KD_vs_NT'].loc[x,]['padj'] if x in diffexp['plusFSK_KD_vs_NT'].index.tolist() else None for x in dat['gene'].tolist()]

dat = dat[['chrom1', 'start1', 'end1','chrom2', 'start2', 'end2',
    'prob_shNT_plusFSK', 'prob_shNT_minusFSK', 'prob_KD_plusFSK', 'Type', 'loop_type',
    'score_delta_plusFSK_vs_minusFSK', 'score_delta_KD_vs_NT_plusFSK', 'gene', 
  'log2FoldChange_plusFSK_vs_minusFSK', 'padj_plusFSK_vs_minusFSK',
 'log2FoldChange_KD_vs_NT_plusFSK', 'padj_KD_vs_NT_plusFSK']]

dat.to_csv('human_loop_prob0.95_delta_gene_table.tsv', sep = '\t', index = None)



# In[ ]:


## showing association bwtween loop change and gene expression change using high confident loops (prob>0.97)
# ### unique gene and unique loop for visulization 
EP_ddf = _get_values_sig_(EP_df.query('prob_shNT_plusFSK > 0.97 or prob_shNT_minusFSK > 0.97'), diff_exp_table=diffexp['NT_plusFSK_vs_minusFSK'], value = 'log2FoldChange', x_range = [], gene_sig = False) #.query('padj < 0.05 and abs(log2FoldChange) > 0.5'))
EP_ddf['loop_type'] = 'EP'

PP_ddf = _get_values_sig_(PP_df.query('prob_shNT_plusFSK > 0.97 or prob_shNT_minusFSK > 0.97'), diff_exp_table = diffexp['NT_plusFSK_vs_minusFSK'], value = 'log2FoldChange', x_range = [], gene_sig = False) #.query('padj < 0.05 and abs(log2FoldChange) > 0.5'))
PP_ddf['loop_type'] = 'PP'

PO_ddf = _get_values_sig_(PO_df.query('prob_shNT_plusFSK > 0.97 or prob_shNT_minusFSK > 0.97'), diff_exp_table = diffexp['NT_plusFSK_vs_minusFSK'], value = 'log2FoldChange', x_range = [], gene_sig = False) #.query('padj < 0.05 and abs(log2FoldChange) > 0.5'))
PO_ddf['loop_type'] = 'PO'


### plot loop delta contact vs. gene expression change
fig, ax = plt.subplots(figsize = (6.5, 5))

## align different loop type range to the same, so that can be better compared
xmax = min([max(EP_ddf['score_delta']), max(PP_ddf['score_delta']), max(PO_ddf['score_delta'])])
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
plt.show()



# ### integration with H2AZ KD RNA-seq and H2AZ ChIP-seq in human

# In[ ]:


## H2AZ dependent loop: FSK upregulated but H2AZKD downregulated
FSK_induce_KD_down_loop = up_loop[up_loop['label'].isin(KD_down_dots['label'])]
FSK_induce_KD_down_EP, FSK_induce_KD_down_PP, FSK_induce_KD_down_PO, FSK_induce_KD_down_OO = _define_EP_PP_(FSK_induce_KD_down_loop, comp = 'plusFSK_KD_vs_NT')

## H2AZ ChIP-seq path
H2AZ_path = {'ChIP_H2AZ_plusFSK_rep1':'/temp_work/ch228298/Coll_Yang/dataprocess/ChIPseq/H2AZ_04032025_danpos/NT_plusFSK_rep1_vs_input/pooled/bam_NT_plusFSK_H2AZ_rep1_unique.sorted.bam.chromSize.bgsub.Fnor.peaks.xls',
            'ChIP_H2AZ_plusFSK_rep2':'/temp_work/ch228298/Coll_Yang/dataprocess/ChIPseq/H2AZ_04032025_danpos/NT_plusFSK_rep2_vs_input/pooled/bam_NT_plusFSK_H2AZ_rep2_unique.sorted.bam.chromSize.bgsub.Fnor.peaks.xls',
            'ChIP_H2AZ_minusFSK_rep1':'/temp_work/ch228298/Coll_Yang/dataprocess/ChIPseq/H2AZ_04032025_danpos/NT_minusFSK_rep1_vs_input/pooled/bam_NT_minusFSK_H2AZ_rep1_unique.sorted.bam.chromSize.bgsub.Fnor.peaks.xls',
            'ChIP_H2AZ_minusFSK_rep2':'/temp_work/ch228298/Coll_Yang/dataprocess/ChIPseq/H2AZ_04032025_danpos/NT_minusFSK_rep2_vs_input/pooled/bam_NT_minusFSK_H2AZ_rep2_unique.sorted.bam.chromSize.bgsub.Fnor.peaks.xls',
            'ChIP_H2AZ_KD_plusFSK_rep1':'/temp_work/ch228298/Coll_Yang/dataprocess/ChIPseq/H2AZ_04032025_danpos/KD_plusFSK_rep1_vs_input/pooled/bam_KD_plusFSK_H2AZ_rep1_unique.sorted.bam.chromSize.bgsub.Fnor.peaks.xls',
            'ChIP_H2AZ_KD_plusFSK_rep2':'/temp_work/ch228298/Coll_Yang/dataprocess/ChIPseq/H2AZ_04032025_danpos/KD_plusFSK_rep2_vs_input/pooled/bam_KD_plusFSK_H2AZ_rep2_unique.sorted.bam.chromSize.bgsub.Fnor.peaks.xls'
            }

H2AZ_peaks = {}
for x in H2AZ_path:
    tmp = pd.read_csv(H2AZ_path[x], sep = '\t').query('width_above_cutoff > 147')
    H2AZ_peaks[x] = tmp

H2AZ_peaks_bed = {}
for x in H2AZ_peaks:
    H2AZ_peaks_bed[x] = pybedtools.BedTool.from_dataframe(H2AZ_peaks[x])
    
# ### highly reproduciable peaks by replicates
H2AZ_plusFSK_reprod_peaks = H2AZ_peaks_bed['ChIP_H2AZ_plusFSK_rep1'].intersect(H2AZ_peaks_bed['ChIP_H2AZ_plusFSK_rep2'], wa=True).to_dataframe().drop_duplicates()
H2AZ_plusFSK_reprod_peaks_bed = pybedtools.BedTool.from_dataframe(H2AZ_plusFSK_reprod_peaks)

H2AZ_minusFSK_reprod_peaks = H2AZ_peaks_bed['ChIP_H2AZ_minusFSK_rep1'].intersect(H2AZ_peaks_bed['ChIP_H2AZ_minusFSK_rep2'], wa=True).to_dataframe().drop_duplicates()
H2AZ_minusFSK_reprod_peaks_bed = pybedtools.BedTool.from_dataframe(H2AZ_minusFSK_reprod_peaks)



# In[ ]:


## ---- percentage of H2AZ occupied loops -----
### identify H2A.Z occupied loops within H2AZ-dependant loops
EP_res = _plot_H2AZ_overlap_v2_(loops=FSK_induce_KD_down_EP, peaks={'H2AZ_peak': H2AZ_plusFSK_reprod_peaks_bed},
                                title = 'FSK_induce_KD_down_EP_H2AZ_overlap',flip = True)
PP_res = _plot_H2AZ_overlap_v2_(loops=FSK_induce_KD_down_PP, peaks={'H2AZ_peak': H2AZ_plusFSK_reprod_peaks_bed},
                                title = 'FSK_induce_KD_down_PP_H2AZ_overlap',flip = True)
PO_res = _plot_H2AZ_overlap_v2_(loops=FSK_induce_KD_down_PO, peaks={'H2AZ_peak': H2AZ_plusFSK_reprod_peaks_bed},
                                title = 'FSK_induce_KD_down_PO_H2AZ_overlap',flip = True)
OO_res = _plot_H2AZ_overlap_v2_(loops=FSK_induce_KD_down_OO, peaks={'H2AZ_peak': H2AZ_plusFSK_reprod_peaks_bed},
                                title = 'FSK_induce_KD_down_OO_H2AZ_overlap',flip = True)
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
plt.show()
plt.close()


# In[ ]:


### ----- gene expression vs loop change in down-regulated genes
KD_union_loop = pd.concat([KD_up_loop, KD_down_loop, KD_stable_loop])
KD_union_loop_EP, KD_union_loop_PP, KD_union_loop_PO, KD_union_loop_OO = _define_EP_PP_(KD_union_loop, comp = 'plusFSK_KD_vs_NT')
KD_union_loop_EP['group'] = 'EP'
KD_union_loop_PP['group'] = 'PP'
KD_union_loop_PO['group'] = 'PO'
KD_union_loop_OO['group'] = 'OO'
KD_union_loop = pd.concat([KD_union_loop_EP, KD_union_loop_PP, KD_union_loop_PO, KD_union_loop_OO])

tmp = pd.DataFrame(quantile_norm(loop_dots[['score_shNT_plusFSK', 'score_shNT_minusFSK', 'score_KD_plusFSK']].copy()), 
                   columns = ['score_shNT_plusFSK', 'score_shNT_minusFSK', 'score_KD_plusFSK'], index = loop_dots['label'].tolist())
KD_union_loop['score_KD_plusFSK'] = tmp.loc[KD_union_loop['label'].tolist(), 'score_KD_plusFSK'].tolist()
KD_union_loop['score_shNT_plusFSK'] = tmp.loc[KD_union_loop['label'].tolist(), 'score_shNT_plusFSK'].tolist()
KD_union_loop['score_delta'] = KD_union_loop['score_KD_plusFSK'] - KD_union_loop['score_shNT_plusFSK']

tmp = pd.DataFrame(quantile_norm(loop_dots[['prob_shNT_plusFSK', 'prob_KD_plusFSK']].copy()), 
                   columns = ['prob_shNT_plusFSK', 'prob_KD_plusFSK'], index = loop_dots['label'].tolist())
KD_union_loop['prob_KD_plusFSK'] = tmp.loc[KD_union_loop['label'].tolist(), 'prob_KD_plusFSK'].tolist()
KD_union_loop['prob_shNT_plusFSK'] = tmp.loc[KD_union_loop['label'].tolist(), 'prob_shNT_plusFSK'].tolist()
KD_union_loop['prob_delta'] = KD_union_loop['prob_KD_plusFSK'] - KD_union_loop['prob_shNT_plusFSK']

### analyze the loop change and differental expression upon H2AZ KD, particularly for H2AZ-dependant loops and genes

KD_union_loop_df = _get_values_sig_(KD_union_loop,
                      diff_exp_table=diffexp['plusFSK_KD_vs_NT'], value = 'stat', x_range = [], gene_sig = False) #.query('padj < 0.05 and abs(log2FoldChange) > 0.5'))

## CL induce and H2AZ KD impaired genes
g1 = diffexp['plusFSK_KD_vs_NT'].query('log2FoldChange < 0')
cut = np.log2(1.5)
g2 = diffexp['NT_plusFSK_vs_minusFSK'].query('log2FoldChange > @cut and padj < 0.05')
g = np.intersect1d(g1.index, g2.index)

## focused on down-regulated genes in this analysis
index1 = [len(np.intersect1d(x.split(';'), g)) > 0 if not pd.isna(x) else False for x in KD_union_loop_df['r1_gene'].tolist()]
index2 = [len(np.intersect1d(x.split(';'), g)) > 0 if not pd.isna(x) else False for x in KD_union_loop_df['r2_gene'].tolist()]
focus_df = KD_union_loop_df.loc[np.array(index1) | np.array(index2),]

focus_df['r1_gene'] = [';'.join(list(np.intersect1d(x.split(';'), g))) if not pd.isna(x) else x for x in focus_df['r1_gene'].tolist()]
focus_df['r2_gene'] = [';'.join(list(np.intersect1d(x.split(';'), g))) if not pd.isna(x) else x for x in focus_df['r2_gene'].tolist()]
focus_df['r1_gene'] = [np.nan if x == '' else x for x in focus_df['r1_gene'].tolist()]
focus_df['r2_gene'] = [np.nan if x == '' else x for x in focus_df['r2_gene'].tolist()]

focus_ddf2 = _get_values_(focus_df.drop(['stat'], axis = 1), diff_exp_table=diffexp['plusFSK_KD_vs_NT'], value = 'stat')
fig, ax = plt.subplots(figsize = (5,5))
sns.scatterplot(data = focus_ddf2.query('group != "OO"'), x = 'stat', y = 'score_delta',
                edgecolor = 'none', s = 5, alpha = .7)
ax.hlines(0, -22, 0, color = 'grey', linestyle = 'dashed')
ax.set(xlabel = 'Differential mRNA expression\n(z score of H2AZ KD vs WT)',
      ylabel = 'Delta contact frequency\n(H2AZ KD - WT)')
ax.set_ylim(-0.015, 0.015)
ax.set_xlim(-22, 0)
plt.tight_layout()
plt.show()
plt.close()



# ### integration with GWAS SNP

# In[ ]:


## consider LD in GWAS snps
hg38_gwas = pd.read_csv('../../MicroC/hg38_GWAScatelog_UCSC.tsv', sep = '\t')

## non-trait associated SNPs for background control
nonsig_snp = hg38_gwas[~hg38_gwas['name'].isin(hg38_gwas.query('pValue < 1e-8')['name'].unique())]
keep = [f"chr{i}" for i in list(range(1,23)) + ["X","Y","M"]]
nonsig_snp = nonsig_snp[nonsig_snp['chrom'].isin(keep)]
nonsig_snp = nonsig_snp[['chrom','chromStart', 'chromEnd']].drop_duplicates()


hg38_gwas = hg38_gwas.query('pValue < 1e-8')## select significant associations
hg38_gwas['label'] = hg38_gwas['chrom']+':'+hg38_gwas['chromStart'].astype('str')+'-'+hg38_gwas['chromEnd'].astype('str')+'_'+hg38_gwas['name']+'_'+hg38_gwas['trait']
hg38_gwas_bed = pybedtools.BedTool.from_dataframe(hg38_gwas[['chrom','chromStart','chromEnd','label']])

print('total gwas snp: %s'%hg38_gwas.shape[0])
## overlap with loop anchors, 100Kb window as it will be used for lead SNP detection in LD
union_loop_anchors = pd.DataFrame(union_loops.iloc[:,:3].values.tolist()+union_loops.iloc[:,3:6].values.tolist()).drop_duplicates()
center = union_loop_anchors[1] + (.5 * (union_loop_anchors[2] - union_loop_anchors[1])).astype('int')
union_loop_anchors[1] = np.clip(center - 50000, a_min = 0, a_max = None)
union_loop_anchors[2] = center + 50000 
union_loop_anchor_bed = pybedtools.BedTool.from_dataframe(union_loop_anchors)
union_loop_anchor_bed = union_loop_anchor_bed.sort()
union_loop_anchors_merge = union_loop_anchor_bed.merge().to_dataframe().drop_duplicates()
union_loop_anchors_merge['name'] = union_loop_anchors_merge['chrom']+':'+union_loop_anchors_merge['start'].astype('str')+'-'+union_loop_anchors_merge['end'].astype('str')


hg38_gwas_loop = hg38_gwas_bed.intersect(pybedtools.BedTool.from_dataframe(union_loop_anchors_merge), wao = True).to_dataframe().drop_duplicates()
hg38_gwas_loop = hg38_gwas_loop.query('score != "."')
print('in up loop anchor gwas snp: %s'%hg38_gwas_loop.shape[0])

# ## pick lead SNP with best significant GWAS in each block and each trait
hg38_gwas_loop = pd.merge(hg38_gwas_loop[['chrom', 'start', 'end','name', 'thickEnd']],
                          hg38_gwas[['label', 'pValue', 'trait']], left_on = 'name', right_on = 'label', how = 'left').drop_duplicates()

hg38_gwas_loop_best = hg38_gwas_loop.groupby(['thickEnd', 'trait']).min().reset_index()
hg38_gwas_loop_best['snp']=[x.split('_')[1] for x in hg38_gwas_loop_best['label'].tolist()]
snp_100kb_windw = hg38_gwas_loop_best['snp'].unique()
out = open('hg38_union_loop_leadSNP_100kwindow.txt', 'w')
out.write('\n'.join(hg38_gwas_loop_best['snp'].unique().tolist()))
out.close()

## --- download LD associations for lead SNPs using LDlink by running the custome script ----
### bash: python download_human_snp_LD.py hg38_union_loop_leadSNP_100kwindow.txt
## --- perform LD expansion by adding LD associated SNPs with lead SNPs -----
### bash: python LD_expand.py
### it will generate hg38_GWAS_LD_expanded.txt file


# In[ ]:


## expanded LD snp sites
hg38_gwas_LD_expended = pd.read_csv('hg38_GWAS_LD_expanded.txt',sep = '\t')
hg38_gwas_LD_expended['trait'] = [x.split('_')[2] for x in hg38_gwas_LD_expended['label'].tolist()]
hg38_gwas_lead = hg38_gwas_loop_best[['chrom', 'start', 'end', 'label', 'snp', 'snp', 'trait']]
hg38_gwas_lead.columns = ['chrom', 'chromStart', 'chromEnd', 'label', 'lead_snp', 'RS_Number',
       'trait']
## prepare H2AZ bind loop anchors for SNP analysis
H2AZ_bind_loops = pd.DataFrame()
for t, r, l in [['EP', EP_res, FSK_induce_KD_down_EP],
         ['PP', PP_res, FSK_induce_KD_down_PP],
         ['PO', PO_res, FSK_induce_KD_down_PO],
         ['OO', OO_res, FSK_induce_KD_down_OO]]:
    tmp = pd.concat([r['H2AZ_peak']['loops']['either'], r['H2AZ_peak']['loops']['both']])
    tmp['anchor1'] = tmp['chrom1']+':'+tmp['start1'].astype('str')+'-'+tmp['end1'].astype('str')
    tmp['anchor2'] = tmp['chrom2']+':'+tmp['start2'].astype('str')+'-'+tmp['end2'].astype('str')
    tmp['Type'] = t
    H2AZ_bind_loops = pd.concat([H2AZ_bind_loops, tmp])

## prepare loop sets for SNP analysis
## H2AZ occupied loop
human_FSK_induce_KD_down_H2AZ_bind_loop_bed = pybedtools.BedTool.from_dataframe(
    pd.DataFrame(H2AZ_bind_loops.iloc[:,:3].values.tolist()+H2AZ_bind_loops.iloc[:,3:6].values.tolist()).drop_duplicates()
)
## background SNP sets
background_nonasso_snp_bed = pybedtools.BedTool.from_dataframe(nonsig_snp[['chrom','chromStart','chromEnd']].drop_duplicates())

from statsmodels import stats as ss
## perform Fisher exact test to evaluate significant overlap between SNP loci and loop anchors
def _fisher_vs_background_(loop_bed, background_bed, gwas_df_bed, gwas_df):
    fishr_res = collections.defaultdict()
    trait_list = loop_bed.intersect(gwas_df_bed, wb = True).to_dataframe()['blockCount'].unique().tolist()
    fishr_res = collections.defaultdict()
    for trait in trait_list:
        trait_bed = pybedtools.BedTool.from_dataframe(gwas_df.query('trait == @trait').iloc[:,0:3])
        # overlap counts
        a = len(trait_bed.intersect(loop_bed, u=True))  # trait in region
        b = len(trait_bed.subtract(loop_bed, A=True))   # trait not in region
        c = len(background_bed.intersect(loop_bed, u=True)) - a
        d = len(background_bed.subtract(loop_bed, A=True)) - b
        N = a + b + c + d
        fold_enrichment = (a * N) / ((a + b) * (a + c))
        oddsratio, p = fisher_exact([[a,b],[c,d]], alternative='greater')
        fishr_res[trait] = [[[a,b],[c,d]], p, oddsratio, fold_enrichment]
    fishr_res = pd.DataFrame(fishr_res, index = ['table', 'pval', 'ratio', 'Fold']).T
    fishr_res['FDR'] = ss.multitest.multipletests(fishr_res['pval'], method='fdr_bh')[1]
    return(fishr_res)

## compute number of overlapped SNPs in GWAS
def _comput_snp_overlap_(loop_bed, gwas_df, promoter_remove = False):
    if promoter_remove:
        print('remove promoter anchors')
        loop_bed = loop_bed.intersect(hg38_promoter_bed, v = True)
    loop_snp = loop_bed.intersect(pybedtools.BedTool.from_dataframe(gwas_df.drop_duplicates()), wao=True).to_dataframe().drop_duplicates()
    loop_snp = loop_snp[loop_snp['name'] != '.']
    # print(loop_snp.head())
    loop_snp = pd.merge(loop_snp, gwas_df[['label', 'trait']].drop_duplicates(), left_on = 'thickStart', right_on = 'label')
    return(loop_snp)
## overlap with SNP
#### H2AZ occupied loop
human_FSK_induce_KD_down_H2AZ_bind_loop_snp = _comput_snp_overlap_(human_FSK_induce_KD_down_H2AZ_bind_loop_bed, 
                                                                   gwas_df = pd.concat([hg38_gwas_LD_expended, hg38_gwas_lead]))

print('% of SNPs in non-promoter anchors: ')
pybedtools.BedTool.from_dataframe(human_FSK_induce_KD_down_H2AZ_bind_loop_snp[['chrom', 'start', 'end', 'itemRgb']]).intersect(
    hg38_promoter_bed, v = True
).to_dataframe().drop_duplicates()['name'].unique().shape[0] * 100 / human_FSK_induce_KD_down_H2AZ_bind_loop_snp['itemRgb'].unique().shape[0]


## perfrom enrichment analysis
## background=non-trait associated snps
gwas_df_bed = pybedtools.BedTool.from_dataframe(gwas_df)
H2AZ_occu_fres_bg = _fisher_vs_background_(human_FSK_induce_KD_down_H2AZ_bind_loop_bed, background_nonasso_snp_bed, gwas_df_bed, gwas_df)


## Visualize obesity related traits using bar chart
traits = ["Waist-to-hip ratio adjusted for BMI",
          "Waist circumference adjusted for body mass index",
          "Metabolically unhealthy in obesity",
"Body fat percentage and type 2 diabetes (pairwise)",
"Body mass index and type 2 diabetes (pairwise)",
"Fasting insulin",
"Fasting glucose",
"Body fat percentage or coronary artery disease (MTAG)",
"Nonalcoholic fatty liver disease (imputed)",
          "Hip circumference adjusted for BMI",
"Weight",
"Triglycerides",
"Coronary artery disease",
"Body mass index",
"Body fat percentage",
"Type 2 diabetes"]

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib import colors

plot_df = H2AZ_occu_fres_bg.loc[traits].reset_index()
plot_df['significance'] = -np.log10(plot_df['FDR'].tolist())
plot_df = plot_df.sort_values('Fold')

# Example setup
values = plot_df['significance']  # e.g., -log10(FDR)
folds = plot_df['Fold']

# Define your custom low and high colors
low_color = "pink"   # light blue (low end)
high_color = "red"  # dark blue (high end)

# Create custom colormap (2-point gradient)
cmap = colors.LinearSegmentedColormap.from_list("custom_cmap", [low_color, high_color])

# Normalize data to color range
norm = colors.Normalize(vmin=0, vmax=14)

# Map normalized values to colors
mapped_colors = cmap(norm(values))

# --- Plot ---
fig, ax = plt.subplots(figsize=(7.5, 4))

bars = ax.barh(
    y=plot_df['index'],
    width=folds,
    color=mapped_colors,
    edgecolor='none'
)

# Add colorbar
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, pad=0.02, shrink = 0.5)
cbar.set_label('-log10(FDR)', rotation=270, labelpad=15)

# Labels and styling
ax.set_xlabel('Fold enrichment')
ax.set_ylabel('')
ax.set_ylim(-0.5, plot_df.shape[0] - 0.5)
sns.despine()
plt.tight_layout()
fig.savefig('GWAS_res_updated_enrichment_fisher_bg_use_nontraitsnp_barplot.pdf')
plt.show()



# ### conservation analysis between mouse loop and human loop
# 

# In[ ]:


## load mouse loop sets
## different mouse loop type
mm_EP_loop = pd.read_csv('../mouse_EP_df.csv', index_col = 0)
mm_PP_loop = pd.read_csv('../mouse_PP_df.csv', index_col = 0)
mm_PO_loop = pd.read_csv('../mouse_PO_df.csv', index_col = 0)
mm_OO_loop = pd.read_csv('../mouse_OO_df.csv', index_col = 0)

## union loops
mm_loop_pk = pk.load(open('../../MicroC/loop_set.pk', 'rb'))
mm_union_loop = mm_loop_pk['union_loops']
mm_union_loop['r1'] = mm_union_loop.iloc[:,0:3].apply(lambda row: '-'.join(row.astype('str').tolist()), axis = 1)
mm_union_loop['r2'] = mm_union_loop.iloc[:,3:6].apply(lambda row: '-'.join(row.astype('str').tolist()), axis = 1)
mm_union_loop.to_csv('union_loop_set_anchors.bed', sep = '\t', index = None, header = None)

## ---- run adaliftover for mouse loop to find human genomic locatons by considering both sequence and epigenomics
#### bash: Rscript adaliftover_run2.R union_loop_set_anchors.bed
### this will generate a file named mm_union_loop_adliftover_to_hg38.txt

## load mouse adaliftover to hg38
mm_to_hg38_set = pd.read_csv('../adaliftover/mm_union_loop_adliftover_to_hg38.txt', sep = '\t', header = None)
mm_to_hg38_set.columns = ['seqnames','start','end','width','strand','epigenome','grammar','score', 'mm_anchor']


# In[ ]:


## loop conservation comparison by comparing liftover human loop location with real human loop sets
lift_term = 'score'
lift_cut = 0.3
gap = 20000
resolution = 5000

mm_loop = mm_union_loop.iloc[:,0:6].copy()

mm_loop['r1'] = mm_loop.iloc[:,0:3].apply(lambda row: '-'.join(row.astype('str').tolist()), axis = 1)
mm_loop['r2'] = mm_loop.iloc[:,3:6].apply(lambda row: '-'.join(row.astype('str').tolist()), axis = 1)

## compare
tmp_set = mm_to_hg38_set[mm_to_hg38_set[lift_term] > lift_cut]
mm_loop['lifted'] = mm_loop['r1'].isin(tmp_set['mm_anchor']) & mm_loop['r2'].isin(tmp_set['mm_anchor'])

## focus on lifted loops
mm_loop_cons = mm_loop.query('lifted == True')
## defined lifted human loop sites
lifted_loops = []
for i, line in mm_loop_cons.iterrows():
    r1, r2 = line['r1'], line['r2']
    r1_set = tmp_set[tmp_set['mm_anchor'] == r1].iloc[0,0:3].tolist()
    r2_set = tmp_set[tmp_set['mm_anchor'] == r2].iloc[0,0:3].tolist()
    lifted_loops.append(r1_set+r2_set+[i])
lifted_loops = pd.DataFrame(lifted_loops, columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'mm_loop'])

## get rid of lifted to different chromosome
lifted_loops = lifted_loops.query('chrom1 == chrom2')
lifted_loops['start1'] = lifted_loops['start1'].astype('int')
lifted_loops['start2'] = lifted_loops['start2'].astype('int')
lifted_loops['end1'] = lifted_loops['end1'].astype('int')
lifted_loops['end2'] = lifted_loops['end2'].astype('int')
lifted_loops = lifted_loops.sort_values(['chrom1', 'start1'])
lifted_loops['r1_center'] = lifted_loops['start1']+int(1/resolution)
lifted_loops['r2_center'] = lifted_loops['start2']+int(1/resolution)

## human loop set
hs_loop = loop_dots_need.iloc[:,0:6].copy()
hs_loop['hs_loop'] = loop_dots_need['label'].tolist()
hs_loop['r1'] = hs_loop.iloc[:,0:3].apply(lambda row: '-'.join(row.astype('str').tolist()), axis = 1)
hs_loop['r2'] = hs_loop.iloc[:,3:6].apply(lambda row: '-'.join(row.astype('str').tolist()), axis = 1)

hs_loop['r1_center'] = hs_loop['start1']+int(1/resolution)
hs_loop['r2_center'] = hs_loop['start2']+int(1/resolution)

## 
signal_conserved = []
for i, line in lifted_loops.iterrows():
    chrom, r1_center, r2_center = line[['chrom1', 'r1_center', 'r2_center']]
    tmp1 = hs_loop[(hs_loop['chrom1'] == chrom) & 
                (hs_loop['start1']-gap+(0.5*resolution) < r1_center) & 
                (hs_loop['end1']+gap-(0.5*resolution) > r1_center) &
                (hs_loop['start2']-gap+(0.5*resolution) < r2_center) & 
                (hs_loop['end2']+gap-(0.5*resolution) > r2_center)]
    tmp2 = hs_loop[(hs_loop['chrom1'] == chrom) & 
                (hs_loop['start1']-gap+(0.5*resolution) < r2_center) & 
                (hs_loop['end1']+gap-(0.5*resolution) > r2_center) &
                (hs_loop['start2']-gap+(0.5*resolution) < r1_center) & 
                (hs_loop['end2']+gap-(0.5*resolution) > r1_center)]
    if tmp1.shape[0] > 0:
        tmp1['mm_loop'] = line['mm_loop']
        signal_conserved.append(tmp1[['hs_loop', 'mm_loop']].values.tolist())
    if tmp2.shape[0] > 0:
        tmp2['mm_loop'] = line['mm_loop']
        signal_conserved.append(tmp2[['hs_loop', 'mm_loop']].values.tolist())

out = open('seq_conservation_data_v2.pk', 'wb')
pk.dump({'hg38_lifted_loops': lifted_loops, 'mouse_loops': mm_loop, 'conserved': signal_conserved,
        'gap': gap, 'lift_term': lift_term, 'lift_cut': lift_cut}, out)
out.close()


## prepare supp table for conserved loop

## AdaLiftOver conservation analysis
lift_term='score'
lift_cut = 0.3

with open('conservation_data_v2.pk', 'rb') as f:
    pk_dat = pk.load(f)

pk_dat['hg38_lifted_loops']['lifted_hs_cord'] = pk_dat['hg38_lifted_loops'].iloc[:,0:6].apply(lambda row: '-'.join(row.astype('str').tolist()), axis = 1).tolist()

conserved_loop_df = []
for x in pk_dat['conserved']:
    conserved_loop_df.extend(x)
conserved_loop_df = pd.DataFrame(conserved_loop_df, columns = ['hs_loop', 'mm_loop'])


tmp = pk_dat['mouse_loops'].sort_values('lifted', ascending = False).drop(['r1', 'r2'], axis = 1)
tmp.columns=['mm_chrom1', 'mm_start1', 'mm_end1', 'mm_chrom2', 'mm_start2', 'mm_end2', 'liftover']
tmp = tmp.merge(pk_dat['hg38_lifted_loops'][['mm_loop', 'lifted_hs_cord']], 
          left_index = True, right_on = 'mm_loop', how = 'left').merge(conserved_loop_df,
                                                                      left_on = 'mm_loop', 
                                                                      right_on = 'mm_loop', how = 'left')

## some liftover to different chromosome need to include in the table for showing
hs_cord_dict = {}
tmp_set = mm_to_hg38_set[mm_to_hg38_set[lift_term] > lift_cut]

for mm_loop in tmp.loc[pd.isna(tmp['lifted_hs_cord'])].query('liftover == True')['mm_loop'].unique().tolist():
    line = pk_dat['mouse_loops'].loc[mm_loop]
    r1, r2 = line['r1'], line['r2']
    r1_set = tmp_set[tmp_set['mm_anchor'] == r1].iloc[0,0:3].tolist()
    r2_set = tmp_set[tmp_set['mm_anchor'] == r2].iloc[0,0:3].tolist()
    hs_cord_dict[mm_loop] = '-'.join([str(int(x)) if 'chr' not in str(x) else x for x in r1_set+r2_set])
tmp['lifted_hs_cord'] = [hs_cord_dict.get(i, j) for i, j in tmp[['mm_loop', 'lifted_hs_cord']].values.tolist()]

## match loop prob
tmp = tmp.merge(loop_dots[['prob_shNT_plusFSK', 'score_shNT_plusFSK', 'prob_shNT_minusFSK',
       'score_shNT_minusFSK', 'prob_KD_plusFSK', 'score_KD_plusFSK', 'label']], 
                left_on = 'hs_loop', right_on = 'label', how = 'left').drop('label', axis = 1)

## find gene and expression change FSK vs Ctnrol
tmp_EP, tmp_PP, tmp_PO, tmp_OO = _define_EP_PP_(loop_dots[loop_dots['label'].isin(tmp['hs_loop'].dropna())], comp = 'NT_plusFSK_vs_minusFSK')
tmp_df = pd.concat([tmp_EP, tmp_PP, tmp_PO, tmp_OO])
tmp_df['hs_Type'] = ['EP']*tmp_EP.shape[0] + ['PP']*tmp_PP.shape[0] + ['PO']*tmp_PO.shape[0] + ['Other']*tmp_OO.shape[0]

## label loop change FC
tmp_df['loop_log2FC_FSK_vs_Cntrl'] = tmp_df.apply(lambda row: np.log2(row['prob_shNT_plusFSK'])-np.log2(np.clip(row['prob_shNT_minusFSK'], a_min=0.0001, a_max = 1)), axis = 1)
tmp_df['loop_log2FC_KD_vs_Cntrl'] = tmp_df.apply(lambda row: np.log2(row['prob_KD_plusFSK'])-np.log2(np.clip(row['prob_shNT_plusFSK'], a_min=0.0001, a_max = 1)), axis = 1)

## rename
tmp_df = tmp_df[['label', 'loop_log2FC_FSK_vs_Cntrl', 'loop_log2FC_KD_vs_Cntrl', 'hs_Type', 'r1_gene', 'r2_gene']]

tmp = tmp.merge(tmp_df, left_on = 'hs_loop', right_on = 'label', how = 'left').drop('label', axis = 1)
tmp = tmp.sort_values('loop_log2FC_FSK_vs_Cntrl').sort_values('hs_loop')

## check percentage of conserved loop
print('% of liftover mouse loop conserved with human loop')
tmp[~pd.isna(tmp['hs_loop'])]['mm_loop'].unique().shape[0] * 100 / tmp[~pd.isna(tmp['lifted_hs_cord'])]['mm_loop'].unique().shape[0]

## prepare table
hs_col = ['prob_shNT_plusFSK', 'score_shNT_plusFSK', 'prob_shNT_minusFSK',
       'score_shNT_minusFSK', 'prob_KD_plusFSK', 'score_KD_plusFSK',
       'loop_log2FC_FSK_vs_Cntrl', 'loop_log2FC_KD_vs_Cntrl', 
       'r1_gene', 'r2_gene', 'r1_gene_log2FC_FSK_vs_Cntrl',
       'r2_gene_log2FC_FSK_vs_Cntrl', 'r1_gene_log2FC_KD_vs_Cntrl',
       'r2_gene_log2FC_KD_vs_Cntrl', 'gene',
       'log2FoldChange_plusFSK_vs_minusFSK', 'padj_plusFSK_vs_minusFSK',
       'log2FoldChange_KD_vs_NT_plusFSK', 'padj_KD_vs_NT_plusFSK']
tmp.columns = ['hs_'+x if x in hs_col else x for x in tmp.columns.tolist()]

## mouse loop
mm_union_loop_comb = pd.concat([mm_EP_loop, mm_PP_loop, mm_PO_loop, mm_OO_loop])
mm_union_loop_comb['loop_type'] = ['EP']*mm_EP_loop.shape[0]+['PP']*mm_PP_loop.shape[0]+['PO']*mm_PO_loop.shape[0]+['OO']*mm_OO_loop.shape[0]

mm_union_loop_comb = mm_union_loop_comb[['prob_shNT_plusCL', 'score_shNT_plusCL', 'prob_shNT_minusCL',
       'score_shNT_minusCL', 'prob_KD_plusCL', 'score_KD_plusCL', 'r1_gene', 'r2_gene', 'Type', 'loop_type']]
mm_union_loop_comb.columns = 'mm_'+mm_union_loop_comb.columns

## combine mouse info to human loop table
mm_union_loop_comb_all = pd.merge(tmp, mm_union_loop_comb, left_on = 'mm_loop', right_index = True, how = 'left')

### to supp table
with pd.ExcelWriter('mouse_loop_human_conservation_supp_table.xlsx') as writer:
    tmp = mm_union_loop_comb_all[['mm_loop', 'mm_r1_gene', 'mm_r2_gene', 'liftover', 'lifted_hs_cord', 'hs_loop', 'hs_r1_gene', 'hs_r2_gene']]
    tmp.to_excel(writer, index = None)
    



