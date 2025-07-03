import os
import pandas as pd
import numpy as np
import pybedtools
from pathlib import Path

# 1) Load your lead SNP list
with open('hg38_union_loop_leadSNP_100kwindow.txt') as f:
    lead_snps = [l.strip() for l in f]

# 2) Gather all valid LD files and read into one DF
ld_records = []
for snp in lead_snps:
    fpath = Path(f'./LD_snp/{snp}.txt')
    if not fpath.exists():
        continue
    # quick error check: only read first KB
    head = fpath.open().read(1024)
    if 'error' in head:
        continue

    tmp = pd.read_csv(fpath, sep='\t')
    tmp = tmp[tmp.RS_Number != snp]
    tmp = tmp[tmp.R2 > 0.8]
    if tmp.empty:
        continue

    # annotate origin and split coords
    tmp = tmp.assign(
        lead_snp = snp,
        chrom    = tmp['Coord'].str.split(':').str[0],
        start    = tmp['Coord'].str.split(':').str[1].astype(int),
    )
    tmp['end'] = tmp['start'] + 1

    # keep only what we need
    ld_records.append(tmp[['chrom','start','end', 'lead_snp','RS_Number']])

if ld_records:
    ld_df = pd.concat(ld_records, ignore_index=True)
else:
    ld_df = pd.DataFrame(columns=['chrom','start','end', 'lead_snp','RS_Number'])

# 3) One-shot intersection vs. your loop anchors
anchors = pd.read_csv('human_union_loop_anchors.txt', sep='\t')
# ensure anchors has columns chrom, start, end, name
anchors_bed = pybedtools.BedTool.from_dataframe(anchors)

# prepare ld bed with two extra fields:
ld_bed = pybedtools.BedTool.from_dataframe(ld_df)

intersected = (
    ld_bed
    .intersect(anchors_bed, wa=True)
    .to_dataframe()
    .drop_duplicates()
)

# pybedtools will call the 4th column 'name' and 5th 'score', so rename:
intersected = intersected.rename(columns={'name':'lead_snp', 'score':'RS_Number'})

# get unique lead→LD pairs
lead_ld_pairs = intersected[['lead_snp','RS_Number']].drop_duplicates()

# 4) Read and filter GWAS once
gwas = pd.read_csv('../../MicroC/hg38_GWAScatelog_UCSC.tsv', sep='\t')
gwas = (
    gwas.query('pValue < 1e-8')
        .drop_duplicates(['chrom','chromStart','chromEnd','name','trait'])
        .assign(label = lambda df: (
            df.chrom.astype(str) + ':' +
            df.chromStart.astype(str) + '-' +
            df.chromEnd.astype(str)   + '_' +
            df.name + '_' + df.trait
        ))
)

# core GWAS output
core = gwas[['chrom','chromStart','chromEnd','label','name','trait']]

# 5) Expand: join lead SNPs to their LD proxies, then pull coords & build labels
expanded = (
    core[['trait','name']]            # trait & lead SNP
        .merge(lead_ld_pairs, left_on='name', right_on='lead_snp')
        .merge(ld_df, on=['lead_snp','RS_Number'])
        .assign(
            label = lambda df: (
                df.chrom.astype(str) + ':' +
                df.start.astype(str) + '-' +
                df.end.astype(str)   + '_' +
                df.RS_Number + '_' + df.trait
            )
        )
        .rename(columns={'start':'chromStart','end':'chromEnd'})
        [['chrom','chromStart','chromEnd','label', 'lead_snp','RS_Number']]
)

expanded.drop_duplicates().to_csv('hg38_GWAS_LD_expanded.txt', sep='\t', index=False)

# 6) final union and write
out = pd.concat(
    [ core[['chrom','chromStart','chromEnd','label']],
      expanded[['chrom','chromStart','chromEnd','label']] ],
    ignore_index=True
).drop_duplicates()

out.to_csv('hg38_GWAS_LD_expanded_and_core_gwas.txt', sep='\t', index=False)
