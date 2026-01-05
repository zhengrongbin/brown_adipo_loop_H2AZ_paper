args = commandArgs(trailingOnly=TRUE)

library(AdaLiftOver)
library(GenomicRanges)
library(BiocParallel)      # ← for parallel map
library(rtracklayer)       # ← import.chain
library(Biostrings)

#loops = read.csv('../union_loop_set_anchors.bed', sep = '\t', header = F)
## read loop set in mouse for liftover
loops = read.csv(args[1], sep = '\t', header = F)
loops$names = apply(loops, 1, function(x) gsub(' ', '', paste0(as.vector(x), collapse = '-')))

## file to save lift over regions
done = read.csv('mm_union_loop_adliftover_to_hg38.txt', sep = '\t', header = F)
loops = subset(loops, !names %in% as.vector(done$V9))
	
data("epigenome_mm10")
data("epigenome_hg38")
data("jaspar_pfm_list")

chain <- rtracklayer::import.chain("mm10.hg38.rbest.chain")

loops$V4 = '.'
gr <- GRanges(seqnames = loops$V1,
           ranges = IRanges(start = loops$V2,
                            end = loops$V3,
                            names = loops$V4))

run_one <- function(gr_i, idx) {
  ## 1) Liftover
  lifted <- adaptive_liftover(gr_i, chain)

  ## 2) Epigenome similarity
  lifted <- compute_similarity_epigenome(gr_i,
                                         lifted,
                                         epigenome_mm10,
                                         epigenome_hg38)

  ## 3) Grammar similarity
  lifted <- compute_similarity_grammar(gr_i,
                                       lifted,
                                       "mm10","hg38",
                                       jaspar_pfm_list)

  ## 4) Filter candidates
  filtered <- gr_candidate_filter(lifted,
                                  best_k = 1L,
                                  threshold = 0)

  ## 5) If any hits, write them
  df <- data.frame(filtered@unlistData)
  if (nrow(df) > 0) {
    df$mm_seq <- paste(loops[idx,1:3], collapse = "-")
    write.table(df,
                file      = 'mm_union_loop_adliftover_to_hg38.txt',
                append    = TRUE,
                row.names = FALSE,
                col.names = FALSE,
                sep       = '\t',
                quote     = FALSE)
    return('1')
  } else {
    return('0')
  }
}

## A function to process one GRanges at index i
process_one <- function(i) {
  gr_i <- gr[i]
  message("Working on index ", i)

  out <- tryCatch({
    run_one(gr[i], i)
  }, error = function(e) {
    message("Error at index ", i, ": ", e$message)
    NULL
  })
  return(out)
}

## Set up a MulticoreParam with, say, 4 workers
param <- MulticoreParam(workers = 16,
                        progressbar = TRUE)

## Run in parallel over all indices
res_list <- bplapply(seq_along(gr),
                     FUN     = process_one,
                     BPPARAM = param)
