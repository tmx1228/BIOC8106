## import the function
source("~/BIOC8061/interval_basics/src/intervals.R")

## set working directory
setwd("~/BIOC8061/interval_basics")

## read the files
filesList = c(
  "data/helas3_ctcf.narrowPeak.gz",
  "data/helas3_jun.narrowPeak.gz",
  "data/hepg2_ctcf.narrowPeak.gz",
  "data/hepg2_jun.narrowPeak.gz"
)

files = lapply(filesList, read.table)

## load files in GRanges format
granges = lapply(files, function(x) {
  GRanges(seqnames=x$V1, ranges=IRanges(x$V2, x$V3))
})

## calculate pairwise Jaccard scores
pairwiseJaccard(granges)
