# Interval basics

## Task 1: Implement functions to count overlaps and calculate Jaccard index

Implement the 3 functions in [intervals.R](src/intervals.R). Descriptions of the functions can be found there. Additional instructions for each function follow:

- `countOverlapsSimple`: your function may use a sequential search to identify overlaps (no need to implement binary search or lower complexity algorithm). You should write this function using base R only, without relying on GenomicRanges or other Bioconductor packages. The inputs should be base R objects (a list for the query, and a data.frame for the database, see examples in the function documentation). (20 pts)

- `calculateJaccardScore`: for this function, you may make use of the Bioconductor GenomicRanges package, which includes fast implementations of interval overlap and union. In particular, you may be interested in the `intersect` and `union` methods that operate on `GenomicRanges::GRanges` objects. The inputs to this function should both be `GenomicRanges::GRanges` objects. (20 pts)

- `pairwiseJaccard`: Input should be a base R list of GRanges objects, and you should return a square matrix with dimension equal to the number of elements in the input list. (10 pts)

## Task 2: Calculate and interpret similarity results

There are 4 narrowPeak files stored under [/data](/data). These are ChIP-seq datasets, for two cell types (HepG2 and HeLaS3) and two transcription factors (Jun and CTCF). I filtered the data to keep things smaller and faster for this assignment.

Use your functions to calculate the pairwise Jaccard similarity for this set of 4 experiments, then answer the following questions:

1. Which two interval sets are the most similar? (10 pts)
2. Which two interval sets are the most different? (10 pts)
3. Based on these results, which factor, CTCF or Jun, would you predict varies more across cell types? (10 pts)
4. Based on these results, do the genomic locations found by ChIP-seq experiments depend more on the cell-type, or on the transcription factor being assayed? (20 pts)

## Testing:

You can test your work by running within the interval_basics directory:

```
Rscript tests/testDriver.R
```
