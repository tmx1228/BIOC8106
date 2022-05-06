#### Outline
- ChIP-seq technology and development
- ChIP-seq data analysis 
    - Strategy
    - Peak calling (MACS) 
    - Peak calling (SICER)

-----------------------------------------------------
#### The goal of ChIP-seq: To determine the locations in the genome associating with a protein factor
![f38a69f1f487fde2c0df94a38092c5e2.png](evernotecid://8180CF8E-FC94-42BE-985A-2C0C8231A501/wwwevernotecom/226309320/ENResource/p7)

Where in the genome does your factor of interest locate?

-----------------------------------------------------

#### ChIP-seq experiment steps:

1. Chromatin ImmunoPrecipitation (ChIP)
2. Protein-DNA crosslinking with formaldehyde (for TF)
3. Chop the chromatin using sonication (TF) or micrococal nuclease (MNase) digestion (histone)
4. Specific factor-targeting antibody
5. Immunoprecipitation
6. DNA purification
7. PCR amplification (~150bp)
8. High-throughput sequencing (Illumina: can only sequencing the end of the DNA fragments)
![d95ffeaa9f23ffe210983d452c89c9a3.png](evernotecid://8180CF8E-FC94-42BE-985A-2C0C8231A501/wwwevernotecom/226309320/ENResource/p8)


-----------------------------------------------------

#### History of the development of ChIP-seq technology
- UV crosslinking (1984) : the protein-DNA interaction can be captured
- Crosslinking + immunoprecipitation (1993) : use antibody to grab the DNA-protein complex
- ChIP-chip (2000) : **genomewide** microarray method was developed, using a **pre-designed** way
- Unbiased chromosomal coverage by tiling array (2004)
- ChIP-seq (2007)

![2c5a5b650b7643484ec0436e4803d405.png](evernotecid://8180CF8E-FC94-42BE-985A-2C0C8231A501/wwwevernotecom/226309320/ENResource/p9)

- ChIP-seq has become a predominant method for profiling chromatin epigenomes


-----------------------------------------------------

#### ChIP-seq data analysis
- Goals: 
• Where in the genome do these sequence reads come from? - Sequence alignment and quality control
• What does the enrichment of sequences mean? - Peak calling
• What can we learn from these data? – Downstream analysis and integration

- Steps：
![2c051698faaec89263d5dcfa9c2039be.png](evernotecid://8180CF8E-FC94-42BE-985A-2C0C8231A501/wwwevernotecom/226309320/ENResource/p10)

    1. Sequencing quality assessment: fastqc. If the quality scores across bases fail, either re-do the experiment or trim the data.
    2. ChIP-seq read mapping: map the fastq file containing the sequence information to the genome; alignment of each sequence read: bowtie, BWA (Burrows–Wheeler Algorithm); usually use the reads can map to a unique/best location in the genome.
    3. Redundancy control: complete identical reads are considered error (for example, induced by PCR)
    Non-redundant rate:# non-redundant reads / # mapped reads
    PBC (PCR Bottleneck Coefficient): # locations w/ 1 read # locations w/ reads
    4. DNA fragment size estimation: 
    peak model (MACS) for TF
    cross-correlation (SICER) for any ChIP-seq (input): calculate the Correlation between two strings with a displacement; Auto-correlation: Cross-correlation with itself
    5. Retrieve DNA fragments
    Full length retrieval (MACS)
    Partial retrieval (sharpen the signal)
    Point retrieval (SICER)
    6. Pile up: Signal map generation

#### ChIP-seq: Study design
Background Control: Input or IgG
- Input chromatin: sonicated/digested chromatin without immunoprecipitation
- IgG: “unspecific” immunoprecipitation

Study Control:
- Control exp sample: ChIP + input 
- Treated exp sample: ChIP + input


#### ChIP-seq: Peak calling
Goal: Identify regions in the genome enriched for sequence reads:
– Compared to genomic background
– Compared to input control

MACS: Model-based Analysis for ChIP-Seq
Read distribution along the genome 
    ~ Poisson distribution (λBG = total tag / genome size)
    ~ Negative binomial distribution (MACS2)

ChIP-seq show local biases in the genome
– Chromatin and sequencing bias
– 200-300bp control windows have too few tags – But can look further

B-H adjustment to correct for FDR : p-value → q-value

Data Visualization
• bedGraph to bigWig 
• macs2 output data
• IGV

Quality Control
• FRiP(FractionofReadsinPeaks)score – 1-10% for TF is normal
• Numberofpeaks
    - Number of peaks with high fold-enrichment, e.g, 5, 10, ... 
    - 2000
• Sequenceconservation
• Fractionofpeakswithinregulatoryregions – 80%

Biological interpretation: ChIP-seq captures a snapshot of binding patterns from a cell population
• TF intrinsic property
• Binding activity
• Cellular heterogeneity

![fec44729da02eaf31c2d6e8a9d9a0d41.png](evernotecid://8180CF8E-FC94-42BE-985A-2C0C8231A501/wwwevernotecom/226309320/ENResource/p11)
