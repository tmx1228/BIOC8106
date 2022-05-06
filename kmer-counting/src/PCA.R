## load function
source("~/BIOC8061/kmer-counting/src/countKmers.R")

## load data
setwd("~/BIOC8061/kmer-counting/fasta")
filesList = c(
  #"all.fa",
  "fos.fa",
  "fosb.fa",
  "fosl1.fa",
  "fosl2.fa",
  "jun.fa",
  "junb.fa",
  "jund.fa",
  "s288c.fa"
  )

## set count matrix
allpos4mer <- c()
for(i1 in c("A","C","G","T")){
  for(i2 in c("A","C","G","T")){
    for(i3 in c("A","C","G","T")){
      for(i4 in c("A","C","G","T")){
        allpos4mer <- c(allpos4mer, paste0(i1,i2,i3,i4))
}}}}
a <- matrix(rep(0, length(filesList)*4**4),nrow=4**4)
colnames(a) <- filesList
rownames(a) <- allpos4mer

## count kmer
for(file in filesList){
  data <- read.delim(file, comment.char = ">") %>% 
    gsub("[[:space:]]", "", .) %>% 
    gsub("[\",c()]","", .)
  kmers.count <- countKmers(data[1],4)
  for(kmer in names(kmers.count)){
    a[kmer,file] = kmers.count[kmer]
  }
}

head(a)

## log2 scale to fit the linear model of PCA
normalized.peaks <- log2(a+1)

## prcomp requires rows to be samples and columns to be features
pca.input <- t(normalized.peaks) %>% .[ , which(apply(., 2, var) != 0)]
pca.input[,1:4]
dim(pca.input)

## use build in normalization method of prcomp: zscore
pca <- prcomp(pca.input, scale. = T)

summary(pca) 
#variance 0.9371 0.02994

## generate data frame for plot
df <- data.frame(pca$x[,1:2])

## plot PCA
pdf("Rplot_PCA.pdf", width = 5, height = 3.5)
ggplot(df, aes(x=PC1, y=PC2, label=rownames(df))) +
  geom_point() +
  theme +
  geom_text_repel(max.overlaps = 15, size = 2) +
  xlab("PC1(93.71% var)")+ylab("PC2(2.99% var)") +
  ggtitle("log2(count + 1)")


## no normalization
pca.input <- t(a)
## use build in normalization method of prcomp: zscore
pca <- prcomp(pca.input, scale. = T)

summary(pca) 
#variance 0.9894 0.00939

## generate data frame for plot
df <- data.frame(pca$x[,1:2])

## plot PCA
ggplot(df, aes(x=PC1, y=PC2, label=rownames(df))) +
  geom_point() +
  theme +
  geom_text_repel(max.overlaps = 15, size = 2) +
  xlab("PC1(98.94% var)")+ylab("PC2(0.94% var)") + 
  ggtitle("raw count")
dev.off()



