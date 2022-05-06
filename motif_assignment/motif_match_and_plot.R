setwd("~/BIOC8061/motif_assignment/")
library(ggplot2)
library(reshape2)

mat <- read.table("rbpj_pwm.txt", header = F)
head(mat)
colnames(mat) <- c("A","C","G","T")
mat
mat.reverse <- mat[nrow(mat):1,]
colnames(mat.reverse) <- c("T","G","C","A")
rownames(mat.reverse) <- 1:8

mat.bg <- data.frame(0.28,0.22,0.22,0.28)
colnames(mat.bg) <- c("A","C","G","T")

fa <- readLines("notch1_binding_sequences.fa")

pdf("fa_scores.pdf")
for (line in fa){
  #line = "TCCTGCACCACCGTGAATTTAATAGCCTGGAGATCCTGTTTCACTTTCTGAGAGAACAGATCCTATAAATTGCATTCTATGTGCCATTCCACTTCCTGTCTGATGGTGATGCGGAAAAACAGCTGTTGTGTGGGAACTGTTCAGAGCTGTGTTTCTCAAACCACAGGTCATGGAAGGGTGGTTAAGTAGGCCCTGTTTCCCTTTCTTCCAAAATAGCACCCAAGGCTCTTCTGACCTTGCCTCCTTCCACCTCTTAGTTGTAAGCCTACTAAGGAGTAAGTTTATCCTACTTCCTTTTGT"
  if (substring(line,1,1) == ">"){
    mytitle <- substring(line,2,nchar(line))
  } else if (substring(line,1,1) != ">"){
    #print(line)
    scores.fa <- data.frame(0,0,0)
    colnames(scores.fa) <- c("pos","score","score.reverse")
    
    for (j in 1:(nchar(line)-7)){
      mer <- substring(line,j,j+7)
      #print(c(mer,score.mer,log2(score.mer),max(0, log2(score.mer))))  
      score.mer = 1
      score.mer.reverse = 1
      #mer="ACGTTCGA"
      
      for (i in 1:8){
        #i=2
        bp <- substring(mer,i,i)
        
        if (bp == "A"){
          score.bp = mat[i,"A"]/mat.bg[1,"A"]
          score.mer = score.mer*score.bp
          score.bp.reverse = mat.reverse[i,"A"]/mat.bg[1,"A"]
          score.mer.reverse = score.mer.reverse*score.bp.reverse
        } else if (bp == "T"){
          score.bp = mat[i,"T"]/mat.bg[1,"T"]
          score.mer = score.mer*score.bp
          score.bp.reverse = mat.reverse[i,"T"]/mat.bg[1,"T"]
          score.mer.reverse = score.mer.reverse*score.bp.reverse
        } else if (bp == "C"){
          score.bp = mat[i,"C"]/mat.bg[1,"C"]
          score.mer = score.mer*score.bp
          score.bp.reverse = mat.reverse[i,"C"]/mat.bg[1,"C"]
          score.mer.reverse = score.mer.reverse*score.bp.reverse
        } else if (bp == "G") {
          score.bp = mat[i,"G"]/mat.bg[1,"G"]
          score.mer = score.mer*score.bp
          score.bp.reverse = mat.reverse[i,"G"]/mat.bg[1,"G"]
          score.mer.reverse = score.mer.reverse*score.bp.reverse
        
        }
        #print(c(mer,score.bp.reverse, score.mer.reverse))
      }
      
      score.mer <- max(0, log2(score.mer))
      score.mer.reverse <- max(0, log2(score.mer.reverse))
      scores.fa[j,] <-  c(j,score.mer,score.mer.reverse)
      
    }
    
    df <- data.frame(scores.fa) %>% melt(., id.var = "pos",variable.name = "reads", value.name = "motif matching score")
    p <- ggplot(df, aes(x = pos, y= `motif matching score`, color=reads)) +
      geom_line() +
      theme +
      ggtitle(mytitle) +
      theme(axis.title.x = element_blank())
    print(p)

  }
    
}

dev.off()
