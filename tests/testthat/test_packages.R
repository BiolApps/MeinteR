rm(list=ls())
packages <- c("ggplot2", "parallel", "plyr","stats4","reshape2", "pqsfinder","BSgenome.Hsapiens.UCSC.hg19", "GenomicRanges", 
              "BSgenome","tidyverse","Biostrings","XVector","GenomeInfoDb","S4Vectors","stats4", "stats","IRanges","BiocGenerics",
              "parallel","base","rtracklayer","reshape2","TFBSTools","DNAshapeR","JASPAR2018" )
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

bioconductor.packages <- c("BSgenome.Hsapiens.UCSC.hg19","GenomicRanges","rtracklayer", "BiocGenerics","S4Vectors", "IRanges","GenomeInfoDb","BSgenome","Biostrings","XVector")
bioc.packages <- bioconductor.packages[!(bioconductor.packages  %in% installed.packages()[,"Package"])]
if (length(bioc.packages)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(bioc.packages)
}