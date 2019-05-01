![Meinter Logo](vignettes/figs/meinter.png)

MeinteR (MEthylation INTERpretation) is an R package that identifies critical differentially methylated sites, based on the following hypothesis: Critical methylation-mediated changes are more likely in genomic regions enriched in cis-acting regulatory elements than in genomic “deserts”. MeinteR calculates the abundance of co-localized elements, such as transcription factor binding sites, tentative splice sites, and other DNA features, such as G-quadruplexes and palindromes that potentially lead to distinct DNA conformational features and rank them with respect to their putative methylation impact.

---

## Package files

[Binary Package](dist/MeinteR_0.99.0.tgz)

[Source Package](dist/MeinteR_0.99.0.tar.gz)

[Package Vignette (PDF)](vignettes/Meinter_vignette.pdf)

[Package Manual (PDF)](vignettes/MeinteR.pdf) 



## Installation

To install MeinteR and all its dependencies install and load devtools. 

**Local installation ** 
First, download the binary distribution of the package, unzip it to a local folder and run the following commands:
```
library(devtools)
package.folder <- "~/MeinteR"
install_local(package.folder)
```
**Install from bitbucket**

MeinteR and its dependencies can be directly installed using the following command:
```
devtools::install_bitbucket("andigoni/meinter", quiet=FALSE)
```

[Use Case 1: Genome-wide association of G-quadruplexes with DNA methylation using TCGA/GEO breast cancer datasets](vignettes/UseCase1.Rmd)

[Use Case 2: Evaluation of genomic signatures on cancer methylation profiles](vignettes/UseCase2.Rmd)

[Use Case 3: Associating genomic DMS signatures with gene expression](vignettes/UseCase3.Rmd)



**Demo run**

```
#Demo Run on sample dataset
library(MeinteR)
rm(list = ls())
#Reorder columns of sample dataset
re.sample <- reorderBed(sample, 1, 2, 3, 5, 4)
#Select sites with delta-beta > 0.50
bed.data <- re.sample[re.sample$score >= 0.50,] 
#Find overlaps with alternative splicing events
altSS <- findAltSplicing(bed.data)
#Find overlaps with splice sites
ss <- findSpliceSites(bed.data, persim = 0.8, offset = 10)
#Find palindromic sequences
pals <- findPals(bed.data)
#Find G-quadruplexes
quads <- findQuads(bed.data, offset = 50)
#Find overlaps with "MA0107.1", "MA0098", "MA115.1", "MA0131.2 transcription factors
tfbs <-
  findTFBS(
    bed.data,
    target = "all",
    tf.ID = c("MA0107.1", "MA0098", "MA115.1", "MA0131.2")
  )
#Find overlaps with conserved transcription factors
ctfbs <-
  findConservedTFBS(bed.data)
#Find overlaps with DNA shape features
shapes <- findShapes(bed.data)
#Weighting scheme
weights = list()
weights[["spls"]] = 1
weights[["ctfbs"]] = 1
weights[["tfbs"]] = 1
weights[["pals"]] = 1
weights[["quads"]] = 1
weights[["shapes"]] = 1
#Mapping MeinteR arguments with function outputs
funList = list()
funList[["spls"]] = ss
funList[["altss"]] = altSS
funList[["tfbs"]] = tfbs
funList[["ctfbs"]] = ctfbs
funList[["pals"]] = pals
funList[["quads"]] = quads
funList[["shapes"]] = shapes
#Calculate genomic index
index <- meinter(re.sample, funList, weights)


```






---





