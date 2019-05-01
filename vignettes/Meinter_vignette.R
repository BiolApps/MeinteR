## ---- message=FALSE, eval=FALSE------------------------------------------
#  library(MeinteR)

## ---- eval=FALSE---------------------------------------------------------
#  nameStudy(study.name="MyProj")

## ---- message=FALSE, eval=FALSE------------------------------------------
#  project.dir <- "~/meinter_dir"

## ---- message=FALSE, eval=FALSE------------------------------------------
#  head(sample)

## ---- eval=FALSE, message=FALSE------------------------------------------
#  re.sample = reorderBed(sample,1,2,3,5,4)
#  head(re.sample)

## ---- eval=FALSE---------------------------------------------------------
#  input.data <- read.csv(file.path(project.dir, "my_data.csv"), sep=",", header = T)

## ---- eval=FALSE---------------------------------------------------------
#  fp <- file.path(project.dir, "GSE37362_annotation.csv")
#  geo.data <- importGEO(gse.acc="GSE37362", annotation.file= fp)

## ---- eval=FALSE, message=FALSE------------------------------------------
#  #Select DMS with delta-beta values equal to 0
#  subsample.1 <- re.sample[re.sample$score == 0,]
#  #Select DMS with delta-beta values greater than or equal to 0.30
#  subsample.2 <- re.sample[re.sample$score >= 0.30,]
#  #Select DMS with absolute delta-beta values greater than or equal to 0.60
#  subsample.3 <- re.sample[abs(re.sample$score) >= 0.60,]

## ---- message=FALSE, eval=FALSE------------------------------------------
#  #Transcription factors of interest
#  tf.ID = c("MA0003.1", "MA0019.1", "MA0004.1", "MA0036.3", "MA0037.3")
#  tfbs <-findTFBS(bed.data=subsample.3,persim=0.8, offset=10, target="PROMOTER",
#                  up.tss=2000, down.tss=100, mcores = 2, tf.ID=tf.ID)

## ---- eval=FALSE---------------------------------------------------------
#  plotTF(tfbs[[1]], topTF=10) #topTF:Number of most frequent transcription factors

## ---- message=FALSE, eval=FALSE------------------------------------------
#  ctfbs <-findConservedTFBS(subsample.3, known.conserved.tfbs.file="~/tfbsConsSites.gz")

## ---- eval=FALSE---------------------------------------------------------
#  scatterConsTF(ctfbs[[2]])

## ---- message=FALSE, eval=FALSE------------------------------------------
#  ss <- findSpliceSites(bed.data=subsample.3, persim=0.8, offset= 10)

## ---- eval=FALSE, message=FALSE------------------------------------------
#  head(ss[[1]])
#  head(ss[[2]])

## ---- eval=FALSE, message=FALSE------------------------------------------
#  altss <- findAltSplicing(subsample.3)

## ---- message=FALSE, eval=FALSE------------------------------------------
#  head(altss[[1]])
#  head(altss[[2]])
#  head(altss[[3]])
#  altss[[4]] #Plot alternative splicing events

## ---- eval=FALSE---------------------------------------------------------
#  pals <- findPals(bed.data=subsample.3, offset=10, min.arm=5, max.loop=5, max.mismatch=1)

## ---- eval=FALSE, message=FALSE------------------------------------------
#  head(pals[[1]],n=1)
#  pals[[2]] # On/Off DMS palindromes
#  head(pals[[3]]) #Number of palindromes per DMS

## ---- eval=FALSE, message=FALSE------------------------------------------
#  #Detect G4s in the 100nt region flanking DMS
#  quads <- findQuads(bed.data=subsample.3, offset=100)
#  quads[[1]] # G4 locus information for each sequence
#  quads[[2]] # G4 on/neighboring input data
#  quads[[3]] # Number of G4 per sequence

## ---- eval=FALSE, message=FALSE------------------------------------------
#  #Detect DNA shapes in the 100nt region flanking DMS in the re.sample data
#  shapes <- findShapes(bed.data=subsample.3, offset=100)

## ----call_genomic_index, eval=FALSE, message=FALSE-----------------------
#  #Calculate genomic index
#  weights = list()
#  weights[["spls"]] = 1
#  weights[["ctfbs"]] = 2
#  weights[["tfbs"]] = 2
#  weights[["pals"]] = 1
#  weights[["quads"]] = 2
#  weights[["shapes"]] = 1
#  funList = list()
#  funList[["spls"]] = ss
#  funList[["altss"]] = altss
#  funList[["tfbs"]] = tfbs
#  funList[["ctfbs"]] = ctfbs
#  funList[["pals"]] = pals
#  funList[["quads"]] = quads
#  funList[["shapes"]] = shapes
#  index <- meinter(subsample.3, funList, weights)
#  head(index) # Highly ranked DMS

## ---- message=FALSE, eval=FALSE------------------------------------------
#  dev.off() #close any open plot device
#  res <- plotCpG(bed.data=re.sample, offset=100)

## ---- message=FALSE, eval=FALSE------------------------------------------
#  plotBeta(bed.data=re.sample)

## ---- message=FALSE, eval=FALSE------------------------------------------
#  project.dir <- "~/GSE37362"
#  setwd(project.dir)
#  library(MeinteR)
#  gse.accession <- "GSE37362"
#  annotation.file <- file.path(project.dir,"GSE37362_annotation.csv")
#  geo.data <- importGEO(gse.acc=gse.accession, annotation.file=annotation.file)

## ---- message=FALSE, eval=FALSE------------------------------------------
#  # Re-order columns and omit rows with empty cells
#  bed.data<-na.omit(reorderBed(geo.data[[1]],3,4,5,2))

## ---- message=FALSE, eval=FALSE------------------------------------------
#  #Select probes with delta-beta < -0.2
#  sub.data <- subset(bed.data, bed.data$score < -0.2) # 713 probes
#  quads <- findQuads(sub.data, offset=100)
#  pals <- findPals(sub.data, offset=100)
#  altspl <- findAltSplicing(sub.data)
#  tfbs <- findTFBS(sub.data, target="PROMOTER", up.tss=5000, down.tss=100)
#  ss <- findSpliceSites(sub.data)
#  shapes <- findShapes(sub.data, offset=100)

## ---- message=FALSE, eval=FALSE------------------------------------------
#  #Transform all the bed.gz files in the ~/GSE69272_RAW folder
#  files <- list.files(path="~/GSE69272_RAW", pattern="*.bed.gz",
#                      full.names=T, recursive=FALSE)
#  all.samples <- lapply(files, function(samples) {
#      loadSeqGEO(file.path=samples, cov=30, chroms="chr19")
#  })
#  

## ---- message=FALSE, eval=FALSE------------------------------------------
#  library(MeinteR)
#  rm(list = ls())
#  re.sample <- reorderBed(sample, 1, 2, 3, 5, 4)
#  bed.data <- re.sample[re.sample$score >= 0.50,]
#  altSS <- findAltSplicing(bed.data)
#  ss <- findSpliceSites(bed.data, persim = 0.8, offset = 10)
#  pals <- findPals(bed.data)
#  quads <- findQuads(bed.data, offset = 50)
#  tfbs <-
#    findTFBS(
#      bed.data,
#      target = "all",
#      tf.ID = c("MA0107.1", "MA0098", "MA115.1", "MA0131.2")
#    )
#  ctfbs <-
#    findConservedTFBS(bed.data, known.conserved.tfbs.file = "~/Downloads/tfbsConsSites.txt.gz")
#  shapes <- findShapes(bed.data)
#  weights = list()
#  weights[["spls"]] = 1
#  weights[["ctfbs"]] = 1
#  weights[["tfbs"]] = 1
#  weights[["pals"]] = 1
#  weights[["quads"]] = 1
#  weights[["shapes"]] = 1
#  funList = list()
#  funList[["spls"]] = ss
#  funList[["altss"]] = altSS
#  funList[["tfbs"]] = tfbs
#  funList[["ctfbs"]] = ctfbs
#  funList[["pals"]] = pals
#  funList[["quads"]] = quads
#  funList[["shapes"]] = shapes
#  index <- meinter(re.sample, funList, weights)

