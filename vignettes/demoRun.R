


#Demo Run on sample dataset
library(MeinteR)
rm(list = ls())
re.sample <- reorderBed(sample, 1, 2, 3, 5, 4)
bed.data <- re.sample[re.sample$score >= 0.50,] 
altSS <- findAltSplicing(bed.data)
ss <- findSpliceSites(bed.data, persim = 0.8, offset = 10)
pals <- findPals(bed.data)
quads <- findQuads(bed.data, offset = 50)
tfbs <-
  findTFBS(
    bed.data,
    target = "all",
    tf.ID = c("MA0107.1", "MA0098", "MA115.1", "MA0131.2")
  )
ctfbs <-
  findConservedTFBS(bed.data, known.conserved.tfbs.file = "~/Downloads/tfbsConsSites.txt.gz")
shapes <- findShapes(bed.data)
weights = list()
weights[["spls"]] = 1
weights[["ctfbs"]] = 1
weights[["tfbs"]] = 1
weights[["pals"]] = 1
weights[["quads"]] = 1
weights[["shapes"]] = 1
funList = list()
funList[["spls"]] = ss
funList[["altss"]] = altSS
funList[["tfbs"]] = tfbs
funList[["ctfbs"]] = ctfbs
funList[["pals"]] = pals
funList[["quads"]] = quads
funList[["shapes"]] = shapes
index <- meinter(re.sample, funList, weights)
