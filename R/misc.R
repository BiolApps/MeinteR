#############################################################################
## Package    : MeinteR
## File       : misc.R
## Functions  : loadFile (exported)
##            : importLimma (exported)
##            : bed2Seq (exported)
##            : nameStudy (exported)
##            : isEmptyDF
##            : reorderBed (exported)
##            : filterByProm
##            : filterByCGI
##            : is.defined
##            : validateBed
##            : plotCpG (exported)
##            : plotBeta (exported)
##            : importGEO  (exported)
##            : loadAnnotationGEO
##            : validateGEO
##            : groupGEO
##            : probesChr
##            : loadSeqGEO (exported)
##            : mval2beta
## Updated    : 20-08-2019
##
## Title      : Miscellaneous functions
#############################################################################


rm(list = ls())
#Package local environmental variables
pkg.env <- new.env()
local(NAME <- "study", envir = pkg.env)

#Install package dependencies
packages <-
  c(
    "ggplot2",
    "plyr",
    "stats4",
    "reshape2",
    "tidyverse",
    "XVector",
    "GenomeInfoDb",
    "S4Vectors",
    "stats",
    "IRanges",
    "BiocGenerics",
    "parallel",
    "base",
    "reshape2",
    "grDevices",
    "graphics",
    "utils"
  )
new.packages <-
  packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new.packages))
  install.packages(new.packages)

bioconductor.packages <-
  c(
    "TFBSTools",
    "TxDb.Hsapiens.UCSC.hg19.knownGene",
    "Biostrings",
    "BSgenome",
    "pqsfinder",
    "Biobase",
    "DNAshapeR",
    "JASPAR2018",
    "BSgenome.Hsapiens.UCSC.hg19",
    "GenomicRanges",
    "rtracklayer",
    "BiocGenerics",
    "S4Vectors",
    "JASPAR2018",
    "IRanges",
    "pqsfinder",
    "GenomeInfoDb",
    "BSgenome",
    "GenomicFeatures",
    "XVector",
    "GEOquery",
    "FDb.InfiniumMethylation.hg19",
    "IlluminaHumanMethylationEPICanno.ilm10b2.hg19"
  )
bioc.packages <-
  bioconductor.packages[!(bioconductor.packages  %in% installed.packages()[, "Package"])]
if (length(bioc.packages)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(bioc.packages)
}



#'Load input data
#'
#' Loads tabular files containing methylation data. The function checks the delimiter and validates the order
#' of the columns (chr, start,end,score,strand).
#'
#' @param   FH   Full path of the tabular methylation data
#' @return  df   A data frame with the tabular methylation data
#' @export

loadFile <- function(FH) {
  title.line <- readLines(FH, n = 1)
  if (grepl(";", title.line)) {
    bed.data <- read.csv(file = FH,
                         header = TRUE,
                         sep = ";")
  }
  if (grepl(",", title.line)) {
    bed.data <- read.csv(file = FH,
                         header = TRUE,
                         sep = ",")
  }
  if (grepl(",\t", title.line)) {
    bed.data <- read.csv(file = FH,
                         header = TRUE,
                         sep = "\t")
  }
  message("Valid column names and ordering: chr,start,end,score,strand")
  message("Use reorderbed function to change the column names and order.")
  return(bed.data)
}

#'Imports the results of differential methylation analyses exported by \code{limma}
#'
#'
#' @param ltop.obj Output of \code{limma} differential analysis. The output data frame produced by the \code{limma::topTable} function.
#' @param platform A string corresponding to the human methylation array. Possible values are "hm27k", "hm450k" and "EPIC" (default:hm450k)
#' @param sortBy The criterion for selecting probes of the MArrayLM object. Possible values are: "logFC", "P.Value", "adj.P.Val"
#' @return  df A well-formatted data frame to be used as input to the MeinteR workflow.
#' @export
 
importLimma <- function(ltop.obj, platform = "hm450k", sortBy="adj.P.Val") {
  platform.list <- c("hm450k", "hm27k", "epic")
  ltop.obj.cols <- c("logFC","AveExpr","t","P.Value","adj.P.Val","B")
  if (!exists("platform")){platform = "hm450k"}
  if (!exists("sortBy")){sortBy="adj.P.Val"}
  `%notin%` <- Negate(`%in%`)
  if (tolower(platform) %notin% (platform.list)) {stop("Invalid platform. Accepted platforms: hm27k, hm450k and epic.")}
  if (sum(colnames(ltop.obj) %notin% (ltop.obj.cols))>0) {stop("Invalid input format. Accepted input: limma::topTable data frame")}
  if (sortBy %notin% (ltop.obj.cols)){stop("Invalid value. Accepted inputs: limma::topTable data frame columns")}
    if (tolower(platform)==platform.list[1]){ 
    hm450k <- as.data.frame(get450k())
    platform.probes <- data.frame(hm450k[, 1:3], strand="+")
    colnames(platform.probes) = c("chr", "start","end","strand")
  }
  if (tolower(platform)==platform.list[2]){ 
    hm27k <- as.data.frame(get27k())
    platform.probes <- data.frame(hm27k[, 1:3], strand="+")
    colnames(platform.probes) = c("chr", "start","end","strand")
    }
  if (tolower(platform)==platform.list[3]){      
    epic = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    platform.probes = data.frame(epic[,1:2], epic[,2]+1, epic[,3])
    colnames(platform.probes) = c("chr", "start","end","strand")
  }
    if (sum(rownames(ltop.obj) %notin% rownames(platform.probes))>0) {
      inv.probes <- which(rownames(ltop.obj) %notin% rownames(platform.probes)==TRUE)
      if (length(inv.probes)>100) {stop("More than 100 invalid probe ids are detected. Rownames should be probe ids e.g. cg03607359.")} else {
      message(paste(c("Invalid probes at lines: ", inv.probes), collapse=" "))}
      stop(paste("Rownames of ltop.obj should be",platform,"probe ids."))}
      lobj <- data.frame("score"=ltop.obj[,colnames(ltop.obj)==sortBy])
      rownames(lobj)=rownames(ltop.obj)
  res <- merge(lobj, platform.probes, by="row.names",all.x=TRUE)
  res <- reorderBed(res,3,4,5,2,6)
  return(res)
  
}



#'Fetch sequences from bed-formatted data frames
#' @param bedline  Valid bed-formatted data frame
#' @param offset Number of nucleotides expanded in each direction ([1,1000])
#' @return A DNAStringset containing the sequences in hg19 genome assembly
#' @export
#'
bed2Seq <- function(bedline, offset) {
  #Depends on Biostrings, BSgenome.Hsapiens.UCSC.hg19
  if (!(dplyr::between(offset, 0, 1000)))
    stop("Offset should be between 1 and 1000.")
  if (sum(is.na(bedline)) > 0) {
    stop("Input dataset contains ", sum(is.na(bedline)), "missing values.")
  }
  strands = c("+", "-")
  if (!exists("bedline$strand")) {
    bedline$strand = "+"
  }
  seq <- getSeq(
    Hsapiens,
    bedline$chr ,
    start = bedline$start - offset,
    end = bedline$end + offset,
    strand = bedline$strand
  )
  seq@ranges@NAMES = paste0(bedline$chr,
                            sep = "_",
                            bedline$start,
                            sep = "_",
                            bedline$end)
  return(seq)
}


#'Set a study name
#'
#'Sets a name to the analysis that appears in the exported plots.
#'
#' @param study.name  A string corresponding to the name of the study
#' @return The name of the study
#' @export
nameStudy <- local(function (study.name) {
  NAME <<- study.name
  return(NAME)
}, envir = pkg.env)


#'Check if data frame is empty
#'
#'Checks if a data frame has no values
#'
#' @param df  The input data frame
#' @return TRUE/FALSE (TRUE is the data frame is empty)
#' @export
isEmptyDF <- function(df) {
  flag <- FALSE
  if (nrow(df) == 0) {
    flag <- TRUE
  }
  return(flag)
}


#' Reorder tabular methylation data to bed format
#'
#' Reorders tabular methylation data to bed-formatted files.
#' Compatible inputs are .txt, .csv data and other textual formats that
#' contain the following mandatory columns: chr, start, end and score.
#'
#'
#' @param input.data A data frame containing input bed-formatted data
#' @param chr.col Column number containing the chromosome name
#' @param start.col Column number containing the chromosome's start position
#' @param end.col Column number containing the chromosome's end position
#' @param score.col Column number containing the methylation score values either beta or delta-beta
#' @param strand.col Column number containing the strand in the use data file ('+' strand is assumed if strand column is missing)
#' @return A valid bed-formatted file (input of the `MeinteR::find*` functions)
#' @export
reorderBed <-
  function(input.data,
           chr.col,
           start.col,
           end.col,
           score.col,
           strand.col = NULL) {
    if (!is.defined(strand.col)) {
      input.data$strand = "+"
      bed <-
        cbind(input.data[, c(chr.col, start.col, end.col, score.col)], input.data$strand)
    } else {
      bed <-
        cbind(input.data[, c(chr.col, start.col, end.col, score.col, strand.col)])
    }
    colnames(bed) = c("chr", "start", "end", "score", "strand")
    bed$start <- as.numeric(as.character(bed$start))
    bed$end <- as.numeric(as.character(bed$end))
    bed$score <- as.numeric(as.character(bed$score))
    return (bed)
  }

#' Filter by promoters
#'
#' Selects genomic coordinates included in promoters based on the UCSC hg19 gene coordinates.
#'
#' @param input.data  A data frame containing input data in bed format
#' @param up.tss  Number of nucleotides upstream transcription start site
#' @param down.tss  Number of nucleotides downstream transcription start site
#' @return A data frame with the CpG sites located in promoter regions

filterByProm <- function(input.data, up.tss, down.tss) {
  if (sum(is.na(input.data)) > 0) {
    stop("Input dataset contains ", sum(is.na(input.data)), " missing values.")
  }
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  genes <- transcriptsBy(txdb, "gene")
  proms <-
    promoters(genes, upstream = up.tss, downstream = down.tss) #IRanges object
  input.GRanges <-
    makeGRangesFromDataFrame(input.data) #GRanges object
  x <-
    as.data.frame(subsetByOverlaps(input.GRanges, proms, type = "within"))
  res <-
    merge(
      input.data,
      x,
      by.x = c("chr", "start", "end", "strand"),
      by.y = c("seqnames", "start", "end", "strand")
    )
  res <- reorderBed(res, 1, 2, 3, 5, 4)
  return(res)
}

#' Filter by CpG islands
#'
#' Selects genomic coordinates included in CpG islands, using the cpgIslands dataset.
#'
#' @param input.data  A data frame containing input data in bed format
#' @return A data frame with the CpG sites located in CpG islands


filterByCGI <- function(input.data) {
  if (sum(is.na(input.data)) > 0) {
    stop("Input dataset contains ", sum(is.na(input.data)), " missing values.")
  }
  CGI <- makeGRangesFromDataFrame(cpgIslands)
  input.GRanges <-
    makeGRangesFromDataFrame(input.data) #GRanges object
  x <-
    as.data.frame(subsetByOverlaps(input.GRanges, CGI, type = "within"))
  res <-
    merge(
      input.data,
      x,
      by.x = c("chr", "start", "end", "strand"),
      by.y = c("seqnames", "start", "end", "strand")
    )
  res <- reorderBed(res, 1, 2, 3, 5, 4)
  return(res)
}


is.defined <- function(x)
  ! is.null(x)


#'Validate format of the input bed data
#'
#'Validates input methylation data. Checks the presence of the chr, start, end, score columns.
#'If column `strand` is not set then `+` strand is assumed.
#'Cleans rows with empty cells and sets numeric format to the start, end and score columns.
#'
#' @param bed.data  A data frame containing input bed-formatted data
#' @param omit.na Omit rows with empty cells (default:TRUE)
#' @return A well-formatted data frame
#' @export
#'
validateBed <- function(bed.data, omit.na = TRUE) {
  if (isEmptyDF(bed.data)) {
    stop("Empty data frame")
  }
  valid.cols <- c("chr", "start", "end", "score", "strand")
  if (!all(colnames(bed.data) %in% (valid.cols))) {
    stop("Invalid column names. Use reorderBed to match bed data with the column names.")
  }
  strands = c("+", "-")
  if (is.null(bed.data$strand)) {
    bed.data$strand <- "+"
  }
  if (isFALSE(bed.data$strand %in% (strands))) {
    bed.data$strand <- "+"
  }
  if (omit.na) {
    bed.data <- na.omit(bed.data)
  }
  if (unique(tolower(substring(bed.data$chr, 1, 3)) != "chr")) {
    bed.data$chr <- paste0("chr", bed.data$chr)
  }
  bed.data$start <- as.numeric(bed.data$start)
  bed.data$end <- as.numeric(bed.data$end)
  bed.data$score <- as.numeric(bed.data$score)
  return(bed.data)
}



#'Plot G+C-content and observed/expected ratio.
#
#'
#' Generates density plots of the G+C content and
#' observed/expected CpG ratio for the input dataset and the human genome CpG islands
#'
#' @param bed.data  A data frame containing input bed-formatted data
#' @param offset  Number of nucleotides expanded in each direction (default:200, min:20 max:1000)
#' @return 1/ A data frame containing the G+C content (percentage of island that is C or G)
#' and ratio of observed (CpG number) to expected(Number of C* Number of G/sequence length)
#' @return 2/ Density plot of the G+C-content
#' @return 3/ Density plot of the observed/expected ratio
#' @export

plotCpG <- function(bed.data, offset = 200) {
  if (!(dplyr::between(offset, 20, 1000)))
    stop("Check offset value (valid range [20,1000]).")
  while (!is.null(dev.list()))
    dev.off()
  seqs = bed2Seq(bed.data, offset)
  res = alphabetFrequency(seqs)
  L = rowSums(res)
  cprop = res[, "C"]
  gprop = res[, "G"]
  GC = (cprop + gprop) / L
  denom = cprop * gprop
  nom = vcountPattern("CG", seqs) * L
  Observed = nom / denom
  median(Observed)
  cpgIslands$obsExp <- as.numeric(as.character(cpgIslands$obsExp))
  cpgIslands$perGc <- as.numeric(as.character(cpgIslands$perGc))
  d.gc.ref <-
    density(cpgIslands$perGc / 100, na.rm = TRUE, kernel = "gaussian")
  d.gc <- density(GC, na.rm = TRUE, kernel = "gaussian")
  y.max <- max(max(d.gc$y), max(d.gc.ref$y))
  plot(
    d.gc.ref,
    xlab = "G+C-content",
    main = "Density of the G+C-content",
    xlim = c(0, 1),
    ylim = c(0, y.max),
    sub = paste0("OFFSET:", offset, "nt, N:", nrow(bed.data))
  )
  grid (NULL, NULL, lty = 6, col = "lightgrey")
  polygon(d.gc.ref, col = 'gold', border = "black")
  polygon(d.gc, col = 'lightskyblue4')
  abline(
    v = 0.5,
    b = 0,
    col = "darkgrey",
    lwd = 3,
    lty = 2
  )
  legend(
    "topright",
    legend = c("CpG islands (hg19)", local(NAME, envir = pkg.env)),
    fill = c('gold', 'lightskyblue4')
  )
  d.obs.ref <-
    density(cpgIslands$obsExp, na.rm = TRUE, kernel = "gaussian")
  d.obs <- density(Observed, na.rm = TRUE, kernel = "gaussian")
  y.max <- max(max(d.obs$y), max(d.obs.ref$y))
  plot(
    d.obs.ref,
    xlab = "Observed/Expected CpG ratio",
    main = "Density of the observed/expected CpG ratios",
    xlim = c(0, 2),
    ylim = c(0, y.max),
    sub = paste0("OFFSET: ", offset, "nt, N: ", nrow(bed.data))
  )
  grid (NULL, NULL, lty = 6, col = "lightgrey")
  polygon(d.obs.ref, col = 'gold', border = "black")
  polygon(d.obs, col = 'lightskyblue4')
  abline(
    v = 0.6,
    b = 0,
    col = "darkgrey",
    lwd = 3,
    lty = 2
  )
  legend(
    "topright",
    legend = c("CpG islands (hg19)", local(NAME, envir = pkg.env)),
    fill = c('gold', 'lightskyblue4')
  )
  res <- data.frame(GC, Observed)
  return(res)
}


#'Plot scores of the input data
#
#' Generates a density plot of the score values listed in the input dataset.
#'
#' @param bed.data A data frame containing input bed-formatted data
#' @export
plotBeta <- function(bed.data) {
  if ("score" %in% colnames(bed.data)) {
    if (!is.numeric(bed.data$score)) {
      stop("Column score must be numeric.")
    }
    if (sum(is.na(bed.data)) > 0) {
      stop("Input dataset contains ", sum(is.na(bed.data)), " missing values.")
    }
    dens <- density(bed.data$score, na.rm = TRUE, kernel = "gaussian")
    plot(dens,
         xlab = "Methylation level",
         main = "Distribution of the methylation values")
    grid (NULL, NULL, lty = 6, col = "lightgrey")
    polygon(dens, col = 'lightskyblue4')
    legend("topright",
           legend =  local(NAME, envir = pkg.env),
           fill = 'lightskyblue4')
  } else {
    stop("Score column does not exist in the input data.")
  }
}



#' Import GEO data series in the workspace
#'
#' Imports GEO data series. The function fetches data matrices corresponding to a pre-defined GSE identifier
#' and builds valid, bed-formatted dataset with delta-beta values between two sample groups, described
#' in a user-defined annotation file.
#'
#' @param gse.acc   A string corresponding to the accession number of the GEO data series
#' @param annotation.file A string corresponding to the full local path to the annotation files containing
#' sample grouping information
#' @return 1/ A bed-formatted data frame with the chromosomal coordinates and of each methylation probe and
#' the corresponding delta-beta values between the two groups
#' @return 2/ Beta values of each sample listed in the annotation file
#' @return 3/ Annotation data frame
#' @return 4/ Mean beta values of group 1
#' @return 5/ Mean beta values of group 2
#' @name importGEO
#' @export
importGEO <- function(gse.acc, annotation.file) {
  if (!is.defined(gse.acc)) {
    stop("FAIL: Data series accession number is missing.")
  }
  if (!is.defined(annotation.file)) {
    stop("Sample annotation is missing.")
  }
  if (grepl("^GSE\\d\\d+$", gse.acc)) {
    message("PASS: Valid accession number for the data series.")
  } else
  {
    stop("Invalid accession number. Currently, GSE accession numbers are supported.")
  }
  anno.df <- loadAnnotationGEO(annotation.file)
  gse <- getGEO(gse.acc, GSEMatrix = TRUE, getGPL = TRUE)
  if (length(gse) > 1) {
    stop(
      "MeinteR currently supports single-matrix data series. Data series ",
      gse.acc,
      " has ",
      length(gse),
      "."
    )
  }
  if (any(is.na(gse[[1]]@featureData@data$ID))) {
    stop("Missing probe identifiers. Retry to download data series.")
  }
  valid <- validateGEO(gse)
  if (valid) {
    "Data series looks good."
  } else
  {
    stop("Failed to import data series in the workspace.")
  }
  if (!all(anno.df$sample %in% gse[[1]]$geo_accession)) {
    stop(
      "Inconsistency between the accession numbers of the annotation file and the GEO data series.
      Check correspondence of the GSM IDs and non-visible characters in the annotation file."
    )
  }
  message("Data series: ", length(gse[[1]]$geo_accession), " samples.")
  message("Annotation file: ", nrow(anno.df), " samples.")
  all.betas <- as.data.frame(exprs(gse[[1]]))
  betas <- all.betas[, anno.df[, 1]]
  if (min(betas, na.rm = TRUE) < 0 |
      max(betas, na.rm = TRUE) > 1) {
    message("M-values are detected and will be transformed to beta values.")
    betas <- mval2beta(betas)
  }
  grp <- groupGEO(anno.df, betas)
  probe.d.beta <- grp[[1]]
  res <- probesChr(gse, probe.d.beta)
  result = list()
  result[[1]] <- res
  result[[2]] <- betas
  result[[3]] <- anno.df
  result[[4]] <- grp[[2]]
  result[[5]] <- grp[[3]]
  return(result)
  }

#Internal function for loading annotation files
loadAnnotationGEO <- function(file.path) {
  if (!file.exists(file.path)) {
    stop(paste0("File ", file.path, " not found. Check full path."))
  }
  valid.cols <- c("SAMPLE", "STATUS")
  anno <-
    read.csv(file = file.path,
             header = TRUE,
             stringsAsFactors = FALSE)
  #Remove leading and/or trailing whitespace from character strings.
  anno <-
    data.frame(lapply(anno, trimws), stringsAsFactors = FALSE)
  if (all(toupper(colnames(anno)) %in% valid.cols)) {
    message("PASS: Columns SAMPLE and STATUS are present.")
  } else {
    stop("Column names <samples> and <status> are not detected in the first line.")
  }
  if (length(colnames(anno)) != 2) {
    stop("Annotation file columns: sample identifier (GSM accession number), sample status")
  }
  if (length(unique(anno$status)) != 2) {
    stop(
      "Other than two conditions are detected. Check non-visible characters in the annotation file."
    )
  }
  return(anno)
}


validateGEO <- function(gse.list) {
  valid = TRUE
  #Validate species
  SPECIES = "9606"
  if (table(gse.list[[1]]$taxid_ch1) == length(gse.list[[1]]$geo_accession)) {
    if (unique(gse.list[[1]]$taxid_ch1) == SPECIES) {
      message("PASS: Human samples detected.")
    } else {
      message("FAIL: The samples are not human")
      valid = FALSE
    }
  } else {
    message("FAIL: Data series contains samples from more than one species.")
    valid = FALSE
  }
  
  # Validate beadchip array
  PLATFORMS <-
    c(
      "GPL13534" = "HM450k",
      "GPL16304" = "HM450k",
      "GPL21145" = "EPIC",
      "GPL8490" = "HM27k"
    )
  if (unique(gse.list[[1]]$platform_id) %in% names(PLATFORMS)) {
    message("PASS: Valid platform")
  } else {
    message(
      "FAIL: Only BeadChip methylation arrays 27k, 450k, EPIC and Illumina HiSeq 2500  are currently supported."
    )
    valid = FALSE
  }
  if (table(gse.list[[1]]$platform_id) == length(gse.list[[1]]$geo_accession)) {
    message("PASS: A unique platform ",
            PLATFORMS[gse.list[[1]]@annotation],
            " detected in the data series.")
  } else {
    message("FAIL: More than one platforms have been identified in the data series.")
    valid = FALSE
  }
  return(valid)
}

groupGEO <- function(anno.df, beta.df) {
  message("Grouping samples...")
  options(scipen = 999)
  #class.g1 if for tumor
  if (toupper(anno.df$status[1])=="NORMAL") {
    class.g1 <- subset(anno.df, anno.df$status == unique(anno.df$status)[2])
    class.g2 <- subset(anno.df, anno.df$status == unique(anno.df$status)[1])
  } else {
  class.g1 <- subset(anno.df, anno.df$status == unique(anno.df$status)[1])
  class.g2 <-subset(anno.df, anno.df$status == unique(anno.df$status)[2])
  }
  beta.g1 <- as.data.frame(beta.df[, class.g1[, 1]])
  beta.g2 <- as.data.frame(beta.df[, class.g2[, 1]])
  d.beta <-
    as.data.frame(rowMeans(beta.g1, na.rm = TRUE) - rowMeans(beta.g2, na.rm =
                                                               TRUE))
  colnames(d.beta) <-
    paste0(unique(class.g1$status), "-", unique(class.g2$status))
  db <- cbind(rownames(d.beta), d.beta)
  colnames(db)[1] <- "probes"
  grp1 <- as.data.frame(rowMeans(beta.g1, na.rm = TRUE))
  grp1$probes <- rownames(grp1)
  rownames(grp1) <- NULL
  colnames(grp1) <- c("tumor", "probes")
  grp2 <- as.data.frame(rowMeans(beta.g2, na.rm = TRUE))
  grp2$probes <- rownames(grp2)
  rownames(grp2) <- NULL
  colnames(grp2) <- c("normal", "probes")
  result <- list()
  result[[1]] <- db
  result[[2]] <- grp1
  result[[3]] <- grp2
  return(result)
}

probesChr <- function(gse, probes.db) {
  if (gse[[1]]@annotation == "GPL8490") {
    hm27 <- as.data.frame(get27k())
    platform.probes <- data.frame(rownames(hm27), hm27[, 1:3])
  } else {
    platform.probes <- cbind(
      gse[[1]]@featureData@data$ID,
      paste0("chr", gse[[1]]@featureData@data$CHR),
      gse[[1]]@featureData@data$MAPINFO,
      gse[[1]]@featureData@data$MAPINFO + 1
    )
  }
  if (ncol(platform.probes) != 4) {
    stop("Data series contains missing data. Exiting...")
  }
  colnames(platform.probes) <- c("probes", "chr", "start", "end")
  bed.data <- merge(probes.db, platform.probes, by = "probes")
  result <- bed.data
  return(result)
}

#' Reformat methylation sequencing data fetched from GEO
#'
#' Transforms sequencing data into bed-formatted files. Valid for per sample usage.
#'
#' @param file.path Local folder of the bed.gz file
#' @param cov Minimum read coverage (default:10)
#' @param chroms A vector of chromosome vector to be included in the analysis (default:ALL)
#' @return A valid bed-formatted data frame
#' @name loadSeqGEO
#' @export
loadSeqGEO <- function(file.path,
                       cov = 10,
                       chroms = NULL) {
  if (!(dplyr::between(cov, 1, 100)))
    stop("Check offset value (valid range [1,100]).")
  file.desc <- gzfile(file.path, open = "r")
  suppressWarnings(gsm <-
                     read.csv(
                       file = file.desc,
                       sep = "\t",
                       header = FALSE
                     ))
  colnames(gsm) = c("chr", "start", "end", "depth", "score", "strand")
  if (is.null(chroms)) {
    chroms = paste0("chr", c(1:22, "X", "Y"))
  }
  f.gsm <- gsm[gsm$chr %in% chroms, ]
  f.gsm$depth <- gsub('[0-9]+/', '', f.gsm$depth)
  f.gsm$depth <- gsub("\\s*\\([^\\)]+\\)", "", f.gsm$depth)
  f.gsm$depth <- as.integer(gsub("\\'", "", f.gsm$depth))
  f.gsm <- f.gsm[f.gsm$depth >= cov, ]
  f.gsm$score <- f.gsm$score / 1000
  f.gsm <- reorderBed(f.gsm, 1, 2, 3, 5, 6)
  close(file.desc)
  return(f.gsm)
}


mval2beta <- function(mval) {
  tr <- 2 ^ (mval)
  beta <- tr / (tr + 1)
  return(beta)
}
