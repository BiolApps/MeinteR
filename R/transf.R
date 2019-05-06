#############################################################################
## Package    : MeinteR
## File       : transf.R
## Functions  : findConservedTFBS (exported)
##            : findTFBS (exported)
##            : scatterConsTF (exported)
##            : plotTF (exported)
##            : meinter (exported)
## Updated    : 30-04-2019
##
## Title      : Find transcription factor binding sites that overlap differentially methylated sites
#############################################################################


#' Find differentially methylated sites overlapping human/mouse/rat conserved transcription factor binding sites.
#'
#' Detects transcription factor binding sites that are conserved in human/mouse/rat alignments and overlap with the input data.
#' A binding site is considered to be conserved across the alignment if its score meets the threshold score for its binding
#' matrix in all three species. The score and threshold are computed with the Transfac Matrix
#' Database (v7.0) created by Biobase. The data are purely computational, and as such not all binding sites listed here are
#' biologically functional binding sites.
#'
#' @param bed.data   A data frame containing input bed-formatted data
#' @param known.conserved.tfbs.file (optional) Full local path to the UCSC conserved transcription factor binding sites.
#' If the table is not available locally then the script will fetch it from UCSC (needs Internet connection).
#' NOTE: It is recommended to download the compressed file (Unzipped file >290MB)
#' @return 1/ Data frame containing overlaps between \code{bed.data} and conserved transcription factor binding sites
#' @return 2/ Frequency table of conserved transcription factors on human genome (input to \code{scatterConsTF} function)
#' @return 3/ A data frame with the number of conserved transcription factor binding sites per sequence (input to \code{meinter} function)
#' @export
findConservedTFBS <-
  function(bed.data, known.conserved.tfbs.file = NULL) {
    if (is.data.frame(bed.data) && nrow(bed.data) == 0) {
      stop("The input dataframe with the differentially methylated sites is not set.")
    }
    if (sum(is.na(bed.data)) > 0) {
      stop("Input dataset contains ", sum(is.na(bed.data)), " missing values.")
    }
    if (!is.defined(known.conserved.tfbs.file)) {
      #Download or load a local copy of the conserved transcription factors
      mySession <-
        browserSession("UCSC", url = "http://genome-euro.ucsc.edu/cgi-bin/")
      genome(mySession) <- "hg19"
      message("Fetching tfbsConsSites table from UCSC... ")
      message("Downloading ~71MB gzipped data. It may take few minutes...")
      known.conserved.tfbs <- getTable(ucscTableQuery(mySession,
                                                      table = "tfbsConsSites"))
    } else {
      message("Loading local copy of UCSC/tfbsConsSites table")
      #Get the extension of the filename
      len <-
        length(unlist(strsplit(known.conserved.tfbs.file, "\\.")[[1]]))
      if (strsplit(known.conserved.tfbs.file, "\\.")[[1]][len] == "gz") {
        # Check if file is gzipped
        known.conserved.tfbs <-
          read.csv(gzfile(known.conserved.tfbs.file),
                   header = TRUE,
                   sep = "\t") #If yes unzip it
      } else {
        known.conserved.tfbs <-
          read.csv(known.conserved.tfbs.file,
                   header = TRUE,
                   sep = "\t")
      }
      colnames(known.conserved.tfbs) = c("bin",
                                         "chrom",
                                         "chromStart",
                                         "chromEnd",
                                         "name",
                                         "score",
                                         "strand",
                                         "zScore")
    }
    
    #Select columns chrom, chromStart,chromEnd,strand, name, zScore
    known.factors <- cTF
    known.tfbs <-
      known.conserved.tfbs[, c("chrom", "chromStart", "chromEnd", "strand", "name", "zScore")]
    g.bed <-
      with(
        bed.data,
        GRanges(
          bed.data$chr,
          IRanges(start = bed.data$start, end = bed.data$end),
          strand = bed.data$strand,
          score = bed.data$score
        )
      )
    tfbs <-
      with(
        known.tfbs,
        GRanges(
          known.tfbs$chrom,
          IRanges(start = known.tfbs$chromStart, end = known.tfbs$chromEnd),
          name = known.tfbs$name,
          score = known.tfbs$zScore
        )
      )
    #Identify overlaps with conserved tfbs
    g.hits = findOverlaps(g.bed, tfbs, ignore.strand = FALSE, select = "all")
    g.df <- cbind(as.data.frame(g.bed[queryHits(g.hits)]),
                  as.data.frame(tfbs[subjectHits(g.hits)]))
    g.df <- g.df[, !duplicated(colnames(g.df))]
    if (isEmptyDF(g.df)) {
      stop("No overlaps between DMS and conserved transcription factor binding sits are found")
    }
    g.factor <-
      merge(g.df, known.factors, by = "name", all.x = TRUE)
    
    names(g.factor)[names(g.factor) == "seqnames"] <- "chr"
    g.freq <- count(g.factor, 'factor')
    m.freq <- merge(g.freq, refFreq, by = "factor")
    colnames(m.freq) <- c("Factor", "Count", "Ref.freq")
    count(g.factor$chr)
    g.factor$coor <-
      paste(g.factor$chr, g.factor$start, g.factor$end, sep = "_")
    freq <- data.frame(count(g.factor$coor))
    colnames(freq) <- c("seqnames", "freq")
    bed.coor <-
      data.frame(paste(bed.data$chr, bed.data$start, bed.data$end, sep = "_"))
    colnames(bed.coor) <- "seqnames"
    meinter <- merge(bed.coor, freq, all.x = TRUE)
    meinter$freq[is.na(meinter$freq)] <- 0
    result <- list()
    result[[1]] <- g.factor
    result[[2]] <- m.freq
    result[[3]] <- meinter
    return(result)
  }

#' Find putative transcription factor binding sites
#'
#' Detects JASPAR's transcription factor binding sites (core collection), co-localized with input data. Both sequence strands
#' are examined. The analysis can be restricted to promoters (use `uptss` and `down.tss` to define promoter length, relative
#' to transcription start site) and CpG islands of the human genome (hg19).
#'
#' @param bed.data  A data frame containing input bed-formatted data
#' @param persim Minimum similarity with transcription factors consensus matrices (default:0.8, range in [0,1])
#' @param offset Number of nucleotides expanded in each direction (default:12, min:5, max:100)
#' @param target Search for transcription factor binding sites on specific regions.
#' `PROMOTER`: selects sites located in promoter regions, `CGI`: selects sites in CpG islands,
#' `ALL`: No filtering is applied (time-consuming for large datasets) (default: "PROMOTER")
#' @param up.tss  Number of nucleotides upstream transcription
#' start site (Only when target="PROMOTER" is set, default: 1000)
#' @param down.tss  Number of nucleotides downstream transcription
#' start site (Only when target="PROMOTER" is set, default: 100)
#' @param mcores Number of cores to be used (default: maximum available)
#' @param tf.ID A vector of JASPAR transcription factors identifiers to search for (default: all)
#' @return 1/ Data frame containing the transcription factors identified in each sequence, their
#' position and binding score (input to `plotTF` function)
#' @return 2/ Data frame of the detected transcription factor binding sites per sequence (input to `meinter` function)
#' @export



findTFBS <-
  function(bed.data,
           persim = 0.8,
           offset = 12,
           target = "PROMOTER",
           up.tss = 1000,
           down.tss = 100,
           mcores = NULL,
           tf.ID = NULL) {
    if (!(dplyr::between(offset, 5, 100)))
      stop("Check offset value (valid range [5,100]).")
    options = c("ALL", "PROMOTER", "CGI")
    DNA.ALPHABET <- c("A", "C", "G", "T")
    JASPAR.SPECIES <- 9606 # 9606:  Human genome support only
    if (!(dplyr::between(persim, 0, 1)))
      stop("Check percent of similarity. Valid range [0,1]")
    
    
    if (sum(is.na(bed.data)) > 0) {
      stop("Input dataset contains ", sum(is.na(bed.data)), " missing values.")
    }
    if (toupper(target) == "PROMOTER") {
      bed.data <-
        filterByProm(bed.data, up.tss = up.tss, down.tss = down.tss)
    }
    if (target == "CGI") {
      bed.data <- filterByCGI(bed.data)
    }
    if (!(toupper(target) %in% options)) {
      stop(
        "Invalid target value. Select among: search all
        dataset (ALL), promoters (PROMOTERS) or CpG islands (CGI) "
      )
    }
    subject <- bed2Seq(bed.data, offset)
    opts = list()
    opts[["species"]] = JASPAR.SPECIES
    opts[["all_versions"]] = FALSE
    opts[["collection"]] = "CORE"
    if (is.defined(tf.ID)) {
      opts[["ID"]] = tf.ID
    }
    PWMatrixList = TFBSTools::toPWM(TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts))
    message("No. of sequences: ", nrow(bed.data), " (", target, ")")
    message("Number of transcription factors: ", length(PWMatrixList))
    message("Start analysing each site ...")
    message(
      "Searching hundreds of trascription factors in thousands of sequences will take several minutes to hours."
    )
    
    if (!is.defined(mcores)) {
      cores <- detectCores()
    } else {
      cores = mcores
    }
    message(paste(cores, "core(s) will be used for this analysis."))
    sitesetL <-
      TFBSTools::searchSeq(
        PWMatrixList,
        subject,
        seqname = "seq1",
        min.score = persim,
        strand = "*",
        mc.cores = cores
      )
    tfbs.frame <- as(sitesetL, "data.frame")
    tfbs.freq <- plyr::count(tfbs.frame, "seqnames")
    tfbs.freq.sep <-
      data.frame(do.call(rbind, strsplit(
        as.character(tfbs.freq$seqnames),
        split = "_",
        fixed = TRUE
      )), tfbs.freq$freq)
    colnames(tfbs.freq.sep) = c("chr", "start", "end", "tfbs")
    tfbs.freq.sep$start <-
      as.integer(as.character(tfbs.freq.sep$start))
    tfbs.freq.sep$end <- as.integer(as.character(tfbs.freq.sep$end))
    suppressMessages(meinter <-
                       join(bed.data[, 1:3], tfbs.freq.sep[, c("chr", "start", "end", "tfbs")]))
    meinter$tfbs[is.na(meinter$tfbs)] <- 0
    result <- list()
    result[[1]] <- tfbs.frame
    result[[2]] <- meinter
    return(result)
    }



#' Create barplot of the identified transcription factor binding sites
#'
#' Generates an overlayed barplot of the results exported by the `findTFBS` function.
#' The bar plot visualises the most frequent transription factors with respect
#' to the total number of occurrences and the number of sequences that contain
#' these trascription factors.
#'
#' @param df   The data frame exported by the `findTFBS` function
#' @param topTF Integer corresponding to the  number of the most frequent trascription
#' factors to be displayed (default:10)
#' @return 1/ A barplot with the `topTF` most frequent transcription factors
#' @return 2/ A barplot with the number of transcription factors per class
#' @return 3/ A scatterplot comparing the observed and expected number of transcription
#' factors per class
#' @export

plotTF <- function(df, topTF = 10) {
  if (topTF < 1)
    stop("topTF > 0")
  if (isEmptyDF(df)) {
    stop('Data frame is empty')
  }
  freq.TF <- as.data.frame(table(df$TF), decreasing = TRUE)
  TF.seqs <- as.data.frame(table(df$TF, df$seqnames))
  TF.seqs.c <-
    as.data.frame(rowSums(table(TF.seqs$Var1, TF.seqs$Freq)[, -1]))
  TF.seqs.c$TF <- rownames(TF.seqs.c)
  colnames(TF.seqs.c) = c("seq.num", "TF")
  colnames(freq.TF) = c("TF", "TF.freq")
  df.merge <- merge(freq.TF, TF.seqs.c, by = "TF")
  df.merge <-
    df.merge[order(-df.merge$TF.freq,-df.merge$seq.num),]
  top <- head(df.merge, topTF)
  if (min(top$seq.num < 1 / 500 * top$TF.freq)) {
    message("The number of sequences is
            too low to appear in the barchart")
  }
  cols <- c("plum1", "plum4")
  p <- ggplot(top, aes(x, y)) +
    geom_bar(
      stat = "identity",
      aes(
        x = top$TF,
        y = top$TF.freq,
        fill = "Total TFBS"
      ),
      width = 0.8,
      colour = "grey",
      alpha = 0.5,
      position = "identity"
    ) +
    geom_bar(
      stat = "identity",
      aes(
        x = top$TF,
        y = top$seq.num,
        fill = "Sequences with TFBS"
      ),
      width = 0.6,
      colour = "grey",
      alpha = 0.5,
      position = "identity"
    ) +
    scale_x_discrete(name = "Transcription factors") +
    scale_y_discrete(
      limits = c(0, top$TF.freq),
      name = "Number of
      transcription factor binding sites"
    ) +
    geom_text(aes(
      x = top$TF,
      y = top$TF.freq,
      label = top$TF.freq
    )) +
    geom_text(aes(
      x = top$TF,
      y = top$seq.num,
      label = top$seq.num
    )) +
    ggtitle(paste0("Transcription factors (", local(NAME, envir = pkg.env), ")")) +
    scale_fill_manual(name = "Legend", values = cols) +
    theme(legend.position = "bottom") +
    theme(plot.title = element_text(
      size = 20,
      hjust = 0.5,
      lineheight = .8,
      face = "bold"
    ))
  #Plot the frequency of each TF class
  cnt = count(df$class)
  q <- ggplot(cnt, aes(x, y)) +
    geom_bar(
      stat = "identity",
      aes(x = cnt$x, y = cnt$freq),
      width = 0.8,
      colour = "grey37",
      alpha = 0.5,
      position = "identity"
    ) +
    scale_x_discrete(name = "Transcription factor classes") +
    scale_y_discrete(limits = c(0, cnt$freq),
                     name = "Number of transcription factor
                     binding sites") +
    geom_text(aes(
      x = cnt$x,
      y = cnt$freq,
      label = cnt$freq
    )) +
    ggtitle(paste0(
      "Transcription factor classes (",
      local(NAME, envir = pkg.env),
      ")"
    )) +
    theme(plot.title = element_text(
      size = 20,
      hjust = 0.5,
      lineheight = .8,
      face = "bold"
    )) + coord_flip()
  fdf <- merge(cnt, TF.class, by.x = "x", by.y = "class")
  r <- ggplot(data = fdf, aes(x = freq, Number)) +
    geom_point() +
    geom_text(
      data = subset(fdf, freq > summary(fdf$freq)[["Mean"]] |
                      Number > summary(fdf$Number)[["Mean"]]),
      aes(label = x),
      hjust = 0,
      vjust = 0
    ) +
    scale_x_continuous(name = "Number of transcription factors per class (Observed)") +
    scale_y_continuous(name = "Number of transcription factors per class (Expected)") +
    geom_smooth(
      method = lm ,
      color = "purple",
      se = TRUE,
      fill = "lightgrey"
    ) +
    ggtitle("Scatterplot of the expected/observed transcription factor classes") +
    theme(plot.title = element_text(
      size = 14,
      hjust = 0.5,
      lineheight = .4,
      face = "bold"
    )) +
    geom_rug(col = "steelblue",
             alpha = 0.1,
             size = 1.5)
  ret <- list()
  ret[[1]] <- p
  ret[[2]] <- q
  ret[[3]] <- r
  return(ret)
  }

#' Create a scatterplot of the identified conserved transcription factors
#'
#' Generates a scatterplot of the results exported by the `findConservedTFBS` function. The scatterplot
#' illustrates the number of binding sites per transription factor relative to the expected
#' frequency on the reference human genome. The trascription factors with high frequency (>= 3rd quantile)
#' in the reference genome or to the analysed data are labeled on the scatterplot.
#'
#' @param df   The data frame exported by the `findConservedTFBS` function
#' @export

scatterConsTF <- function(df) {
  ggplot(data = df, aes(x = Count, Ref.freq)) +
    geom_point() +
    geom_text(
      data = subset(
        df,
        Count > summary(df$Count)[["3rd Qu."]] |
          Ref.freq > summary(df$Ref.freq)[["3rd Qu."]]
      ),
      aes(label = Factor),
      hjust = 0,
      vjust = 0
    ) +
    scale_x_continuous(name = "Number of detected binding sites") +
    scale_y_continuous(name = "Relative frequency in human genome (hg19)") +
    geom_smooth(method = "glm" ,
                color = "red",
                se = TRUE) +
    ggtitle(paste0(
      "Conserved trascription factors (",
      local(NAME, envir = pkg.env),
      ")"
    )) +
    theme(plot.title = element_text(
      size = 20,
      hjust = 0.5,
      lineheight = .8,
      face = "bold"
    )) +
    geom_rug(col = "steelblue",
             alpha = 0.1,
             size = 1.5)
}


#' Calculate the genomic index of methylation sites based on the Meinter's `find*` functions' outputs
#'
#' Calculates the genomic index given a set of features pre-analysed using MeinteR's `find*` functions.
#' First, the function builds the local genomic signature of each site and the it calculates the genomic
#' index using a weighting scheme.


#' @param bed.data  A data frame containing input bed-formatted data
#' @param funList List of `find*` functions outputs. At least one core function is needed to calculate
#' the genomic index. Valid element names of the list: `spls`-`findSpliceSites`, `altss`-`findAltSplicing`,
#' `ctfbs`-`findConservedTFBS`, `tfbs`-`findTFBS`, `pals`-`findPals`, `quads`-`findQuads`, `shapes`-`findShapes`
#' @param weights A list of positive values corresponding to feature weights [0,10]. Same list elements with
#' `funList` list
#' @return  A data frame with the genomic index of the input data
#' @export
#'
meinter <- function(bed.data, funList, weights)
{
  if (length(funList) == 0)
    stop("funList is empty.")
  if (length(weights) == 0)
    stop("Weight list is empty.")
  P.VAL = 0.05
  MAX.RANGE = 10
  if (!exists("bed.data")) {
    stop("Parameter bed.data not defined")
  }
  mtx <-
    data.frame(paste(bed.data$chr, bed.data$start, bed.data$end, sep = "_"),
               bed.data$score)
  colnames(mtx) <- c("seqname", "score")
  v.fts <- c()
  SS <- FALSE
  if (!"spls" %in% names(funList)) {
    message("EXCLUDED FEATURE: Splice sites")
    mtx$spls = 0
  } else {
    if (!"spls" %in% names(weights) ||
        !(dplyr::between(weights[["spls"]], 0, MAX.RANGE)))
      stop("Check weight for splice sites. Valid range [0,10].")
    else
      SS <- TRUE
    v.fts <- c(v.fts, "spls")
    colnames(funList[["spls"]][[1]])[colnames(funList[["spls"]][[1]]) ==
                                       "freq"] <- "spls"
    mtx <-
      merge(x = mtx, y = funList[["spls"]][[1]], by = "seqname")
  }
  
  ALTSS <- FALSE
  if (!"altss" %in% names(funList)) {
    message("EXCLUDED FEATURE: Alternative splicing events")
    if (!SS)
      mtx$spls = NULL
  } else {
    if (!"spls" %in% names(weights) ||
        !(dplyr::between(weights[["spls"]], 0, MAX.RANGE)))
      stop("Check weight for alternative splicing.Valid range [0,10].")
    else
      ALTSS <- TRUE
    if (!("spls" %in% v.fts))
      v.fts <- c(v.fts, "spls")
    funList[["altss"]][[3]]$seqname <-
      paste(funList[["altss"]][[3]]$chr, funList[["altss"]][[3]]$start, funList[["altss"]][[3]]$end, sep =
              "_")
    colnames(funList[["altss"]][[3]])[colnames(funList[["altss"]][[3]]) ==
                                        "obs"] <- "altss"
    mtx <-
      merge(x = mtx, y = funList[["altss"]][[3]][, 4:5], by = "seqname")
    mtx$spls <- as.numeric((mtx$altss + mtx$spls) > 0)
    mtx$altss = NULL
  }
  
  if (!ALTSS && !SS) {
    mtx$spls <- NULL
  }
  
  TFBS <- FALSE
  if (!"tfbs" %in% names(funList))
    message("EXCLUDED FEATURE: Transcription factors")
  else {
    if (!"tfbs" %in% names(weights) ||
        !(dplyr::between(weights[["tfbs"]], 0, MAX.RANGE)))
      stop("Check weight for transcription factor bindings. Valid range [0,10].")
    else
      TFBS <- TRUE
    v.fts <- c(v.fts, "tfbs")
    funList[["tfbs"]][[2]]$seqname <-
      paste(funList[["tfbs"]][[2]]$chr, funList[["tfbs"]][[2]]$start, funList[["tfbs"]][[2]]$end, sep =
              "_")
    message(
      "Genomic index is available for ",
      nrow(funList[["tfbs"]][[2]]),
      " genomic loci out of ",
      nrow(bed.data),
      " sites."
    )
    mtx <-
      merge(
        y = mtx,
        x = funList[["tfbs"]][[2]][, 4:5],
        all.x = TRUE,
        by = "seqname"
      )
    mtx$tfbs <-
      sapply(mtx$tfbs, function(x) {
        round(x / max(mtx$tfbs, na.rm = TRUE), 1)
      })
  }
  
  CTFBS <- FALSE
  if (!"ctfbs" %in% names(funList))
    message("EXCLUDED FEATURE: Conserved transcription factors")
  else {
    if (!"ctfbs" %in% names(weights) ||
        !(dplyr::between(weights[["ctfbs"]], 0, MAX.RANGE)))
      stop("Check weight for conserved transcription factor bindings. Valid range [0,10].")
    else
      CTFBS <- TRUE
    v.fts <- c(v.fts, "ctfbs")
    colnames(funList[["ctfbs"]][[3]])[colnames(funList[["ctfbs"]][[3]]) ==
                                        "freq"] <- "ctfbs"
    mtx <-
      merge(
        x = mtx,
        y = funList[["ctfbs"]][[3]],
        by.x = "seqname",
        by.y = "seqnames"
      )
    mtx$ctfbs <-
      sapply(mtx$ctfbs, function(x) {
        round(x / max(mtx$ctfbs), 1)
      })
  }
  
  PALS <- FALSE
  if (!"pals" %in% names(funList))
    message("EXCLUDED FEATURE: Palindromes")
  else {
    if (!"pals" %in% names(weights) ||
        !(dplyr::between(weights[["pals"]], 0, MAX.RANGE)))
      stop("Check weight for palindromic sequences. Valid range [0,10].")
    else
      PALS <- TRUE
    v.fts <- c(v.fts, "pals")
    #  colnames(funList[["pals"]][[3]])[colnames(funList[["pals"]][[3]])=="freq"] <- "pals"
    mtx <-
      merge(x = mtx, y = funList[["pals"]][[3]], by = "seqname")
    mtx$pals <-
      sapply(mtx$pals, function(x) {
        round(x / max(mtx$pals), 1)
      })
  }
  
  QUADS <- FALSE
  if (!"quads" %in% names(funList))
    message("EXCLUDED FEATURE: G-quadruplexes")
  else {
    if (!"quads" %in% names(weights) ||
        !(dplyr::between(weights[["quads"]], 0, MAX.RANGE)))
      stop("Check weight for G-quadruplex sequences. Valid range [0,10].")
    else
      QUADS <- TRUE
    v.fts <- c(v.fts, "quads")
    mtx <-
      merge(x = mtx, y = funList[["quads"]][[2]], by = "seqname")
    mtx$quads <-
      sapply(mtx$quads, function(x) {
        round(x / max(mtx$quads), 1)
      })
  }
  
  SHAPES <- FALSE
  if (!"shapes" %in% names(funList))
    message("EXCLUDED FEATURE: Shape features")
  else {
    if (!"shapes" %in% names(weights) ||
        !(dplyr::between(weights[["shapes"]], 0, MAX.RANGE)))
      stop("Check weight for DNA conformations. Valid range [0,10].")
    else
      SHAPES <- TRUE
    v.fts <- c(v.fts, "shapes")
    tmp.shapes <-
      cbind(funList[["shapes"]][[1]], funList[["shapes"]][[2]][, 2], funList[["shapes"]][[3]][, 2], funList[["shapes"]][[4]][, 2])
    # at least one conformation statistically significant
    tmp.shapes$shapes <-
      apply(tmp.shapes[, 2:ncol(tmp.shapes)], 1, function(x) {
        sum(x < P.VAL)
      })
    mtx <- merge(x = mtx, y = tmp.shapes[, c(1, 6)], by = "seqname")
    rm(tmp.shapes)
  }
  idx <- vector()
  for (i in v.fts) {
    ptl <- mtx[[i]] * weights[[i]]
    idx <- apply(cbind(idx, ptl), 1, sum, na.rm = TRUE)
  }
  tmp <- strsplit(as.character(mtx$seqname), split = "_")
  loc <- do.call(rbind, lapply(tmp, rbind))
  res <- as.data.frame(cbind(loc, mtx$score, idx))
  colnames(res) <- c("chr", "start", "end", "score", "g.index")
  res$chr <- as.character(res$chr)
  indx <- sapply(res, is.factor)
  res[indx] <-
    lapply(res[indx], function(x)
      as.numeric(as.character(x)))
  res <- res[order(res$g.index, decreasing = TRUE),]
  return(res)
}
