#############################################################################
## Package    : MeinteR
## File       : splicing.R
## Functions  : findAltSplicing (exported)
##            : findSpliceSites (exported)
##
## Updated    : 30-04-2019
##
## Title      : Alternative splicing events and canonical exons overlapping differentially methylated sites.
#############################################################################

#'Find alternative splicing events
#
#'
#' Idetifies known alternative splicing events co-localized with input data.
#'
#' @param bed.data A data frame containing input bed-formatted data
#' @param known.alt.splice (optional) Full local path to the UCSC knownalt table.
#' If the table is not available locally then the script will fetch known alternative splicing
#' events from UCSC (needs Internet connection).
#' @return 1/ A data frame with the identified alternative splicing event overlaps (hg19)
#' @return 2/ A summary table with the frequency of each alternative splicing event compared to the
#' reference frequency
#' @return 3/ A data frame with the number of alternative splicing events per sequence (input to \code{meinter} function)
#' @return 4/ An overlayed bar chart object
#' @export
findAltSplicing <- function(bed.data, known.alt.splice = NULL) {
  options(warn = -1)
  if (!is.defined(bed.data)) {
    stop("Bed data is not defined.")
  }
  if (sum(is.na(bed.data)) > 0) {
    stop("Input dataset contains ", sum(is.na(bed.data)), " missing values.")
  }
  if (!is.defined(known.alt.splice)) {
    mySession <-
      browserSession("UCSC", url = "http://genome-euro.ucsc.edu/cgi-bin/")
    genome(mySession) <- "hg19"
    message("Fetching knownAlt table from UCSC table browser...")
    known.alt <-
      getTable(ucscTableQuery(mySession, table = "knownAlt"))
    message("DONE")
  } else {
    len <-
      length(unlist(strsplit(known.alt.splice, "\\.")[[1]])) #Get the extension of the filename
    if (strsplit(known.alt.splice, "\\.")[[1]][len] == "gz") {
      # Check if file is gzipped
      known.alt <-
        read.csv(gzfile(known.alt.splice),
                 header = TRUE,
                 sep = "\t") #If yes unzip it
    } else {
      known.alt <- read.csv(known.alt.splice, header = TRUE, sep = "\t")
    }
  }
  #Check if target cytosine is part of an alternative splicing events and report the merged
  #lines bed + altsplice record
  g.bed = with(
    bed.data,
    GRanges(
      bed.data$chr,
      IRanges(start = bed.data$start, end = bed.data$end),
      strand = bed.data$strand,
      score = bed.data$score
    )
  )
  alt = with(
    known.alt,
    GRanges(
      known.alt$chrom,
      IRanges(start = known.alt$chromStart, end = known.alt$chromEnd),
      strand = known.alt$strand,
      name = known.alt$name
    )
  )
  g.hits = findOverlaps(g.bed, alt, ignore.strand = FALSE, select = "all")
  if (length(g.hits) == 6) {
    stop("No alternative splicing events were found in the dataset.")
  }
  g.merged <-
    cbind(as.data.frame(g.bed[queryHits(g.hits)]), as.data.frame(alt[subjectHits(g.hits)]))
  g.count <- count(g.merged, 'name')
  colnames(g.count)[2] <- "Count"
  ref.count <- count(known.alt, 'name')
  colnames(ref.count)[2] <- "ref.count"
  g.count$freq <- (g.count$Count / sum(g.count$Count)) * 100
  ref.count$ref.freq <-
    (ref.count$ref.count / sum(ref.count$ref.count)) * 100
  
  total.freq <- merge(g.count, ref.count, by = "name", all = TRUE)
  total.freq[is.na(total.freq)] <- 0
  total.freq$freq <- round(total.freq$freq, 2)
  total.freq$ref.freq <- round(total.freq$ref.freq, 2)
  
  cols <- c("steelblue2", "steelblue4")
  plot.title = paste0(
    local(NAME, envir = pkg.env),
    " events: ",
    nrow(g.merged),
    ",",
    " Reference: ",
    nrow(known.alt),
    " events"
  )
  p <- ggplot(total.freq, aes(x, y)) +
    geom_bar(
      stat = "identity",
      aes(
        x = total.freq$name,
        y = total.freq$freq,
        fill = local(NAME, envir = pkg.env)
      ),
      width = 0.8,
      colour = "grey",
      alpha = 0.5,
      position = "identity"
    ) +
    geom_bar(
      stat = "identity",
      aes(
        x = total.freq$name,
        y = total.freq$ref.freq,
        fill = "Reference (hg19)"
      ),
      width = 0.6,
      colour = "grey",
      alpha = 0.5,
      position = "identity"
    ) +
    scale_x_discrete(name = "Alternative splicing events") +
    scale_y_discrete(limits = c(0, total.freq$freq), name = "Frequency (%)") +
    geom_text(aes(
      x = total.freq$name,
      y = total.freq$freq,
      label = total.freq$freq
    )) +
    geom_text(aes(
      x = total.freq$name,
      y = total.freq$ref.freq,
      label = total.freq$ref.freq
    )) +
    ggtitle(plot.title) +
    scale_fill_manual(name = "Legend", values = cols) +
    theme(legend.position = "bottom") +
    theme(plot.title = element_text(
      family = "Garamond",
      size = 20,
      hjust = 0.5,
      lineheight = .8,
      face = "bold"
    ))
  colnames(g.merged) = c(
    "chr",
    "start",
    "end",
    "width",
    "strand",
    "score",
    "alt.chr",
    "alt.start",
    "alt.end",
    "alt.width",
    "alt.strand",
    "name"
  )
  g.merged$obs <- 1
  suppressMessages(meinter <-
                     join(bed.data[, 1:3], g.merged[, c("chr", "start", "end", "obs")]))
  meinter$obs[is.na(meinter$obs)] <- 0
  result <- list()
  result[[1]] <- g.merged
  result[[2]] <- total.freq
  result[[3]] <- unique(meinter)
  result[[4]] <- p
  return(result)
}


#'Find splice sites
#'
#' Detects potential splice sites in the proximal region of the input genomic coordinates.
#' The function implements the prediction model proposed by Shapiro and Senapathy (Shapiro MB, Senapathy P.
#' Nucleic Acids Research. 1987;15(17):7155-7174.)
#'
#' @param bed.data A data frame containing input bed-formatted data
#' @param persim Similarity with the splice site consensus (default:0.8, range between [0,1])
#' @param offset Number of nucleotides expanded in each direction (default:10, min:5, max:50)
#' @return 1/ A detailed table with the location of the detected splice sites in each
#' sequence and the corresponding similarity score
#' @return 2/ A summary table with the number of splice sites detected in each sequence (input `meinter` function)
#' @export

findSpliceSites <- function(bed.data,
                            persim = 0.8,
                            offset = 10) {
  if (!(dplyr::between(persim, 0, 1)))
    stop("Check similarity level (valid range [0,1]).")
  if (!(dplyr::between(offset, 1, 50)))
    stop("Check offset value (valid range [5,50]).")
  subject <- bed2Seq(bed.data, offset)
  pfm5 <- PFMatrix(
    ID = "SS5",
    name = "Donor",
    matrixClass = "SpliceSite",
    strand = "+",
    bg = c(
      A = 0.25,
      C = 0.25,
      G = 0.25,
      T = 0.25
    ),
    tags = list(family = "ShapiroSenapathy",
                medline = "3658675"),
    profileMatrix = matrix(
      c(
        32L,
        58L,
        10L,
        0L,
        0L,
        57L,
        71L,
        5L,
        16L,
        37L,
        13L,
        4L,
        0L,
        0L,
        2L,
        8L,
        6L,
        15L,
        19L,
        15L,
        78L,
        100L,
        0L,
        39L,
        12L,
        84L,
        22L,
        12L,
        15L,
        8L,
        0L,
        100L,
        2L,
        9L,
        5L,
        47L
      ),
      byrow = TRUE,
      nrow = 4,
      dimnames = list(c("A", "C", "G", "T"))
    )
  )
  
  pfm3 <- PFMatrix(
    ID = "SS3",
    name = "Acceptor",
    matrixClass = "SpliceSite",
    strand = "+",
    bg = c(
      A = 0.25,
      C = 0.25,
      G = 0.25,
      T = 0.25
    ),
    tags = list(family = "ShapiroSenapathy", medline = "3658675"),
    profileMatrix = matrix(
      c(
        9L,
        9L,
        7L,
        7L,
        10L,
        10L,
        7L,
        9L,
        6L,
        6L,
        23L,
        3L,
        100L,
        0L,
        28L,
        31L,
        33L,
        31L,
        35L,
        35L,
        35L,
        43L,
        41L,
        39L,
        40L,
        29L,
        74L,
        0L,
        0L,
        13L,
        15L,
        13L,
        11L,
        7L,
        7L,
        11L,
        7L,
        8L,
        6L,
        8L,
        23L,
        1L,
        0L,
        100L,
        49L,
        45L,
        45L,
        51L,
        51L,
        47L,
        44L,
        42L,
        42L,
        48L,
        46L,
        24L,
        22L,
        0L,
        0L,
        10L
      ),
      byrow = TRUE,
      nrow = 4,
      dimnames = list(c("A", "C", "G", "T"))
    )
  )
  pwm5 <- toPWM(pfm5)
  pwm3 <- toPWM(pfm3)
  
  siteset5 <- searchSeq(pwm5, subject, min.score = persim, strand = "*")
  siteset3 <- searchSeq(pwm3, subject, min.score = persim, strand = "*")
  ss5 <- writeGFF3(siteset5)
  ss3 <- writeGFF3(siteset3)
  
  res.verbose <- rbind(ss5, ss3)[, c(1, 4:7, 9)]
  rownames(res.verbose) <- 1:nrow(res.verbose)
  res.sum <- plyr::count(res.verbose, "seqname")
  res.sum$seqname <- as.character(res.sum$seqname)
  subject.df <- as.data.frame(subject)
  subject.df$seqname <- row.names(subject.df)
  res.sum.f <- merge(subject.df, res.sum, by = "seqname", all.x = TRUE)
  res.sum.f$freq[is.na(res.sum.f$freq)] <- 0
  result <- list()
  result[[1]] <- res.sum.f[, -2]
  result[[2]] <- res.verbose
  return(result)
}