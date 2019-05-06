#######################################################################
## Package    : MeinteR
## File       : quadromes.R
## Functions  : findPals (exported)
##            : findQuads (exported)
## Updated    : 30-04-2019
##
## Title      : Functions for palindrome and G-quadruplex detection.
#######################################################################

#' Find palindromes in a bed-formatted dataset
#'
#' Deetects whether the target cytosine overlaps with a palindromic sequences or it is located inbetween of the two
#' arms of a palindromic sequence i.e. in the loop formed by the palindrome.
#'
#' @param bed.data A data frame containing input bed-formatted data
#' @param offset Number of nucleotides expanded in each direction (default:10, max:200)
#' @param min.arm Minimum length of each arm (default:5)
#' @param max.loop Maximum length of the loop between the two arms of the palindrome
#' @param max.mismatch The maximum number of mismatching letters
#' allowed between the two arms of the palindromes
#' @return 1/ DNAString subject with the identified palindromes
#' @return 2/ Number of palindromes falling on/neighboring input data
#' @return 3/ Number of palindromes per sequence (input to `meinter` function)
#' @export
findPals <- function(bed.data,
                     offset = 10,
                     min.arm = 5,
                     max.loop = 5,
                     max.mismatch = 1) {
  if (offset > 200 | offset < 0) {
    stop('Valid offset: [1,200]')
  }
  if (max.mismatch < 0 |
      max.mismatch > 5) {
    stop('Allowed mismatches: [0,10]')
  }
  if (sum(is.na(bed.data)) > 0) {
    stop("Input dataset contains ", sum(is.na(bed.data)), " missing values.")
  }
  seqs = bed2Seq(bed.data, offset)
  res = lapply(seqs, function(x) {
    findPalindromes(
      x,
      min.armlength = min.arm,
      max.looplength = max.loop,
      max.mismatch = max.mismatch
    )
  })
  pal.start = pal.width = arm.length = vector()
  mc.loc = offset + 1
  for (i in seq_len(length(res))) {
    pal.start <- c(pal.start, res[[i]]@ranges@start)
    pal.width <- c(pal.width, res[[i]]@ranges@width)
    arm.length <-
      c(
        arm.length,
        palindromeArmLength(
          res[[i]],
          min.armlength = min.arm,
          max.looplength = max.loop,
          max.mismatch = max.mismatch
        )
      )
  }
  l.arm = r.arm = vector()
  off.pal = on.pal = 0
  for (i in seq_len(length(pal.start))) {
    if (pal.start[i] > mc.loc) {
      off.pal = off.pal + 1
    } else if (pal.start[i] + pal.width[i] < mc.loc) {
      off.pal = off.pal + 1
    } else {
      on.pal = on.pal + 1
    }
  }
  seqname <- c()
  pals <- c()
  for (name in names(res)) {
    seqname <- c(seqname, name)
    pals <- c(pals, length(res[[name]]@ranges))
  }
  result <- list()
  result[[1]] <- res
  result[[2]] <- data.frame(on.pal, off.pal)
  result[[3]] <- data.frame(seqname, pals)
  return(result)
}

#' Find quadruplexes in sequences centered at CpG sites
#'
#' This function will detect DNA sequence patterns that likely fold into G-quadruplex structures.
#'
#' @param bed.data    A data frame containing input bed-formatted data
#' @param offset Number of nucleotides expanded in each direction (default:100, max:1000)
#' @return A DNAString subject with the identified G-quadruplexes, their length and relative coordinates
#' @return Number of G-quadruplexes per sequence (input to `meinter` function)
#' @export

findQuads <- function(bed.data, offset = 100) {
  if (offset > 1000 | offset < 0) {
    stop('Valid offset: [1..1000]')
  }
  if (sum(is.na(bed.data)) > 0) {
    stop("Input dataset contains ", sum(is.na(bed.data)), " missing values.")
  }
  message("Number of sequences to examine: ", nrow(bed.data))
  if (nrow(bed.data) > 10000 ||
      offset > 500) {
    warning("Large datasets and/or broad genomic regions might few minutes to be analysed...")
  }
  seqs = bed2Seq(bed.data, offset)
  
  res <- lapply(seqs, function(x) {
    pqsview <- pqsfinder(x)
    list(pqsview, q = pqsview@elementMetadata@nrows)
  })
  seqname <- c()
  quads <- c()
  for (name in names(res)) {
    seqname <- c(seqname, name)
    quads <- c(quads, res[[name]]$q)
  }
  result <- list()
  result[[1]] <- res
  result[[2]] <- data.frame(seqname, quads)
  return(result)
}
