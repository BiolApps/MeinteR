\name{cpgIslands}
\alias{cpgIslands}
\docType{data}
\title{CpG islands (data)
}
\description{
List of 27,718 human CpG islands with their corresponding GC-content and observed/expected ratio in chromosomes chr1..22,X,Y (hg19). The list is obtained from the cpgIslandExt table of UCSC Table Browser. 
}

\format{
    \describe{
    \item{\code{chrom}}{CpG island chromosome}
    \item{\code{chromStart}}{CpG island chromosome start position}
    \item{\code{chromEnd}}{CpG island chromosome end position}
    \item{\code{perGC}}{GC-content i.e. Percentage of the CpG island that is C or G}
    \item{\code{obsExp}}{Observed/expected ratio i.e.: Number of CpG * N / (Number of C * Number of G), N = Sequence length}
}
}
\keyword{datasets}
