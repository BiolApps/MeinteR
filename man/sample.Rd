\name{sample}
\alias{sample}
\docType{data}
\title{Sample DNA methylation dataset (sample data)
}
\description{
A data frame with 5840 methylation sites containing the chromosomal position of the methylation sites and the corresponding delta beta values.
}
\format{
  The dataset has the following 5 variables:
  \describe{
    \item{\code{Chromosome}}{A factor with valid values chr1 to chr22, chrX, chrY, chrM}
    \item{\code{Start}}{Start position of the methylation site}
    \item{\code{End}}{End position of the methylation site}
    \item{\code{Strand}}{Strand of the methylation site, either \code{+} or \code{-}}
    \item{\code{Differences}}{A numeric vector with the (group1 - group2) methylation value for each CpG site.}
  }
}
\examples{
#Distribution of the methylation values
plot(sample[,5])
}
\keyword{datasets}
