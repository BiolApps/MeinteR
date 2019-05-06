\name{test.data}
\alias{test.data}
\docType{data}
\title{Test dataset with chromosomal position of the methylated sites and the corresponding delta beta values (sample data)
}
\description{
A valid and well-formatted sample dataset containing 401 differentially methylated data with |delta beta-values|>0.3.
}
\format{
  The dataset has the following 5 variables:
  \describe{
    \item{\code{chr}}{A factor with valid values chr1 to chr22, chrX, chrY and chrM}
    \item{\code{start}}{Start position of the methylation site}
    \item{\code{end}}{End position of the methylation site}
    \item{\code{score}}{A numeric vector with the (group1 - group2) methylation value for each methylation site.}
    \item{\code{strand}}{Strand of the methylation site, either \code{+} or \code{-}}
  }
}
\examples{
#Distribution of the methylation values
plot(test.data$score)
}
\keyword{datasets}
