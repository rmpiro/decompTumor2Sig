% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluateDecompositionQuality.R
\name{evaluateDecompositionQuality}
\alias{evaluateDecompositionQuality}
\title{Evaluate tumor decomposition quality.}
\usage{
evaluateDecompositionQuality(exposure, signatures, genome,
plot=FALSE)
}
\arguments{
\item{exposure}{(Mandatory) A single vector containing the estimated
signature contributions, or exposures, of a single tumor as provided by
\code{decomposeTumorGenomes}. The number of elements of the
exposure vector must correspond to the number of signatures (see below).}

\item{signatures}{(Mandatory) The list of signatures (vectors, data
frames or matrices) for which the exposures were obtained. Each of the
list objects represents one mutational signature. Vectors are used for
Alexandrov signatures, data frames or matrices for Shiraishi signatures.}

\item{genome}{(Mandatory) A single tumor genome in form of mutation
frequencies specified either in the Alexandrov or the Shiraishi format
(must match the format used for \code{signatures}, see above).}

\item{plot}{(Optional) If \code{FALSE} (default), the numerical results
(see below) will be returned. If \code{TRUE}, the reconstructed mutation
frequencies will be plotted against the original, observed mutation
frequencies and the numerical results will be integrated as text labels
in the plot.}
}
\value{
A named list object containing measurements for the Pearson
correlation coefficient between the reconstructed and observed mutation
frequencies, and the explained variance; or alternatively, a plot with
these measurements (see option \code{plot} above).
}
\description{
`evaluateDecompositionQuality()` evaluates the quality of the decomposition
into exposures of a single tumor. The function evaluates the quality of a
decomposition obtained from the function \code{decomposeTumorGenomes}
by comparing the re-composed (=re-constructed) tumor genome mutation
frequencies to those actually observed in the tumor genome. Tumor genome
mutation frequencies are reconstructed using
\code{composeGenomesFromExposures} and the results can optionally be plotted.
}
\examples{

### get Alexandrov signatures from COSMIC
signatures <- readAlexandrovSignatures()

### load preprocessed breast cancer genomes (object 'genomes') from
### Nik-Zainal et al (PMID: 22608084) 
gfile <- system.file("extdata",
         "Nik-Zainal_PMID_22608084-genomes-Alexandrov_3bases.Rdata", 
         package="decompTumor2Sig")
load(gfile)

### compute exposures
exposures <- decomposeTumorGenomes(genomes, signatures, verbose=FALSE)

### evaluate the decomposition by comparing to the original data
evaluateDecompositionQuality(exposures[[1]], signatures, genomes[[1]])

}
\references{
\url{http://rmpiro.net/decompTumor2Sig/}\cr
Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
signatures active in individual tumors. BMC Bioinformatics
20(Suppl 4):152.\cr
}
\seealso{
\code{\link{decompTumor2Sig}}\cr
\code{\link{decomposeTumorGenomes}}\cr
\code{\link{composeGenomesFromExposures}}\cr
\code{\link{computeExplainedVariance}}
}
\author{
Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
M. Piro\cr E-Mail: <rmpiro@gmail.com> or <rosariomichael.piro@polimi.it>
}
