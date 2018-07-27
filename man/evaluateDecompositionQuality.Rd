\name{evaluateDecompositionQuality}
\alias{evaluateDecompositionQuality}
\title{
  evaluateDecompositionQuality
}
\description{
  Evaluates the quality of a decomposition with
  \code{decomposeTumorGenomes} by comparing the re-composed
  (=re-constructed) tumor genome mutation frequencies to those actually
  observed in the tumor genome. Tumor genome mutation frequencies are
  reconstructed using \code{composeGenomesFromExposures} and the results
  can optionally be plotted.
}
\usage{
evaluateDecompositionQuality(exposure, signatures, genome, plot=FALSE)
}
\arguments{
  \item{exposure}{
    (Mandatory) A single vector containing the estimated signature
    contributions/exposures of a single tumor as provided by the function
    \code{decomposeTumorGenomes}. The number of elements of the exposure
    vector must correspond to the number of signatures (see below). 
  }
  \item{signatures}{
    (Mandatory) The list of signatures (vectors, data frames or
    matrices) for which the exposures were obtained. Each of the
    list objects represents one mutational signature. Vectors are used
    for Alexandrov signatures, data frames or matrices for Shiraishi
    signatures.
  }
  \item{genome}{
    (Mandatory) A single tumor genome in form of mutation frequencies
    specified either in the Alexandrov or the Shiraishi format (must
    match the format used for \code{signatures}, see above). 
  }
  \item{plot}{
    (Optional) If \code{FALSE} (default), the numerical results (see
    below) will be returned. If \code{TRUE}, the reconstructed mutation
    frequencies will be plotted against the original, observed mutation
    frequencies and the numerical results will be integrated as text
    labels in the plot.
  }
}
%\details{
%}
\value{
  A named list object containing measurements for the Pearson
  corrleation coefficient between the reconstructed and observerd
  mutation frequencies, and the explained variance; or alternatively, a
  plot with these measurements (see option \code{plot} above).
}
\references{
Krueger, Piro (2017) Identification of Mutational Signatures Active in
Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
Italy. PeerJ Preprints 5:e3257v1, 2017.\cr
http://rmpiro.net/decompTumor2Sig/
}
\author{
  Rosario M. Piro\cr
  Freie Universitaet Berlin\cr
  Maintainer: Rosario M. Piro\cr
  E-Mail: <r.piro@fu-berlin.de> or <rmpiro@gmail.com>
}
%\note{
%}


\seealso{
  decompTumor2Sig\cr
  decomposeTumorGenomes\cr
  composeGenomesFromExposures\cr
  computeExplainedVariance\cr
}
\examples{
### get Alexandrov signatures from COSMIC
signatures <- loadAlexandrovSignatures()

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
\keyword{mutational signatures}
\keyword{mutation signatures}
\keyword{mutations}
\keyword{tumor genomes}
