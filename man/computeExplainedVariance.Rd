\name{computeExplainedVariance}
\alias{computeExplainedVariance}
\title{
  computeExplainedVariance
}
\description{
  For a single genome or a set of genomes, the function computes the
  explained variance(s) of the estimated signature contributions/exposures. 
}
\usage{
computeExplainedVariance(exposures, signatures, genomes)
}
\arguments{
  \item{exposures}{
    (Mandatory) A single vector or list of vectors containing the
    estimated signature contributions/exposures as provided by the
    function \code{decomposeTumorGenomes}. A list of vectors is used if
    the explained variance shall be computed for multiple genomes. The
    number of exposure vectors must correspond to the number of genomes
    (see below). The number of elements of each exposure vector must
    correspond to the number of signatures (see below).
  }
  \item{signatures}{
    (Mandatory) The list of signatures (vectors, data frames or
    matrices) for which the exposures were obtained. Each of the
    list objects represents one mutational signature. Vectors are used
    for Alexandrov signatures, data frames or matrices for Shiraishi
    signatures.
  }
  \item{genomes}{
    (Mandatory) Can be either a vector, a data frame or a matrix (for an
    individual tumor genome), or a list of one of these object types (for
    multiple tumors). Each tumor genome must be of the same form as
    the '\code{signatures}' (see above).
  }
}
%\details{
%}
\value{
  A numeric vector of explained variances, one for each genome.
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
  plotExplainedVariance
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

### compute explained variance for the tumor genomes
computeExplainedVariance(exposures, signatures, genomes)
}
\keyword{mutational signatures}
\keyword{mutation signatures}
\keyword{mutations}
\keyword{tumor genomes}
