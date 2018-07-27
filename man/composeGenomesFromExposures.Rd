\name{composeGenomesFromExposures}
\alias{composeGenomesFromExposures}
\title{
  composeGenomesFromExposures
}
\description{
  (Re-)compose/construct tumor genome characteristics, or mutation
  frequencies, from the signatures and their corresponding exposures, or
  contributions. The (re-)composition is performed by computing the
  weighted sum of the mutational signatures, where the weights
  are to the exposures (=contributions) of the corresponding
  signatures. This can, for example, be used to verify that a
  decomposition obtained from \code{decomposeTumorGenomes} is
  meaningful. 
}
\usage{
composeGenomesFromExposures(exposures, signatures)
}
\arguments{
  \item{exposures}{
    (Mandatory) A single vector or list of vectors containing the
    estimated signature contributions/exposures as provided by the
    function \code{decomposeTumorGenomes}. A list of vectors is used if
    the (re-)composition shall be performed for multiple genomes. The
    number of elements of each exposure vector must correspond to the
    number of signatures (see below). 
  }
  \item{signatures}{
    (Mandatory) The list of signatures (vectors, data frames or
    matrices) for which the exposures were obtained. Each of the
    list objects represents one mutational signature. Vectors are used
    for Alexandrov signatures, data frames or matrices for Shiraishi
    signatures.
  }
}
%\details{
%}
\value{
  A list of "predicted" genomes, i.e., their mutational patterns
  computed as weighted sums of the mutational signatures, where the
  weights correspond to the exposures/contributions of the
  corresponding signatures.
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

### re-compose (predict) tumor genome features from exposures
predGenomes <- composeGenomesFromExposures(exposures, signatures)
}
\keyword{mutational signatures}
\keyword{mutation signatures}
\keyword{mutations}
\keyword{tumor genomes}
