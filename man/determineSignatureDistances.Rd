\name{determineSignatureDistances}
\alias{determineSignatureDistances}
\title{
  determineSignatureDistances
}
\description{
  Determines all similarities between a given target signature
  (of type Alexandrov or Shiraishi) with a set of other signatures (of
  the same type). This can help to compare signatures that have been
  determined in different ways or from different datasets. 
  Different distance measures can be used (see below).
}
\usage{
determineSignatureDistances(target, signatures, method="euclidean")
}
\arguments{
  \item{target}{
    (Mandatory) A single signature of the Alexandrov (vector) or
    Shiraishi type (data frame or matrix). 
  }
  \item{signatures}{
    (Mandatory) The list of signatures for which the distances to the
    \code{target} need to be computed. The signatures must be of the
    same type as \code{target}. 
  }
  \item{method}{
    (Optional) The distance measure to be used. This can be one of the
    following: "\code{frobenius}" for Frobenius distance between
    matrices (only for Shiraishi signatures); "\code{rss}" for the
    residual sum of squares (squared error); or any distance measure
    available for the function \code{dist()} of the \code{stats}
    package. Default: "\code{euclidean}".
  }
}
%\details{
%}
\value{
  A signature-named vector containing all distances. This vector has the
  same order as the signature list, so it is not sorted according to
  distance.
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
  mapSignatureSets\cr
}
\examples{
### get Alexandrov signatures from COSMIC
signAlexandrov <- loadAlexandrovSignatures()

### convert them to Shiraishi signatures
signAlex2Shi <- convertAlexandrov2Shiraishi(signAlexandrov)

### define an arbitrary signature just for testing
### (similar to signature 1)
testSig <- matrix(c(0.1,   0, 0.7, 0.1, 0.1,   0,
                    0.3, 0.2, 0.3, 0.2,   0,   0,
                    0.2, 0.1, 0.5, 0.2,   0,   0), nrow=3, byrow=TRUE) 

### compute distances of the test signature to the converted
### Alexandrov signatures from COSMIC
determineSignatureDistances(testSig, signAlex2Shi, method="frobenius")
}
\keyword{mutational signatures}
\keyword{mutation signatures}
\keyword{mutations}
\keyword{tumor genomes}
