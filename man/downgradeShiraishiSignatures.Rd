\name{downgradeShiraishiSignatures}
\alias{downgradeShiraishiSignatures}
\title{
  downgradeShiraishiSignatures
}
\description{
  Downgrades/trims Shiraishi signatures by dropping flanking bases
  (reducing the length of the sequence pattern and/or the transcription
  direction. This can be easily down because the flanking bases and the
  transcription direction are considered as independent features
  according to the Shiraishi model of mutational signatures.
}
\usage{
downgradeShiraishiSignatures(signatures, numBases=NULL, removeTrDir=FALSE)
}
\arguments{
  \item{signatures}{
    (Mandatory) A list of Shiraishi signatures that need to be
    downgraded/trimmed.
  }
  \item{numBases}{
    (Conditionally optional) The total number of bases (mutated base
    plus flanking bases around the mutated base) that should be
    kept. All further flanking bases farther away from the mutated bases
    are dropped. If specified, \code{numBases} must be odd an smaller
    than the \code{signatures}' current number of bases. If \code{NULL},
    no flanking bases will be dropped. At least one of \code{numBases}
    or \code{removeTrDir} must be specified. 
  }
  \item{removeTrDir}{
    (Conditionally optional) Logical value that specifies whether
    information on the transcript direction should be dropped (if
    present at all). At least one of \code{numBases} or
    \code{removeTrDir} must be specified.
  }
}
%\details{
%}
\value{
  A list of Shiraishi signatures that have been accordingly downgraded.
}
\references{
Krueger, Piro (2017) Identification of Mutational Signatures Active in
Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
Italy. PeerJ Preprints 5:e3257v1, 2017.\cr
http://rmpiro.net/decompTumor2Sig/
}
\author{
  Rosario M. Piro, Sandra Krueger\cr
  Freie Universitaet Berlin\cr
  Maintainer: Rosario M. Piro\cr
  E-Mail: <r.piro@fu-berlin.de> or <rmpiro@gmail.com>
}
%\note{
%}

\seealso{
  decompTumor2Sig\cr
}
\examples{
### Load 15 Shiraishi signatures obtained from 435 tumor genomes from
### Alexandrov et al. (number of bases: 5, transcription direction: yes)
sfile <- system.file("extdata",
         "Alexandrov_PMID_23945592_435_tumors-pmsignature-15sig.Rdata", 
         package="decompTumor2Sig")
load(sfile)

### downgrade the signatures to include only 3 bases and drop the
### trancription direction 
downgradeShiraishiSignatures(signatures, numBases=3, removeTrDir=TRUE)
}
\keyword{mutational signatures}
\keyword{mutation signatures}
\keyword{mutations}
\keyword{tumor genomes}
