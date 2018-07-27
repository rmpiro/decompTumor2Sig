\name{convertAlexandrov2Shiraishi}
\alias{convertAlexandrov2Shiraishi}
\title{
  convertAlexandrov2Shiraishi
}
\description{
  Converts a set Alexandrov signatures (loaded with
  \code{loadAlexandrovSignatures}) to the Shiraishi format, summing the 
  respective frequencies of base changes, and upstream and downstream
  flanking bases. The resulting Shiraishi signatures don't provide
  information on the transcription strand, as this is not part of the
  Alexandrov signatures.

  [Attention: this conversion is experimental and the applicability of
  Shiraishi signatures derived from Alexandrov signatures has not been
  explored!]
}
\usage{
convertAlexandrov2Shiraishi(signatures)
}
\arguments{
  \item{signatures}{
    (Mandatory) A list of Alexandrov signatures with named elements as
  produced by \code{loadAlexandrovSignatures}.
  }
}
%\details{
%}
\value{
  A list of Shiraishi signatures that can be used for decomposeTumorGenomes.
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
  loadAlexandrovSignatures\cr
  loadShiraishiSignatures\cr
}
\examples{
### get Alexandrov signatures from COSMIC
signAlexandrov <- loadAlexandrovSignatures()

### convert them to the Shiraishi model
signShiraishi <- convertAlexandrov2Shiraishi(signAlexandrov)
}
\keyword{mutational signatures}
\keyword{mutation signatures}
\keyword{mutations}
\keyword{tumor genomes}
