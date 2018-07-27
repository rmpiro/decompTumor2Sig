\name{loadAlexandrovSignatures}
\alias{loadAlexandrovSignatures}
\title{
  loadAlexandrovSignatures
}
\description{
  Loads a set Alexandrov signatures from a flat file or URL. Signatures
  must be specified in the tab-separated format used by the COSMIC
  website; see http://cancer.sanger.ac.uk/cosmic/signatures -> Download
  signatures. 

  Example:

  \tabular{llllll}{
    Subst. \tab Trinucleotide \tab Mutation Type \tab Signature 1 \tab Signature 2 \tab...\cr
    C>A \tab ACA \tab A[C>A]A \tab 0.011098326166 \tab 0.000682708227 \tab...\cr 
    C>A \tab ACC \tab A[C>A]C \tab 0.009149340734 \tab 0.000619107232 \tab...\cr
    C>A \tab ACG \tab A[C>A]G \tab 0.001490070468 \tab 0.000099278956 \tab...\cr
    C>A \tab ACT \tab A[C>A]T \tab 0.006233885236 \tab 0.000323891363 \tab...\cr
    [...]\cr
    T>G \tab TTG \tab T[T>G]G \tab 0.002031076880 \tab 0.000206615168 \tab...\cr
    T>G \tab TTT \tab T[T>G]T \tab 0.004030128160 \tab 0.000023598204 \tab...\cr
  }
}
\usage{
%loadAlexandrovSignatures(file="http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt")
loadAlexandrovSignatures(file)
}
\arguments{
  \item{file}{
    (Mandatory) Can be a single file name or an URL for
    download. Default (COSMIC): 
    "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
  }
}
%\details{
%}
\value{
  A list of Alexandrov signatures that can be used for
  \code{decomposeTumorGenomes}. 
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
  loadShiraishiSignatures\cr
}
\examples{
### get Alexandrov signatures from COSMIC
signatures <- loadAlexandrovSignatures()
}
\keyword{mutational signatures}
\keyword{mutation signatures}
\keyword{mutations}
\keyword{tumor genomes}
