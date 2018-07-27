\name{plotDecomposedContribution}
\alias{plotDecomposedContribution}
\title{
  plotDecomposedContribution
}
\description{
  Plots a the decomposition/exposures of a genome to the mutational
  signatures (mutational processes), that is, the contributions of the
  signatures to the mutations observed in a tumor genome.

  These decompositions can be obtained running
  \code{decomposeTumorGenomes()}. 
}
\usage{
plotDecomposedContribution(decomposition, signatures=NULL, removeNA=TRUE)
}
\arguments{
  \item{decomposition}{
    (Mandatory) A decomposition vector (exposure vector) obtained for a
    single tumor genome.
  }
  \item{signatures}{
    (Optional) A list object containing the signatures used to compute
    the decomposition. If specified, the signature labels used in the
    plot will be taken from the element names of the list; otherwise
    signatures will be named from sign_1 to sign_N.
  }
  \item{removeNA}{
    (Optional) If \code{TRUE} (default), signatures with an NA as
    exposure will not be included on the x-axis of the the plot.
    Exposures can NA if they have been determined with a greedy search.
  }
}
%\details{
%}
\value{
  No return value. The function creates a plot of the decomposed tumor
  genome (i.e., contributions of the single signatures). 
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

### plot signature composition of the first genome
plotDecomposedContribution(exposures[[1]], signatures=NULL)
}
\keyword{mutational signatures}
\keyword{mutation signatures}
\keyword{mutations}
\keyword{tumor genomes}
