\name{plotMutationDistribution}
\alias{plotMutationDistribution}
\title{
  plotMutationDistribution
}
\description{
  Plots a single signature or mutation frequency data for a single
  genome. This works for signatures or genome data of both the Shiraishi
  and the Alexandrov type.

  [IMPORTANT: The function requires the '\code{pmsignature}' package to
  be installed (Shiraishi et al. PLoS Genet 11(12):e1005657, 2015)!]
}
\usage{
plotMutationDistribution(mutData)
}
\arguments{
  \item{mutData}{
    (Mandatory) The signature or genome mutation frequency data to be 
    plotted. This can either be a matrix (Shiraishi-type model) or a
    numeric vector (Alexandrov-type model).
  }
}
%\details{
%}
\value{
  No return value. The creates a plot on the standard graphical output.
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
}
\examples{
### Attention: using plotMutationDistribution requires the package
### pmsignature to be installed!


### get Alexandrov signatures from COSMIC
signatures <- loadAlexandrovSignatures()

### plot the first Alexandrov signature
\donttest{
plotMutationDistribution(signatures[[1]])
}

### load four Shiraishi signatures for breast cancer genomes from 
### Nik-Zainal et al (PMID: 22608084) from flat files
sigfiles <- system.file("extdata",
         paste0("Nik-Zainal_PMID_22608084-pmsignature-sig",1:4,".tsv"),
         package="decompTumor2Sig")
signatures <- loadShiraishiSignatures(sigfiles)

### plot the first Shiraishi signature
\donttest{
plotMutationDistribution(signatures[[1]])
}
}
\keyword{mutational signatures}
\keyword{mutation signatures}
\keyword{mutations}
\keyword{tumor genomes}
