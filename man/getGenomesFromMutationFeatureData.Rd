\name{getGenomesFromMutationFeatureData}
\alias{getGenomesFromMutationFeatureData}
\title{
  getGenomesFromMutationFeatureData
}
\description{
  Takes a \code{MutationFeatureData} object (mutation count data) as
  read by the '\code{pmsignature}' package (e.g., by
  \code{pmsignature::readMPFile}) and extracts the mutation counts. For
  passing them to \code{decomposeTumorGenomes}, the mutation counts must
  be normalized to mutation fractions, which is done by
  default. IMPORTANT: set \code{normalize} to \code{FALSE} only if you
  are interested in full integer counts, but do not pass unnormalized
  counts to \code{decomposeTumorGenomes}! 
}
\usage{
getGenomesFromMutationFeatureData(countData, normalize=TRUE)
}
\arguments{
  \item{countData}{
    (Mandatory) A \code{MutationFeatureData} object as constructed, for
    example, by \code{pmsignature::readMPFile}.
  }
  \item{normalize}{
    (Optional) Boolean value to specify whether to normalize the
    mutation count data to mutation fractions between 0 and 1. This is
    the default and NECESSARY in case you want to pass the return value
    to \code{decomposeTumorGenomes}. Set \code{normalize} to
    \code{FALSE} only if you are interested in full integer counts, but
    do not pass unnormalized counts to \code{decomposeTumorGenomes}!
  }
}
%\details{
%}
\value{
  A list of (normalized) mutation counts, one object per genome. The
  format is the same table used by the corresponding Shiraishi
  signatures. 
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
decompTumor2Sig
}
\examples{
### get breast cancer genomes from 
### Nik-Zainal et al (PMID: 22608084) in the format produced by
### pmsignature (PMID: 26630308)
pmsigdata <- system.file("extdata", 
         "Nik-Zainal_PMID_22608084-pmsignature-G.Rdata", 
         package="decompTumor2Sig")
load(pmsigdata)

### extract the genomes from the pmsignature G object
genomes <- getGenomesFromMutationFeatureData(G, normalize=TRUE)
}
\keyword{mutational signatures}
\keyword{mutation signatures}
\keyword{mutations}
\keyword{tumor genomes}
