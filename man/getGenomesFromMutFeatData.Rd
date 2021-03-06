% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGenomesFromMutFeatData.R
\name{getGenomesFromMutFeatData}
\alias{getGenomesFromMutFeatData}
\title{Get genomes (mutation frequencies) from \code{MutationFeatureData}.}
\usage{
getGenomesFromMutFeatData(mutFeatData, normalize=TRUE)
}
\arguments{
\item{mutFeatData}{(Mandatory) A \code{MutationFeatureData} object as
constructed, for example, by \code{pmsignature::readMPFile}.}

\item{normalize}{(Optional) Boolean value to specify whether to normalize
the mutation count data to mutation fractions between 0 and 1. This is
the default and NECESSARY in case you want to pass the return value to
\code{decomposeTumorGenomes}. Set \code{normalize} to \code{FALSE} only
if you are interested in full integer counts, but do not pass unnormalized
counts to \code{decomposeTumorGenomes}!}
}
\value{
A list of mutation frequencies (or mutation counts if not
normalized), one object per genome. The format is either according to the
Shiraishi or the Alexandrov model, depending on how the mutation data was
loaded with \code{pmsignature}.
}
\description{
`getGenomesFromMutFeatData()` takes a \code{MutationFeatureData} object
(mutation count data) as read by the '\code{pmsignature}' package (e.g.,
by \code{pmsignature::readMPFile}, version 0.3.0) and extracts the mutation
counts of the genomes therein. For passing the genomes to
\code{decomposeTumorGenomes}, the mutation counts must be normalized to
mutation frequencies, which is done by default. [IMPORTANT: set
\code{normalize} to \code{FALSE} only if you are interested in full
integer counts, but do not pass unnormalized counts to
\code{decomposeTumorGenomes}!]
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
genomes <- getGenomesFromMutFeatData(G, normalize=TRUE)

}
\references{
\url{http://rmpiro.net/decompTumor2Sig/}\cr
Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
signatures active in individual tumors. BMC Bioinformatics
20(Suppl 4):152.\cr
}
\seealso{
\code{\link{decompTumor2Sig}}
}
\author{
Rosario M. Piro, Politecnico di Milano\cr
Sandra Krueger, Freie Universitaet Berlin\cr Maintainer: Rosario
M. Piro\cr E-Mail: <rmpiro@gmail.com> or <rosariomichael.piro@polimi.it>
}
