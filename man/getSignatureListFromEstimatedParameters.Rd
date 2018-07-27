\name{getSignatureListFromEstimatedParameters}
\alias{getSignatureListFromEstimatedParameters}
\title{
  getSignatureListFromEstimatedParameters
}
\description{
  Takes an \code{EstimatedParameters} object (signatures data) as
  computed by the '\code{pmsignature}' package (by
  \code{pmsignature::getPMSignature}) and extracts the signature
  information. This can then be passed to \code{decomposeTumorGenomes}. 
}
\usage{
getSignatureListFromEstimatedParameters(Param)
}
\arguments{
  \item{Param}{
    (Mandatory) An \code{EstimatedParameters} object as the one produced
    by the \code{pmsignature} package's signature contruction method
    \code{pmsignature::getPMSignature}.
  }
}
%\details{
%}
\value{
  A list of Shiraishi signatures, one object per signature. The format
  is a table with mutliple rows, the first for the base substitution,
  then 2*N rows for N flanking bases in each direction, and finally an
  (optional) row for the transcription strand, if it has been taken into
  account. Please see the \code{pmsignature} package or the
  \code{decompTumor2Sig} vignette for more information.
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
### load signatures for breast cancer genomes from 
### Nik-Zainal et al (PMID: 22608084) in the format produced by
### pmsignature (PMID: 26630308)
pmsigdata <- system.file("extdata",
         "Nik-Zainal_PMID_22608084-pmsignature-Param.Rdata", 
         package="decompTumor2Sig")
load(pmsigdata)

### extract the signatures from the pmsignature Param object
\donttest{
signatures <- getSignatureListFromEstimatedParameters(Param)
}
}
\keyword{mutational signatures}
\keyword{mutation signatures}
\keyword{mutations}
\keyword{tumor genomes}
