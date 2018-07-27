\name{plotExplainedVariance}
\alias{plotExplainedVariance}
\title{
  plotExplainedVariance
}
\description{
  For a single genome and a given set of signatures, the function plots
  the variance of the genome's mutation load that can be explained with
  an increasing number of signatures (increasing subset of signatures).
  For each number K of signatures, the highest variance explained by any
  possible subset of K signatures will be plotted. This can help to
  evaluate what minimum threshlod for the explained variance should be
  used to decompose tumor genomes with the function
  \code{decomposeTumorGenomes}.
}
\usage{
plotExplainedVariance(genome, signatures, minExplainedVariance=NULL,
                      minNumSignatures=2, maxNumSignatures=NULL)
}
\arguments{
  \item{genome}{
    (Mandatory) The mutation load of a single genome in Alexandrov- of
    Shiraishi-format, i.e. as vector or matrix. The format must be the
    same as the one used for the \code{signatures} (see below).
  }
  \item{signatures}{
    (Mandatory) The list of signatures (vectors, data frames or
    matrices) which are to be evaluated. Each of the list objects
    represents one mutational signature. Vectors are used for Alexandrov
    signatures, data frames or matrices for Shiraishi signatures.
  }
  \item{minExplainedVariance}{
    (Optional) If a numeric value between 0 and 1 is specified, the plot
    highlights the smallest subset of signatures that is sufficient to
    explain at least the specified fraction of the variance of the
    genome's mutation load. If, for example, \code{minExplainedVariance}
    is 0.99 the smallest subset of signatures that explains at least
    99\% of the variance will be highlighted.
  }
  \item{minNumSignatures}{
    (Optional) The plot will be generated only for
    K>=\code{minNumSignatures}. 
  }
  \item{maxNumSignatures}{
    (Optional) The plot will be generated only for
    K<=\code{minNumSignatures}. 
  }
}
%\details{
%}
\value{
  No return value. The function creates a plot of the explained variance
  as a function of the number of signatures.
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
  computeExplainedVariance
}
\examples{
### get 15 pre-processed Shiraishi signatures computed (object 'signatures') 
### from 435 tumor genomes Alexandrov et al (PMID: 23945592)
### using the pmsignature package
sfile <- system.file("extdata",
         "Alexandrov_PMID_23945592_435_tumors-pmsignature-15sig.Rdata", 
         package="decompTumor2Sig")
load(sfile)

### load preprocessed breast cancer genomes (object 'genomes') from
### Nik-Zainal et al (PMID: 22608084) 
gfile <- system.file("extdata",
         "Nik-Zainal_PMID_22608084-genomes-Shiraishi_5bases_trDir.Rdata", 
         package="decompTumor2Sig")
load(gfile)

### plot the explained variance for 2 to 6 signatures of the first genome
plotExplainedVariance(genomes[[1]], signatures,
         minExplainedVariance=0.98, minNumSignatures=2, maxNumSignatures=6)
}
\keyword{mutational signatures}
\keyword{mutation signatures}
\keyword{mutations}
\keyword{tumor genomes}
