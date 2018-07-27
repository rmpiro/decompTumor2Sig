\name{decomposeTumorGenomes}
\alias{decomposeTumorGenomes}
\title{
  decomposeTumorGenomes
}
\description{
  Takes a set of mutational signatures and mutation features from one or
  more tumor genomes and computes weights/contributions for each of the
  signatures in each individual genome. Alternatively, the function can
  determine for each genome only a subset of signatures and their
  contributions which are sufficient to exceed a user-given minimum
  threshold for the explained variance of the genome's mutation load.
}
\usage{
decomposeTumorGenomes(genomes, signatures, minExplainedVariance=NULL,
                      minNumSignatures=2, maxNumSignatures=NULL,
                      greedySearch=FALSE,
                      constrainToMaxContribution=FALSE, tolerance=0.1,
                      verbose=FALSE)
}
\arguments{
  \item{genomes}{
    (Mandatory) Can be either a vector, a data frame or a matrix (for an
    individual tumor genome), or a list of one of these object types (for
    multiple tumors). Each tumor genome must be of the same form as
    the '\code{signatures}' (see below).
  }
  \item{signatures}{
    (Mandatory) A list of vectors, data frames or matrices. Each of the
    objects represents one mutational signature. Vectors are used for
    Alexandrov signatures, data frames or matrices for Shiraishi
    signatures. 
  }
  \item{minExplainedVariance}{
    (Optional) If \code{NULL} (default), exactly \code{maxNumSignatures}
    (see below; default: all) will be taken for decomposing each
    genome. If a numeric value between 0 and 1 is specified for
    \code{minExplainedVariance}, for each genome the function will
    select the smallest number of signatures which is sufficient to
    explain at least the specified fraction of the variance of the
    genome's mutation load. E.g., if \code{minExplainedVariance}=0.99
    the smallest subset of signatures that explains at least 99\% of the
    variance is taken. Please note: depending on the number of
    signatures, this may take quite a while because for each number K of
    signatures, all possible subsets composed of K signatures will be
    tested to identify the subset that explains the highest part of the
    variance. If not enough variance is explained, K will be incremented
    by one. Notes: 1) to speed up the search, the parameters
    \code{minNumSignatures}, \code{maxNumSignatures} and
    \code{greedySearch} can be used (see below); 2) for genomes for
    which none of the possible subsets of signatures explains enough
    variance, the returned exposure vector will be set to \code{NULL}. 
  }
  \item{minNumSignatures}{
    (Optional) Used if \code{minExplainedVariance} is specified (see
    above). To find the smallest subset of signatures which explain the
    variance, at least \code{minNumSignatures} will be taken. This can
    be used to reduce the search space in a time-consuming search over a
    large number of signatures.
  }
  \item{maxNumSignatures}{
    (Optional) Used if \code{minExplainedVariance} is specified (see
    above). To find the smallest subset of signatures which explain the
    variance, at most \code{maxNumSignatures} will be taken. This can be
    used to reduce the search space in a time-consuming search over a
    large number of signatures. If \code{maxNumSignatures} is
    \code{NULL} (default), all signatures will be taken as the maximum.
  }
  \item{greedySearch}{
    (Optional) Used if \code{minExplainedVariance} is specified (see
    above). If \code{greedySearch} is \code{TRUE} then not all possible
    combinations of \code{minNumSignatures} to \code{maxNumSignatures}
    signatures will be checked. Instead, first all possible combinations
    for exactly \code{minNumSignatures} will be checked to select the
    best starting set, then iteratively the next best signature will be
    added (maximum increase in explained variability) until
    \code{minExplainedVariance} of the variance can be explained (or 
    \code{maxNumSignatures} is exceeded). NOTE: this is highly
    recommended for large sets of signatures (>15)! 
  }
  \item{constrainToMaxContribution}{
    (Optional) [Note: this is experimental and is usually not needed!] 
    If \code{TRUE}, the maximum contribution that can be attributed
    to a signature will be constraint by the variant feature counts
    (e.g., specific flanking bases) observed in the individual tumor
    genome. If, for example, 30\% of all observed variants have a
    specific feature and 60\% of the variants produced by a mutational
    process/signature will manifest the feature, then the signature can
    have contributed up to 0.3/0.6 (=0.5 or 50\%) of the observed
    variants. The lowest possible contribution over all signature
    features will be taken as the allowed maximum contribution of the
    signature. This allowed maximum will additionally be increased by
    the value specified as \code{tolerance} (see below). For the
    illustrated example and \code{tolerance}=0.1 a contribution of up to
    0.5+0.1 = 0.6 (or 60\%) of the signature would be allowed. 
  }
  \item{tolerance}{
    (Optional) If \code{constrainToMaxContribution} is \code{TRUE}, the
    maximum contribution computed for a signature is increased by this
    value (see above). If the parameter
    \code{constrainToMaxContribution} is \code{FALSE}, the tolerance
    value is ignored. Default: 0.1. 
  }
  \item{verbose}{
    (Optional) If \code{TRUE} some information about the processed
    genome and used number of signatures will be printed.
  }
}
%\details{
%}
\value{
  A list of signature weight vectors (also called 'exposures'), one for
  each tumor genome. E.g., the first vector element of the first list
  object is the weight/contribution of the first signature in the first
  tumor genome. IMPORTANT: If \code{minExplainedVariance} is specified,
  then the exposures of a genome will NOT be returned if the minimum
  explained variance is not reached within the requested minimum and
  maximum numbers of signatures (\code{minNumSignatures} and
  \code{maxNumSignatures})! The corresponding exposure vector will be
  set to \code{NULL}. 
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
### get Alexandrov signatures from COSMIC
signatures <- loadAlexandrovSignatures()

### load reference genome
refGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

### load breast cancer genomes from Nik-Zainal et al (PMID: 22608084) 
gfile <- system.file("extdata",
         "Nik-Zainal_PMID_22608084-VCF-convertedfromMPF.vcf.gz", 
         package="decompTumor2Sig")
genomes <- loadGenomesFromVCF(gfile, numBases=3, type="Alexandrov",
         trDir=FALSE, refGenome=refGenome, verbose=FALSE)

### compute exposures
exposures <- decomposeTumorGenomes(genomes, signatures, verbose=FALSE)

### (for further examples on searching subsets, please see the vignette)
}
\keyword{mutational signatures}
\keyword{mutation signatures}
\keyword{mutations}
\keyword{tumor genomes}
