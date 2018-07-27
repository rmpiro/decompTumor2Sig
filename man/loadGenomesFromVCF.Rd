\name{loadGenomesFromVCF}
\alias{loadGenomesFromVCF}
\title{
  loadGenomesFromVCF
}
\description{
  Loads a single tumor genome (sample) or a set of genomes from a VCF
  file and determines the mutation frequencies according to a specific
  model of mutational signatures (Alexandrov or Shiraishi).
}
\usage{
loadGenomesFromVCF(file, numBases=5, type="Shiraishi", trDir=TRUE,
refGenome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
transcriptAnno=TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
verbose=TRUE) 
}
\arguments{
  \item{file}{
    (Mandatory) The name of the VCF file (can be compressed with
    \code{gzip}). 
  }
  \item{numBases}{
    (Mandatory) Total number of bases (mutated base and flanking bases)
    to be used for sequence patterns. Default: 5
  }
  \item{type}{
    (Mandatory) Signature model or type ('\code{Alexandrov}' or
    '\code{Shiraishi}'). Default: \code{Shiraishi}
  }
  \item{trDir}{
    (Mandatory) Specifies whether the transcription direction is taken
    into account in the signature model. If so, only mutations within
    genes can be considered. Default: \code{TRUE}
  }
  \item{refGenome}{
    (Mandatory) The reference genome (\code{BSgenome}) needed to extract
    sequence patterns. Default: \code{BSgenome} object for hg19.
  }
  \item{transcriptAnno}{
    (Optional) Transcript annotation (\code{TxDb} object) used to
    determine the transcription direction. This is required only if
    \code{trDir} is \code{TRUE}. Default: \code{TxDb} object for hg19.
  }
  \item{verbose}{
    (Optional) Print information about reading and processing the
    mutation data. Default: \code{TRUE}
  }
}
%\details{
%}
\value{
   A list containing the genomes in terms of frequencies of the mutated
   sequence patterns. This list of genomes can be used for
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
  loadGenomesFromVCF\cr
  getGenomesFromMutationFeatureData\cr
}
\examples{
### load reference genome and transcript annotation (if direction is needed)
refGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
transcriptAnno <-
  TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene

### load breast cancer genomes from Nik-Zainal et al (PMID: 22608084) 
gfile <- system.file("extdata",
         "Nik-Zainal_PMID_22608084-VCF-convertedfromMPF.vcf.gz", 
         package="decompTumor2Sig")
genomes <- loadGenomesFromVCF(gfile, numBases=5, type="Shiraishi",
         trDir=TRUE, refGenome=refGenome, transcriptAnno=transcriptAnno,
         verbose=FALSE)
}
\keyword{mutational signatures}
\keyword{mutation signatures}
\keyword{mutations}
\keyword{tumor genomes}
