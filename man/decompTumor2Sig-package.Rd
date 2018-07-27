\name{decompTumor2Sig-package}
\alias{decompTumor2Sig-package}
\alias{decompTumor2Sig}
\docType{package}
\title{
decompTumor2Sig
}
\description{
  The decompTumor2Sig package uses quadratic programming to decompose
  the somatic mutation catalog from an individual tumor sample (or
  multiple individual tumor samples) into a set of given mutational
  signatures (either of the "Alexandrov model" by Alexandrov et al,
  Nature 500(7463):415-421, 2013), or the "Shiraishi model" by Shiraishi
  et al, PLoS Genet 11(12):e1005657, 2015), thus computing weights (or
  "exposures") that reflect the contributions of the signatures to the
  mutation load of the tumor.

  The package additionally provides helper functions to extract genomes
  (mutation catalogs) and signatures from objects provided by
  the \code{pmsignature} package (Shiraishi et al, 2015) or to read them
  from files. 
}
\details{
\tabular{ll}{
Package: \tab decompTumor2Sig\cr
Type: \tab Package\cr
Version: \tab 1.3.0\cr
Date: \tab 2018-07-26\cr
License: \tab GPL (>=2)\cr
}
The package provides the following functions:
\tabular{ll}{
  decomposeTumorGenomes():\tab determines the weights/contributions of a\cr
  \tab set of SIGNATURES to each of a set of \cr
  \tab individual tumor GENOMES.\cr 
  loadGenomesFromVCF():\tab loads a genome or set of genomes from a\cr
  \tab Variant Call Format (VCF) file.\cr 
  loadGenomesFromMPF():\tab loads a genome or set of genomes from a \cr
  \tab Mutation Position Format (MPF) file.\cr
  convertGenomesFromVRanges():\tab convert a genome or set of genomes \cr
  from a \code{VariantAnnotation::VRanges} object.\cr
  getGenomesFromMutationFeatureData():\tab extracts the GENOMES from a \cr
  \tab \code{MutationFeatureData} object as provided \cr
  \tab  by, for example, \code{pmsignature::readMPFile}.\cr
  loadShiraishiSignatures():\tab loads Shiraishi signatures from flat files.\cr 
  getSignatureListFromEstimatedParameters():\tab extracts a set of 
  SIGNATURES from an \cr
  \tab \code{EstimatedParameters} object as computed by \cr
  \tab \code{pmsignature::getPMSignature}.\cr
  loadAlexandrovSignatures():\tab loads Alexandrov signatures in the \cr
  \tab COSMIC format from a flat file or URL.\cr
  convertAlexandrov2Shiraishi():\tab converts a set of Alexandrov \cr
  \tab signatures to Shiraishi signatures.\cr
  computeExplainedVariance():\tab determines the variance explained by \cr
  \tab estimated signature contributions/exposures.\cr
  plotExplainedVariance():\tab plots the variance of the genome's \cr
  \tab mutation load that can be explained with an increasing \cr
  \tab number of signatures.\cr
  plotDecomposedContribution():\tab plot a the decomposition/exposures \cr
  \tab of a genome to the mutational signatures.\cr
  plotMutationDistribution():\tab plot a single signature or mutation \cr
  \tab frequency data for a single genome.\cr
  composeGenomesFromExposures():\tab (re-)construct tumor genome mutation \cr 
  \tab frequencies from the signatures and their\cr
  \tab corrsponding exposures/contributions.\cr
  evaluateDecompositionQuality():\tab evaluates the quality of a \cr
  \tab decomposition by comparing the re-composed\cr
  \tab (=re-constructed) tumor genome mutation\cr
  \tab frequencies to those actually observed in\cr
  \tab the tumor genome.\cr
  determineSignatureDistances():\tab for a given target signature \cr
  \tab compute its distances to each of a set\cr
  \tab of other signatures.\cr
  mapSignatureSets():\tab find a mapping from one signature set to\cr
  \tab another.\cr
  downgradeShiraishiSignatures():\tab downgrades Shiraishi signatures\cr
  \tab by removing flanking bases and/or the transcription direction.\cr
  
}
}
\author{
  Rosario M. Piro, Sandra Krueger\cr
  Freie Universitaet Berlin\cr
  Maintainer: Rosario M. Piro\cr
  E-Mail: <r.piro@fu-berlin.de> or <rmpiro@gmail.com>
}
\references{
  Krueger, Piro (2017) Identification of Mutational Signatures Active in
  Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
  Personalized Medicine in the Big Data Era, October 16-18, 2017,
  Palermo, Italy. PeerJ Preprints 5:e3257v1, 2017.\cr
  http://rmpiro.net/decompTumor2Sig/
}
\keyword{mutational signatures}
\keyword{mutation signatures}
\keyword{mutations}
\keyword{tumor genomes}
