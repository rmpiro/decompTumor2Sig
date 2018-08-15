#' decompTumor2Sig - Decomposition of individual tumors into mutational
#' signatures by signature refitting
#' 
#' The decompTumor2Sig package uses quadratic programming to decompose the
#' somatic mutation catalog from an individual tumor sample (or multiple
#' individual tumor samples) into a set of given mutational signatures (either
#' of the "Alexandrov model" by Alexandrov et al, Nature 500(7463):415-421,
#' 2013, or the "Shiraishi model" by Shiraishi et al, PLoS Genet
#' 11(12):e1005657, 2015), thus computing weights (or "exposures") that reflect
#' the contributions of the signatures to the mutation load of the tumor.
#' 
#' \tabular{ll}{
#' Package: \tab decompTumor2Sig\cr
#' Type: \tab Package\cr
#' Version: \tab 1.3.1\cr
#' Date: \tab 2018-08-15\cr
#' License: \tab GPL (>=2)\cr
#' }
#' 
#' The package provides the following functions:
#' 
#' \tabular{ll}{
#' 
#' composeGenomesFromExposures():\tab (re-)construct tumor genome mutation\cr
#'                         \tab frequencies from the signatures and\cr
#'                         \tab their corresponding exposures, or\cr
#'                         \tab contributions.\cr
#'
#' computeExplainedVariance():\tab determine the variance explained by \cr
#'                         \tab estimated signature contributions\cr
#'                         \tab (i.e., exposures to signatures).\cr
#'
#' convertAlexandrov2Shiraishi():\tab convert a set of Alexandrov \cr
#'                         \tab signatures to Shiraishi signatures.\cr
#'
#' convertGenomesFromVRanges():\tab convert a genome or set of genomes\cr
#'                         \tab from a \code{VariantAnnotation::VRanges}\cr
#'                         \tab object.\cr
#'
#' decomposeTumorGenomes():\tab determine the weights/contributions of\cr
#'                         \tab a set of signatures to each of a set of\cr
#'                         \tab individual tumor genomes.\cr
#' 
#' determineSignatureDistances():\tab for a given signature \cr
#'                         \tab compute its distances to each of a set\cr
#'                         \tab of target signatures.\cr
#'
#' downgradeShiraishiSignatures():\tab downgrade Shiraishi signatures\cr
#'                         \tab by removing flanking bases and/or the\cr
#'                         \tab transcription direction.\cr
#' 
#' evaluateDecompositionQuality():\tab evaluate the quality of a\cr
#'                         \tab decomposition by comparing the\cr
#'                         \tab re-composed (=re-constructed) tumor\cr
#'                         \tab mutation frequencies to those actually\cr
#'                         \tab observed in the tumor genome.\cr
#'
#' getGenomesFromMutFeatData():\tab extract the genomes from a \cr
#'                         \tab \code{MutationFeatureData} object as \cr
#'                         \tab provided by, for example,\cr
#'                         \tab \code{pmsignature::readMPFile}.\cr
#'
#' getSignaturesFromEstParam():\tab extract a set of signatures from an \cr
#'                         \tab \code{EstimatedParameters} object as\cr
#'                         \tab returned by function \code{getPMSignature}\cr
#'                         \tab of the \code{pmsignature} package.\cr
#'
#' mapSignatureSets():     \tab find a mapping from one signature\cr
#'                         \tab set to another.\cr
#'
#' plotDecomposedContribution():\tab plot the decomposition of a\cr
#'                         \tab genome into mutational signatures\cr
#'                         \tab (i.e., the contributions of, or\cr
#'                         \tab exposures to, the signatures).\cr
#'
#' plotExplainedVariance():\tab plot the variance of a genome's\cr
#'                         \tab mutation patterns which can be\cr
#'                         \tab explained with an increasing number\cr
#'                         \tab of signatures.\cr
#'
#' plotMutationDistribution():\tab plot a single signature or the\cr
#'                         \tab mutation frequency data for a single\cr
#'                         \tab genome.\cr
#'
#' readAlexandrovSignatures():\tab read Alexandrov signatures in the\cr
#'                         \tab COSMIC format from a flat file or URL.\cr
#'
#' readGenomesFromMPF():   \tab read a genome or set of genomes from a \cr
#'                         \tab Mutation Position Format (MPF) file.\cr
#'
#' readGenomesFromVCF():   \tab read a genome or set of genomes from a\cr
#'                         \tab Variant Call Format (VCF) file.\cr
#'
#' readShiraishiSignatures():\tab read Shiraishi signatures from\cr
#'                         \tab flat files.\cr
#'
#' }
#' 
#' @name decompTumor2Sig-package
#' @aliases decompTumor2Sig-package decompTumor2Sig
#' @docType package
#' @author Rosario M. Piro [aut, cre], Sandra Krueger [ctb]\cr
#' Freie Universitaet Berlin\cr
#' Maintainer: Rosario M. Piro\cr
#' E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
NULL