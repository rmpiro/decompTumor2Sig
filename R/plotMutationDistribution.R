#' Plot mutation frequency data of a mutational signature or tumor genome.
#'
#' `plotMutationDistribution()` plots a single signature or the mutation
#' frequency data for a single genome. This works for signatures or genome
#' data of both the Shiraishi and the Alexandrov type. [IMPORTANT: The
#' function requires the '\code{pmsignature}' package to be installed
#' (Shiraishi et al. PLoS Genet 11(12):e1005657, 2015).]
#'
#' @usage plotMutationDistribution(mutData)
#' @param mutData (Mandatory) The signature or genome mutation frequency data
#' to be plotted. This can either be a matrix (Shiraishi model) or a
#' numeric vector (Alexandrov model).
#' @return Returns (or draws) a plot according to the Alexandrov or Shiraishi
#' model of mutational signatures.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @seealso \code{\link{decompTumor2Sig}}
#' @examples
#' 
#' ### Attention: using plotMutationDistribution requires the package
#' ### pmsignature to be installed!
#' 
#' 
#' ### get Alexandrov signatures from COSMIC
#' signatures <- readAlexandrovSignatures()
#' 
#' ### plot the first Alexandrov signature
#' \donttest{
#' plotMutationDistribution(signatures[[1]])
#' }
#' 
#' ### read four Shiraishi signatures for breast cancer genomes from 
#' ### Nik-Zainal et al (PMID: 22608084) from flat files
#' sigfiles <- system.file("extdata",
#'          paste0("Nik-Zainal_PMID_22608084-pmsignature-sig",1:4,".tsv"),
#'          package="decompTumor2Sig")
#' signatures <- readShiraishiSignatures(sigfiles)
#' 
#' ### plot the first Shiraishi signature
#' \donttest{
#' plotMutationDistribution(signatures[[1]])
#' }
#' 
#' @importFrom methods new
#' @export plotMutationDistribution
plotMutationDistribution <- function(mutData) {

    if (!requireNamespace("pmsignature", quietly=TRUE)) {
        stop(paste("Function plotMutationDistribution requires the",
                   "package pmsignature to be installed."))
    }

    if (isSignatureSet(mutData)) {
        if (length(mutData) > 1) {
            warning(paste("object 'mutData' passed to plotMutationDistribution",
                          "is a list with multiple elements. Taking only the",
                          "first element and trying to plot it"))
        }
        
        mutData <- mutData[[1]]  # take first object from list
    }

    if (!is.probability.object(mutData)) {
        stop("'mutData' must be a single genome or signature!")
    }
    
    # build an "EstimatedParameters" object for visualization with pmsignature
    pmParam <- new(Class = "EstimatedParameters")

    # try to determine the number of bases and whether the transcription
    # direction was considered
    dataFeatures <- determineTypeNumBasesAndTrDir(mutData)
    if (is.null(dataFeatures$type)) {
        stop(paste("Object 'mutData' must be a single Shiraishi- or",
                   "Alexandrov signature!"))
    }

    haveType <- dataFeatures$type

    pmParam <- setNumFlankingBases(pmParam, as.integer(dataFeatures$numBases))
    pmParam <- setTrDir(pmParam, dataFeatures$trDir)


    if(haveType == "Shiraishi") {
        # Shiraishi-type
        pmParam <- setSigType(pmParam, "independent") # Shiraishi-type
        sigFeatDist <- array(0, dim=c(1, nrow(mutData), ncol(mutData)) )
        
    } else if (haveType == "Alexandrov") {
        # Alexandrov-type
        pmParam <- setSigType(pmParam, "full") # Alexandrov-type
        sigFeatDist <- array(0, dim=c(1, length(mutData), 1) ) 
    }

    sigFeatDist[1,,] <- mutData  # data to visualize
    pmParam <- setSigFeatDist(pmParam, sigFeatDist)

    pmsignature::visPMSignature(pmParam, 1, isScale=TRUE)
}

