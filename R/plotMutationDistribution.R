#' Plot mutation frequency data of a mutational signature or tumor genome.
#'
#' `plotMutationDistribution()` plots a single signature or the mutation
#' frequency data for a single genome. This works for signatures or genome
#' data of both the Shiraishi and the Alexandrov type.
#'
#' @usage plotMutationDistribution(mutData, colors = NULL, strip = NULL)
#' @param mutData (Mandatory) The signature or genome mutation frequency data
#' to be plotted. This can either be a matrix (Shiraishi model) or a
#' numeric vector (Alexandrov model).
#' @param colors Vector of colors to be used for the base change data. For
#' Alexandrov-type data, this vector must contain six elements (one per base
#' change). For Shiraishi-type data, this vector must contain four elements
#' (one per base). If \code{NULL} (default), for Alexandrov-type data, the
#' colors are set to those used by the COSMIC website; for Shiraishi-type
#' data, the consensus base colors for sequence logos will be used.
#' @param strip Background color for strip labels; used only for
#' Alexandrov-type data. If \code{NULL} (default), "papayawhip" will be used.
#' @return Returns (or draws) a plot according to the Alexandrov or Shiraishi
#' model of mutational signatures.
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2018) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics (accepted for
#' publication).\cr
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
#' plotMutationDistribution(signatures[[1]])
#' 
#' ### read four Shiraishi signatures for breast cancer genomes from 
#' ### Nik-Zainal et al (PMID: 22608084) from flat files
#' sigfiles <- system.file("extdata",
#'          paste0("Nik-Zainal_PMID_22608084-pmsignature-sig",1:4,".tsv"),
#'          package="decompTumor2Sig")
#' signatures <- readShiraishiSignatures(sigfiles)
#' 
#' ### plot the first Shiraishi signature
#' plotMutationDistribution(signatures[[1]])
#' 
#' @importFrom methods new
#' @export plotMutationDistribution
plotMutationDistribution <- function(mutData, colors = NULL, strip = NULL) {

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

    
    # try to determine the signature model, the number of bases, and
    # whether the transcription direction was considered
    dataFeatures <- determineTypeNumBasesAndTrDir(mutData)
    if (is.null(dataFeatures$type)) {
        stop(paste("Object 'mutData' must be a single Shiraishi- or",
                   "Alexandrov signature!"))
    }

    sigType <- dataFeatures$type
    numBases <- dataFeatures$numBases
    trDir <- dataFeatures$trDir
    
    if(sigType == "Shiraishi") {
        # Shiraishi-type

        plotShiraishiModel(mutData, numBases, trDir, colors)
        
    } else if (sigType == "Alexandrov") {
        # Alexandrov-type

        plotAlexandrovModel(mutData, numBases, trDir, colors, strip)
    }

}

