#' Compose tumor genomes from exposures.
#'
#' `composeGenomesFromExposures()` re-composes (or predicts) tumor genomes
#' (i.e., their mutation frequencies) from the given mutational signatures and
#' their corresponding exposures, or contributions. The (re-)composition is
#' performed by computing the weighted sum of the mutational signatures, where
#' the weights are the exposures (=contributions) of the corresponding
#' signatures. This can, for example, be used to verify that a decomposition
#' obtained from \code{decomposeTumorGenomes} is meaningful.
#'
#' @usage composeGenomesFromExposures(exposures, signatures)
#' @param exposures (Mandatory) A single vector or list of vectors containing
#' the estimated signature contributions/exposures as computed by the function
#' \code{decomposeTumorGenomes}. A list of vectors is used if the
#' (re-)composition shall be performed for multiple genomes. The number of
#' elements of each exposure vector must correspond to the number of
#' \code{signatures}. 
#' @param signatures (Mandatory) The list of signatures (vectors, data frames
#' or matrices) for which the exposures were obtained. Each of the list
#' objects represents one mutational signature. Vectors are used for
#' Alexandrov signatures, data frames or matrices for Shiraishi signatures.
#' @return A list of "predicted" genomes, i.e., the frequencies of their
#' mutational patterns computed as weighted sums of the mutational signatures,
#' where the weights correspond to the contributions of, i.e., exposures to,
#' the corresponding signatures.
#'
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics 
#' 20(Suppl 4):152.\cr
#' @seealso \code{\link{decompTumor2Sig}}\cr
#' \code{\link{decomposeTumorGenomes}}
#' @examples
#' 
#' ### get Alexandrov signatures from COSMIC
#' signatures <- readAlexandrovSignatures()
#' 
#' ### load preprocessed breast cancer genomes (object 'genomes') from
#' ### Nik-Zainal et al (PMID: 22608084) 
#' gfile <- system.file("extdata",
#'          "Nik-Zainal_PMID_22608084-genomes-Alexandrov_3bases.Rdata", 
#'          package="decompTumor2Sig")
#' load(gfile)
#' 
#' ### compute exposures
#' exposures <- decomposeTumorGenomes(genomes, signatures, verbose=FALSE)
#' 
#' ### re-compose (predict) tumor genome features from exposures
#' predGenomes <- composeGenomesFromExposures(exposures, signatures)
#' 
#' @export composeGenomesFromExposures
composeGenomesFromExposures <- function(exposures, signatures) {

    if (is.probability.vector(exposures)) {
        # this is only one genome, use a list nonetheless for later iteration
        exposures <- list(exposures)
    }

    if (!isExposureSet(exposures)) {
        stop("Parameter exposures must be a list of probability vectors.")
    }
    
    if (!isSignatureSet(signatures)) {
        stop("Parameter signatures must be a set (list) of signatures.")
    }

    predGenomes <- list()
    
    for (e in seq_along(exposures)) {

        if(is.null(exposures[[e]])) {
            # didn't exceed minExplainedVariance ... can't reconstruct this
            predGenomes[[e]] <- NULL
            next
        }
        
        exp <- exposures[[e]]

        if (!is.numeric(exp)) {
            stop("exposures must be numeric vectors.")
        }

        
        # make sure we replace NA by 0 (signatures that were not selected in a
        # search, i.e., they of course have 0 contribution/exposure
        exp[which(is.na(exp))] <- 0
        
        # check that we have as many exposure values as signatures
        if(length(exp)!=length(signatures)) {
            stop(paste("Number of exposure values different from",
                       "number of signatures."))
        }

        # multiply exposures with the corresponding signatures and sum results
        for (ii in seq_along(exp)) {

            if (ii==1) {
                # first: set pred
                pred <- exp[ii]*signatures[[ii]]
            } else {
                # add the following
                pred <- pred + exp[ii]*signatures[[ii]]
            }
        }

        tnbt <- determineTypeNumBasesAndTrDir(pred)
        if (tnbt$type == "Shiraishi") {
            # Shiraishi model; add column and row names
            pred <- setNames4ShiraishiTable(pred)
        } else if (tnbt$type == "Alexandrov") {
            # Alexandrov model; add vector names
            names(pred) <-
                buildSortedAlexandrovSignaturePatternList(tnbt$numBases,
                                                          tnbt$trDir)
        }
        
        predGenomes[[e]] <- pred
    }

    # names of the genomes
    if (!is.null(names(exposures))) {
        names(predGenomes) <- names(exposures)
    }
    
    return(predGenomes)
}

