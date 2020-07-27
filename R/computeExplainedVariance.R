#' Compute the explained variance.
#'
#' `computeExplainedVariance()` computes for one or more tumors the variance
#' which is explained by the estimated contributions (exposures) of a set of
#' signatures when compared to the observed genomes.
#'
#' @usage computeExplainedVariance(exposures, signatures, genomes)
#' @param exposures (Mandatory) A single vector or list of vectors containing
#' the estimated signature contributions/exposures as provided by the function
#' \code{decomposeTumorGenomes}. A list of vectors is used if the explained
#' variance shall be computed for multiple genomes. The number of exposure
#' vectors must correspond to the number of \code{genomes}. The number
#' of elements of each exposure vector must correspond to the number of
#' \code{signatures}.
#' @param signatures (Mandatory) The list of signatures (vectors, data frames
#' or matrices) for which the exposures were obtained. Each of the list
#' objects represents one mutational signature. Vectors are used for
#' Alexandrov signatures, data frames or matrices for Shiraishi signatures.
#' @param genomes (Mandatory) Can be either a vector, a data frame or a matrix
#' (for an individual tumor genome), or a list of one of these object types
#' (for multiple tumors). Each tumor genome must be of the same form as the
#' \code{signatures}.
#' @return A numeric vector of explained variances, one for each genome.
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics 
#' 20(Suppl 4):152.\cr
#' @seealso \code{\link{decompTumor2Sig}}\cr
#' \code{\link{decomposeTumorGenomes}}\cr
#' \code{\link{plotExplainedVariance}}
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
#' ### compute explained variance for the tumor genomes
#' computeExplainedVariance(exposures, signatures, genomes)
#' 
#' @export computeExplainedVariance
computeExplainedVariance <- function(exposures, signatures, genomes) {

    if (is.probability.object(genomes)) {
        # this is only one genome, use a list nonetheless for later iteration
        genomes <- list(genomes)
    }

    if (!isSignatureSet(genomes)) { # same type of object as signatures
        stop("Parameter genomes must be a set (list) of genomes.")
    }
    
    if (is.probability.vector(exposures)) {
        # same as for genomes
        exposures <- list(exposures)
    }

    if (!isExposureSet(exposures)) {
        stop("Parameter exposures must be a list of probability vectors.")
    }

    if (!isSignatureSet(signatures)) {
        stop("Parameter signatures must be a set (list) of signatures.")
    }

    if (!sameSignatureFormat(signatures, genomes)) {
        stop("Signatures and genomes must be of the same format.")
    }

    # check that we have as many exposure estimates as genomes
    if(length(genomes)!=length(exposures)) {
        stop("Number of exposure vectors different from number of genomes.")
    }

    expvar <- NULL

    # for each genome, compare estimated exposures
    for (g in seq_along(genomes)){

        if (is.null(exposures[[g]])) {
            # can be NULL if minExplainedVariance wasn't reached in a
            # subset search

            expvar <- c(expvar, NA)
            next
        }

        # set NA in exposures to 0
        exposures[[g]][is.na(exposures[[g]])] <- 0
        
        # first compute estimated genome from exposures and signatures
        if (length(signatures) != length(exposures[[g]])) {
            stop("Number of exposures different from number of signatures.")
        }

        first <- TRUE        
        for (s in seq_along(signatures)) {

            if (first) {
                estgenome <- signatures[[s]] * as.numeric(exposures[[g]][s])
                first <- FALSE
            } else {
                estgenome <- estgenome + (signatures[[s]] *
                                              as.numeric(exposures[[g]][s]))
            }
        }

        # now compute the residual sum of squares (RSS) between estimated and
        # observed genome
        #rss <- sum((estgenome - genomes[[g]])^2)
        rss <- computeRSS(estgenome, genomes[[g]])

        # ... and the total sum of squares (TSS) between the observed genome
        # and the "average/unvaried" genome.
        # Important: for the Alexandrov format, the unvaried genome
        # corresponds to the average over all mutation frequencies (which is
        # precisely 1/96 for all mutation frequencies because the 96 frequencies
        # add up to 1!!!)
        # For Shiraishi signatures the "unvaried" genome is represented by 
        # the following exemplary matrix:
        # [ 1/6  1/6  1/6  1/6  1/6  1/6    <- 16.7% for each of 6 base changes
        #   1/4  1/4  1/4  1/4   0    0     <- 25% for each flanking base
        #   ...
        #   1/2  1/2   0    0    0    0 ]   <- 50% for each transcription dir.

        if (isAlexandrovSet(genomes[g])) {

            tss <- sum( (genomes[[g]] - mean(genomes[[g]]))^2 )

        } else if (isShiraishiSet(genomes[g])) {

            unvargenome <- genomes[[g]]
            unvargenome[1,] <- rep(1/6, 6)
            for (ii in seq(2,nrow(genomes[[g]]))) {
                unvargenome[ii,] <- c(rep(1/4, 4), rep(0, 2))
            }
            if ( (nrow(genomes[[g]])%%2) == 0) {
                # even number of rows; last one must be transcription direction
                unvargenome[nrow(genomes[[g]]),] <- c(rep(1/2, 2), rep(0, 4))
            }

            tss <- sum( (genomes[[g]] - unvargenome)^2 )
            
        } else { # something unknown, complain
            stop(paste("genomes must be either of the Alexandrov (numeric",
                       "vectors) or the Shiraishi format (matrices or data",
                       "frames)!"))
        }

        # finally, compute the explained variance
        evar <- 1 - (rss / tss)

        # for the formulas, see for example
        # https://www.rdocumentation.org/packages/SomaticSignatures/
        #                              versions/2.8.3/topics/numberSignatures

        expvar <- c(expvar, evar)
    }

    # finally, name the explained variances like the genomes
    if (!is.null(names(genomes))) {
        names(expvar) <- names(genomes)
    }

    return(expvar)
}

