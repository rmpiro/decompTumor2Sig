####################
# public functions #
####################


#' Decompose tumor genomes into mutational signatures
#'
#' `decomposeTumorGenomes()` is the core function of this package. It
#' decomposes tumor genomes into a given set of mutational signatures by
#' computing their contributions (exposures) to the mutational load via
#' quadratic programming. The function takes a set of mutational signatures
#' and the mutation features of one or more tumor genomes and computes
#' weights, i.e., contributions for each of the signatures in each
#' individual genome. Alternatively, the function can determine for each
#' genome only a subset of signatures whose contributions are sufficient
#' to exceed a user-given minimum threshold for the explained variance
#' of the genome's mutation patterns.
#'
#' @usage decomposeTumorGenomes(genomes, signatures, minExplainedVariance=NULL,
#' minNumSignatures=2, maxNumSignatures=NULL, greedySearch=FALSE,
#' constrainToMaxContribution=FALSE, tolerance=0.1, verbose=FALSE)
#' @param genomes (Mandatory) Can be either a vector, a data frame or a
#' matrix (for an individual tumor genome), or a list of one of these
#' object types (for multiple tumors). Each tumor genome must be of the
#' same form as the \code{signatures}.
#' @param signatures (Mandatory) A list of vectors, data frames or matrices.
#' Each of the objects represents one mutational signature. Vectors are
#' used for Alexandrov signatures, data frames or matrices for Shiraishi
#' signatures. 
#' @param minExplainedVariance (Optional) If \code{NULL} (default), exactly
#' \code{maxNumSignatures} (see below; default: all) will be taken for
#' decomposing each genome. If a numeric value between 0 and 1 is specified
#' for \code{minExplainedVariance}, for each genome the function will select
#' the smallest number of signatures which is sufficient to explain at least
#' the specified fraction of the variance of the genome's mutation patterns.
#' E.g., if \code{minExplainedVariance}=0.99 the smallest subset of
#' signatures that explains at least 99\% of the variance is taken.
#' Please note: depending on the number of signatures, this may take quite
#' a while because by default for each number K of signatures, all possible
#' subsets composed of K signatures will be tested to identify the subset that
#' explains the highest part of the variance. If not enough variance is
#' explained, K will be incremented by one. Notes: 1) to speed up the search,
#' the parameters \code{minNumSignatures}, \code{maxNumSignatures} and
#' \code{greedySearch} can be used; 2) for genomes for which
#' none of the possible subsets of signatures explains enough variance, the
#' returned exposure vector will be set to \code{NULL}. 
#' @param minNumSignatures (Optional) Used if \code{minExplainedVariance} is
#' specified (see above). To find the smallest subset of signatures which
#' explain the variance, at least this number of signatures will be taken. This
#' can be used to reduce the search space in a time-consuming search over a
#' large number of signatures.
#' @param maxNumSignatures (Optional) If \code{minExplainedVariance} is
#' specified to find the smallest subset of signatures which
#' explain the variance, at most \code{maxNumSignatures} will be taken. This 
#' can be used to reduce the search space in a time-consuming search over a
#' large number of signatures. If \code{minExplainedVariance} is \code{NULL},
#' then exactly \code{maxNumSignatures} signatures will be used. The default
#' for \code{maxNumSignatures} is \code{NULL} (all signatures).
#' @param greedySearch (Optional) Used only in case \code{minExplainedVariance}
#' has been specified. If \code{greedySearch} is \code{TRUE} then not all
#' possible combinations of \code{minNumSignatures} to \code{maxNumSignatures}
#' signatures will be checked. Instead, first all possible combinations for
#' exactly \code{minNumSignatures} will be checked to select the best starting
#' set, then iteratively the next best signature will be added (maximum
#' increase in explained variability) until \code{minExplainedVariance} of the
#' variance can be explained (or \code{maxNumSignatures} is exceeded).
#' NOTE: this approximate search is highly recommended for large sets of
#' signatures (>15)! 
#' @param constrainToMaxContribution (Optional) [Note: this is EXPERIMENTAL
#' and is usually not needed!] If \code{TRUE}, the maximum contribution that
#' can be attributed to a signature will be constraint by the variant feature
#' counts (e.g., specific flanking bases) observed in the individual tumor
#' genome. If, for example, 30\% of all observed variants have a specific
#' feature and 60\% of the variants produced by a mutational process/signature
#' will manifest the feature, then the signature can have contributed up to
#' 0.3/0.6 (=0.5 or 50\%) of the observed variants. The lowest possible
#' contribution over all signature features will be taken as the allowed
#' maximum contribution of the signature. This allowed maximum will
#' additionally be increased by the value specified as \code{tolerance}
#' (see below). For the illustrated example and \code{tolerance}=0.1 a
#' contribution of up to 0.5+0.1 = 0.6 (or 60\%) of the signature would be
#' allowed. 
#' @param tolerance (Optional) If \code{constrainToMaxContribution} is
#' \code{TRUE}, the maximum contribution computed for a signature is increased
#' by this value (see above). If the parameter \code{constrainToMaxContribution}
#' is \code{FALSE}, the tolerance value is ignored. Default: 0.1. 
#' @param verbose (Optional) If \code{TRUE} some information about the
#' processed genome and used number of signatures will be printed.
#' @return A list of signature weight vectors (also called 'exposures'), one
#' for each tumor genome. E.g., the first vector element of the first list
#' object is the weight/contribution of the first signature to the first
#' tumor genome. IMPORTANT: If \code{minExplainedVariance} is specified, then
#' the exposures of a genome will NOT be returned if the minimum explained
#' variance is not reached within the requested minimum and maximum numbers
#' of signatures (\code{minNumSignatures} and \code{maxNumSignatures})! The
#' corresponding exposure vector will be set to \code{NULL}. 
#' @author Rosario M. Piro, Sandra Krueger\cr Freie Universitaet Berlin\cr
#' Maintainer: Rosario M. Piro\cr E-Mail: <rmpiro@@gmail.com> or
#' <r.piro@@fu-berlin.de>
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
#' ### get Alexandrov signatures from COSMIC
#' signatures <- readAlexandrovSignatures()
#' 
#' ### load reference genome
#' refGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' 
#' ### read breast cancer genomes from Nik-Zainal et al (PMID: 22608084) 
#' gfile <- system.file("extdata",
#'          "Nik-Zainal_PMID_22608084-VCF-convertedfromMPF.vcf.gz", 
#'          package="decompTumor2Sig")
#' genomes <- readGenomesFromVCF(gfile, numBases=3, type="Alexandrov",
#'          trDir=FALSE, refGenome=refGenome, verbose=FALSE)
#' 
#' ### compute exposures
#' exposures <- decomposeTumorGenomes(genomes, signatures, verbose=FALSE)
#' 
#' ### (for further examples on searching subsets, please see the vignette)
#' 
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @export decomposeTumorGenomes
decomposeTumorGenomes <- function(genomes, signatures,
                                  minExplainedVariance=NULL,
                                  minNumSignatures=2, maxNumSignatures=NULL,
                                  greedySearch=FALSE, 
                                  constrainToMaxContribution=FALSE,
                                  tolerance=0.1, verbose=FALSE) {
    
    # input: gnomes = either a list of genome mutation counts as matrices or
    #                  vectors, or a single genome (matrix or vector)
    #        signatures = a list of signatures (matrices or vectors)
    #        [need to have the same format]

    # constrainToMaxContribution and tolerance:
    # We want to constrain the maximum contribution each signature can have
    # (given some tolerance): For each row in sigMa (each signature), we
    # divide the counts from the genome by the corresponding fractions in
    # the signature and take the minimum over the signature, the signature
    # can't have a higher contribution than this ratio, but we add a tolerance
    # value (given that the data is noisy!)
    # Example: if in the genome 30% of variants have a specific feature
    # (e.g., a specific flanking base), i.e., the counts score is 0.3, and
    # 60% (0.6) of the variants produced by a process/signature have this 
    # feature, then the signature can have contributed up to 0.3/0.6 = 0.5 (50%)
    # of the genome's variants. 
    # With the default tolerance of 0.1, we would allow the quadratic
    # programming approach to assign up to 0.5+0.1 (60%) of the variants to
    # this signature

    # minExplainedVariance and minNumSignatures and maxNumSignatures:
    # If minExplainedVariance is specified, the minimum number of signatures
    # required will be determined such that the minimum threshold for
    # explained variance is satisfied. The tested numbers of signatures range
    # from minNumSignatures (default: 2) to maxNumSignatures (default: all)

    # if greedySearch=TRUE then _not_ all possible combinations of
    # minNumSignatures to maxNumSignatures signatures will be checked. Instead,
    # first all possible combinations for exactly minNumSignatures will be
    # checked, then the next best signature will be added (maximum increase
    # in explained variability) until minExplainedVariance is reached
    # (or maxNumSignatures is exceeded).
    
    # If minExplainedVariance is NULL, then exactly maxNumSignatures
    # signatures will be taken (default: all)

    # Returns: exposure(s), but only if they satisfy minExplainedVariance or
    # if not minimum explained variance is requested.

    
    if (!isSignatureSet(signatures)) {
        stop("Parameter 'signatures' must be a set (list) of signatures.")
    }

    # if maximum number of signatures is not defined, set it to the total
    # number of signatures
    if(is.null(maxNumSignatures)) {
        maxNumSignatures <- length(signatures)
    }

    # if we don't require to find the minimum number of signatures for which
    # we exceed the minimum explained variance, we simply take the maximum
    # number of signatures (usually the default: all)
    # This way we make sure that we use exactly k = maxNumSignatures
    if (is.null(minExplainedVariance)) {
        minNumSignatures <- maxNumSignatures
        
    } else if(!is.numeric(minExplainedVariance)
              || minExplainedVariance<0 || minExplainedVariance>1) {

        stop("minExplainedVariance must be NULL or between 0 and 1!")
    }

    
    # is the signatures are unnamed, name them by enumerating them
    if(is.null(names(signatures))) {
        names(signatures) <- paste0("sign_",seq_along(signatures))
    }

    
    if (!is.logical(greedySearch)) {
         stop("greedySearch must be logical (TRUE or FALSE)!")
    }

    if (!is.logical(constrainToMaxContribution)) {
         stop("constrainToMaxContribution must be logical (TRUE or FALSE)!")
    }

    if (constrainToMaxContribution && (tolerance < 0 ||  tolerance > 1)) {
        stop(paste("tolerance must be between 0 and 1 when constraining the",
                   "maximum contribution of signatures (default: 0.1)!"))
    }



    if (is.probability.object(genomes)) {
        # this is only one genome, use a list nonetheless for later iteration
        genomes <- list(genomes)
    }

    if (!isSignatureSet(genomes)) { # same type of object as signatures
        stop("Parameter 'genomes' must be a set (list) of genomes.")
    }
    
    if (is.null(names(genomes))) { # numbering of genomes if they are not named!
                                   # (we need the names for the results)
        names(genomes) <- paste0("genome_", as.character(seq_along(genomes)))
    }


    if (!sameSignatureFormat(signatures, genomes)) {
        stop("Signatures and genomes must be of the same format.")
    }

    
    # initialize an empty list with all elements set to NULL
    decompositions  <- vector("list", length(genomes))  
    
    for (g in seq_along(genomes)) {
        # evaluate each genome individually

        if(verbose) {
            cat(paste0("Decomposing genome ",g," (",names(genomes)[g],")"))
        }
        
        counts <- genomes[[g]]
       
        # iterating for all possible numbers of signatures from
        # minNumSignatures to maxNumSignatures
        decompTmpList <- list()

        if(!greedySearch) {
            # default: full search; all possible combinations!
            for (k in seq(minNumSignatures, maxNumSignatures)) {

                if(verbose) {
                    cat(paste(" with", k, "signatures ...\n"))
                }
            
                decompTmpList[[length(decompTmpList)+1]] <-
                    getBestDecomp4Ksignatures(counts, signatures, k,
                                              constrainToMaxContribution,
                                              tolerance)

                if (!is.null(minExplainedVariance)) {
                    # test the decomposition to check if we found one that
                    # explains at least minExplainedVariance of the variance
                    # of the observed mutation data of the genome

                    if (decompTmpList[[length(decompTmpList)]]$explVar
                        >= minExplainedVariance) {
                        
                        break
                    }
                }
            }
        } else {
            # greedy search! start with best combination of minNumSignatures;
            # then add one at a time
            k <- minNumSignatures

            if(verbose) {
                cat(paste(" with", k, "signatures ...\n"))
            }

            decompTmpList[[length(decompTmpList)+1]] <-
                getBestDecomp4Ksignatures(counts, signatures, k,
                                          constrainToMaxContribution,
                                          tolerance)

            while(decompTmpList[[length(decompTmpList)]]$explVar
                       < minExplainedVariance
                  && k < maxNumSignatures) {

                if(verbose) {
                    cat(paste(" adding signature", k+1, "...\n"))
                }

                haveSubset <- decompTmpList[[length(decompTmpList)]]$sigList
                decompTmpList[[length(decompTmpList)+1]] <-
                    addBestSignatureToSubset(counts, signatures, haveSubset,
                                             constrainToMaxContribution,
                                             tolerance)
                k <- k + 1
            }
        }

        haveExplVar <- decompTmpList[[length(decompTmpList)]]$explVar
        
        if(verbose) {
            cat(paste(" explained variance:", haveExplVar))
        }
        
        # The last decomposition we obtained should be the correct one, either
        # because the k signatures exceed the threshold minExplainedVariance,
        # or because no such threshold was required
        # However, if no combination of signatures reached the specified
        # threshold, we do not report the result!

        if (is.null(minExplainedVariance)
            || haveExplVar >= minExplainedVariance) {
            # threshold not specified or satisfied

            # return _all_ signatures but set unused signatures to NA

            # empty vector of NAs
            decompositions[[g]] <- rep(NA, length(signatures))
            names(decompositions[[g]]) <- names(signatures)

            # set used signatures to their exposure/contribution
            listlen <- length(decompTmpList)
            decompositions[[g]][decompTmpList[[listlen]]$sigList] <-
                decompTmpList[[listlen]]$decomposition

            if(verbose) {
                cat("\n")
            }
        } else {
            # threshold specified but not satisfied; do not report the results
            # (set to NULL instead)

            if(verbose) {
                cat(paste(" <", minExplainedVariance,"(rejected)\n"))
            }
        }

        names(decompositions)[g] <- names(genomes)[g]

    }
    
    return(decompositions)
}

