#' Get genomes (mutation frequencies) from \code{MutationFeatureData}.
#'
#' `getGenomesFromMutFeatData()` takes a \code{MutationFeatureData} object
#' (mutation count data) as read by the '\code{pmsignature}' package (e.g.,
#' by \code{pmsignature::readMPFile}, version 0.3.0) and extracts the mutation
#' counts of the genomes therein. For passing the genomes to
#' \code{decomposeTumorGenomes}, the mutation counts must be normalized to
#' mutation frequencies, which is done by default. [IMPORTANT: set
#' \code{normalize} to \code{FALSE} only if you are interested in full
#' integer counts, but do not pass unnormalized counts to
#' \code{decomposeTumorGenomes}!] 
#' 
#' @usage getGenomesFromMutFeatData(mutFeatData, normalize=TRUE)
#' @param mutFeatData (Mandatory) A \code{MutationFeatureData} object as
#' constructed, for example, by \code{pmsignature::readMPFile}.
#' @param normalize (Optional) Boolean value to specify whether to normalize
#' the mutation count data to mutation fractions between 0 and 1. This is
#' the default and NECESSARY in case you want to pass the return value to
#' \code{decomposeTumorGenomes}. Set \code{normalize} to \code{FALSE} only
#' if you are interested in full integer counts, but do not pass unnormalized
#' counts to \code{decomposeTumorGenomes}!
#' @return  A list of mutation frequencies (or mutation counts if not
#' normalized), one object per genome. The format is either according to the
#' Shiraishi or the Alexandrov model, depending on how the mutation data was
#' loaded with \code{pmsignature}. 
#' @author Rosario M. Piro and Sandra Krueger\cr Freie Universitaet Berlin\cr
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
#' ### get breast cancer genomes from 
#' ### Nik-Zainal et al (PMID: 22608084) in the format produced by
#' ### pmsignature (PMID: 26630308)
#' pmsigdata <- system.file("extdata", 
#'          "Nik-Zainal_PMID_22608084-pmsignature-G.Rdata", 
#'          package="decompTumor2Sig")
#' load(pmsigdata)
#' 
#' ### extract the genomes from the pmsignature G object
#' genomes <- getGenomesFromMutFeatData(G, normalize=TRUE)
#' 
#' @importFrom methods is
#' @export getGenomesFromMutFeatData
getGenomesFromMutFeatData <- function(mutFeatData, normalize=TRUE) {
    # convert MutationFeatureData (e.g. read with pmsignature's function
    # readMPFile) to a list of genome matrices

    # IMPORTANT: set normalize=FALSE only if you want to see full counts,
    # but USE ONLY NORMALIZED count tables for predicting/computing exposures,
    # i.e., the contributions of mutational signatures (because these are
    # normalized, too)

    if (!isMutFeatDataObject(mutFeatData)) {
        stop(paste("Function getGenomesFromMutFeatData requires the",
                   "'mutFeatData' input to be compatible an object of type",
                   "MutationFeatureData (as produced by pmsignature's",
                   "function readMPFile; version 0.3.0)"))
    }
    
    tableList <- list()

    # get list of samples (numbers)
    samples <- unique(getCountData(mutFeatData)[2,])

    # get number of flanking bases and transcription direction
    numBases <- getNumFlankingBases(mutFeatData)
    trDir <- haveTrDir(mutFeatData)

    
    # what signature model do we have?
    isShiraishi <- TRUE
    if (getSigType(mutFeatData) == "full") {
        isShiraishi <- FALSE
    }

    if (isShiraishi) {
        # if the user has set trDir = TRUE, we need flankinBasesNum+1 rows,
        # otherwise only the flankingBasesNum
        if (haveTrDir(mutFeatData)) {
            rows <- as.numeric(getNumFlankingBases(mutFeatData)+1)
        } else {
            rows <- as.numeric(getNumFlankingBases(mutFeatData)) 
        }
    } else {
        rows <- 1   # for Alexandrov we expect only vectors, not matrices
    }

    for (l in seq_along(samples)) {  # for each sample l

        # initialize the new count table/matrix (or vector) with 0
        if (isShiraishi) {
            countTable <- matrix(rep(0, 6*rows), nrow=rows, ncol=6) 
        } else {
            countTable <- rep(0, 6 * 4^(numBases-1) * (as.numeric(trDir)+1))

        }
        
        # in oneSample we extract the count data relative to sample l
        oneSample <-
            getCountData(mutFeatData)[,which(getCountData(mutFeatData)[2,]==l)]

        for (i in seq_len(ncol(oneSample))) { # for each i in the vector

            # index of the mutation type in featureVectorList
            mutTypeIndex <- oneSample[1,i] 

            # get the features for this mutation type
            # (corresponds to the column indices to be incremented for each
            # row in the count table)
            v <- getFeatVectList(mutFeatData)[,mutTypeIndex] 

            if (isShiraishi) {
                for (k in seq_len(rows)) { #for each k in the length of our rows
                
                    # in our countTable at position k and v[k] we store the
                    # number of occurrences which is stored in the third row
                    # and i-th column of oneSample
                    countTable[k,v[k]] <- countTable[k,v[k]] + oneSample[3,i]
                }
            } else {
                # only one row, so v is already the index and not a list
                # of indices
                countTable[v] <- countTable[v] + oneSample[3,i]
            }
        }

        # normalize the count table/vector?
        if (normalize) {
            if (isShiraishi) {
                countTable <- countTable/(sum(countTable)/nrow(countTable))
            } else {
                countTable <- countTable/sum(countTable)
            }
        }

        # finally, name the rows and/or columns of the count table/vector
        if (isShiraishi) {
            countTable <- setNames4ShiraishiTable(countTable)
        } else {
            names(countTable) <-
                buildSortedAlexandrovSignaturePatternList(numBases=numBases,
                                                          trDir=trDir)
        }
        
        tableList[[l]] <- countTable  # put each countTable in the tableList

    }

    # keep the samples' names/IDs!
    names(tableList) <- getSampleList(mutFeatData)

    return(tableList)
}

