#' Get genomes (mutation frequencies) from \code{MutationFeatureData}.
#'
#' `getGenomesFromMutFeatData()` takes a \code{MutationFeatureData} object
#' (mutation count data) as read by the '\code{pmsignature}' package (e.g.,
#' by \code{pmsignature::readMPFile}) and extracts the mutation counts of the
#' genomes therein. For passing the genomes to \code{decomposeTumorGenomes},
#' the mutation counts must be normalized to mutation frequencies, which is
#' done by default. [IMPORTANT: set \code{normalize} to \code{FALSE} only if
#' you are interested in full integer counts, but do not pass unnormalized
#' counts to \code{decomposeTumorGenomes}!] 
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
#' @return  A list of (normalized) mutation counts, one object per genome. The
#' format is the same table used by the corresponding Shiraishi signatures. 
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
#' ### get breast cancer genomes from 
#' ### Nik-Zainal et al (PMID: 22608084) in the format produced by
#' ### pmsignature (PMID: 26630308)
#' pmsigdata <- system.file("extdata", 
#'          "Nik-Zainal_PMID_22608084-pmsignature-G.Rdata", 
#'          package="decompTumor2Sig")
#' load(pmsigdata)
#' 
#' ### extract the genomes from the pmsignature G object
#' \donttest{
#' genomes <- getGenomesFromMutFeatData(G, normalize=TRUE)
#' }
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

    if (!requireNamespace("pmsignature", quietly=TRUE)) {
        stop(paste("Function getGenomesFromMutFeatData requires the",
                   "package pmsignature to be installed."))
    }
    
    if (!is(mutFeatData, "MutationFeatureData")) {
        stop(paste("mutFeatData must be an object of type MutationFeatureData",
                   "(as produced by pmsignature's readMPFile)"))
    }
    
    # if the user has set trDir = TRUE, we need flankinBasesNum+1 rows,
    # otherwise only the flankingBasesNum
    if (haveTrDir(mutFeatData)) {
        rows <- as.numeric(getNumFlankingBases(mutFeatData)+1)
    } else {
        rows <- as.numeric(getNumFlankingBases(mutFeatData)) 
    }

    tableList <- list()

    # get list of samples (numbers)
    samples <- unique(getCountData(mutFeatData)[2,])   

    for (l in seq_along(samples)) {  # for each sample l

        # initialize the new count table/matrix with 0
        countTable <- matrix(rep(0, 6*rows), nrow=rows, ncol=6) 

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

            
            for (k in seq_len(rows)) { #for each k in the length of our rows
                
                # in our countTable at position k and v[k] we store the
                # number of occurrences which is stored in the third row
                # and i-th column of oneSample
                countTable[k,v[k]]= countTable[k,v[k]] + oneSample[3,i] 
            }
        }

        # normalize the count table?
        if (normalize) {
            countTable <- countTable/(sum(countTable)/nrow(countTable))
        }
        
        tableList[[l]] <- countTable  # put each countTable in the tableList

    }

    # keep the samples' names/IDs!
    names(tableList) <- getSampleList(mutFeatData)

    return(tableList)
}

