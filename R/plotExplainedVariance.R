#' Plot the explained variance as a function of the number of signatures
#'
#' `plotExplainedVariance()` plots the explained variance of a single tumor
#' genome's mutation patterns as a function of the number of signatures
#' (increasing subsets of signatures) used for decomposition. For each
#' number K of signatures, the highest variance explained by possible
#' subsets of K signatures will be plotted (full or greedy search, see below).
#' This can help to evaluate what minimum threshold for the explained variance
#' can be used to decompose tumor genomes with the function
#' \code{decomposeTumorGenomes}.
#'
#' @usage plotExplainedVariance(genome, signatures, minExplainedVariance=NULL,
#' minNumSignatures=2, maxNumSignatures=NULL, greedySearch=FALSE)
#' @param genome (Mandatory) The mutation load of a single genome in
#' Alexandrov- of Shiraishi-format, i.e. as vector or matrix. The format
#' must be the same as the one used for the \code{signatures} (see below).
#' @param signatures (Mandatory) The list of signatures (vectors,
#' data frames or matrices) which are to be evaluated. Each of the list
#' objects represents one mutational signature. Vectors are used for
#' Alexandrov signatures, data frames or matrices for Shiraishi signatures.
#' @param minExplainedVariance (Optional) If a numeric value between 0 and 1
#' is specified, the plot highlights the smallest subset of signatures which
#' is sufficient to explain at least the specified fraction of the variance
#' of the genome's mutation patterns. If, for example,
#' \code{minExplainedVariance} is 0.99 the smallest subset of signatures
#' that explains at least 99\% of the variance will be highlighted.
#' @param minNumSignatures (Optional) The plot will be generated only for
#' K>=\code{minNumSignatures}. 
#' @param maxNumSignatures (Optional) The plot will be generated only for
#' K<=\code{maxNumSignatures}. 
#' @param greedySearch (Optional) If \code{greedySearch} is set to \code{TRUE}
#' then not all possible combinations of \code{minNumSignatures} to
#' \code{maxNumSignatures} signatures will be checked. Instead, first all
#' possible combinations for exactly \code{minNumSignatures} will be checked
#' to select the best starting set, then iteratively the next best signature
#' will be added (maximum increase in explained variability) until
#' \code{maxNumSignatures} is reached). NOTE: while this is only an
#' approximation, it is highly recommended for large sets of signatures (>15)!
#' @return Returns (or draws) a plot of the explained variance as a function
#' of the number of signatures.
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @seealso \code{\link{decompTumor2Sig}}\cr
#' \code{\link{decomposeTumorGenomes}}\cr
#' \code{\link{computeExplainedVariance}}
#' @examples
#' 
#' ### get 15 pre-processed Shiraishi signatures computed (object 'signatures') 
#' ### from 435 tumor genomes Alexandrov et al (PMID: 23945592)
#' ### using the pmsignature package
#' sfile <- system.file("extdata",
#'          "Alexandrov_PMID_23945592_435_tumors-pmsignature-15sig.Rdata", 
#'          package="decompTumor2Sig")
#' load(sfile)
#' 
#' ### load preprocessed breast cancer genomes (object 'genomes') from
#' ### Nik-Zainal et al (PMID: 22608084) 
#' gfile <- system.file("extdata",
#'          "Nik-Zainal_PMID_22608084-genomes-Shiraishi_5bases_trDir.Rdata", 
#'          package="decompTumor2Sig")
#' load(gfile)
#' 
#' ### plot the explained variance for 2 to 6 signatures of the first genome
#' plotExplainedVariance(genomes[[1]], signatures,
#'          minExplainedVariance=0.98, minNumSignatures=2, maxNumSignatures=6)
#' 
#' @importFrom graphics abline plot points text
#' @export plotExplainedVariance
plotExplainedVariance <- function(genome, signatures,
                                  minExplainedVariance=NULL,
                                  minNumSignatures=2, maxNumSignatures=NULL,
                                  greedySearch=FALSE) {
    # input: gnome = mutation counts for a single genome as matrix or vector
    #        signatures = a list of signatures (matrices or vectors)
    #        [need to have the same format]

    if (!isSignatureSet(signatures)) {
        stop("Parameter 'signatures' must be a set (list) of signatures!")
    }

    if (isSignatureSet(genome)) { # it's a list of genomes
        if (length(genome) == 1) {  # accept if only one!
            genome <- genome[[1]]
        } else { # more than one genome
            stop(paste("plotExplainedVariance can plot the explained",
                       "variance for only one genome!"))
        }
    }

    if (!is.probability.object(genome)) {
        stop("'genome' must be a genome in Alexandrov or Shiraishi format!")
    }


    # check the genome format; same as signature format?
    if (!sameSignatureFormat(list(genome), signatures)) {
        stop("Formats of genome and signatures must match!")
    }

    # is the signatures are unnamed, name them by enumerating them
    if(is.null(names(signatures))) {
        names(signatures) <- paste0("sign_",seq_along(signatures))
    }

    # if maximum number of signatures is not defined, set it to the total
    # number of signatures
    if(is.null(maxNumSignatures)) {
        maxNumSignatures <- length(signatures)
    }

    # greedySearch must be logical
    if (!is.logical(greedySearch)) {
         stop("greedySearch must be logical (TRUE or FALSE)!")
    }


    bestDecompositions <- list()
    
    # iterating for all possible numbers of signatures from minNumSignatures
    # to maxNumSignatures

    if(!greedySearch) {
        # not greedy; test ALL combinations

        for (k in seq(minNumSignatures,maxNumSignatures)) {
    
            bestDecompositions[[length(bestDecompositions)+1]] <-
                getBestDecomp4Ksignatures(genome, signatures, k)
        }
    } else {
        # greedy; test all combinations of minNumSignatures, then increase

        k <- minNumSignatures

        bestDecompositions[[length(bestDecompositions)+1]] <-
            getBestDecomp4Ksignatures(genome, signatures, k)

        while(k < maxNumSignatures) {
            # need to add another signature

            haveSubset <-
                bestDecompositions[[length(bestDecompositions)]]$sigList

            bestDecompositions[[length(bestDecompositions)+1]] <-
                addBestSignatureToSubset(genome, signatures, haveSubset)

            k <- k + 1
        }
    }

    # now plot: 
    # k on x-axis; 
    # expl. var on y-axis; 
    # min. k with expl. var >= thres in red

    k <- vapply(bestDecompositions, function(x) { x$k }, FUN.VALUE=numeric(1) )
    explVar <- vapply(bestDecompositions, function(x) { x$explVar },
                      FUN.VALUE=numeric(1) )
    sigList <- lapply(bestDecompositions, function(x) { x$sigList } )  # list

    minExplVarIndex <- NULL
    if (!is.null(minExplainedVariance)) {
        whichExceed <- which(explVar >= minExplainedVariance)
        if (length(whichExceed) > 0) {
            minExplVarIndex <- whichExceed[1]  # first to exceed the threshold
        }
    }

    yrange <- c(max(min(explVar),0), 1)
    
    # basic plot
    plot(k, explVar, xlab="number of signatures",
         ylab="highest explained variance", ylim=yrange)
    abline(h=1, lty=2, col="black")
    
    # indicate threshold for minimum explained variance?
    if (!is.null(minExplVarIndex)) {
        
        # indicate desired threshold for explained variance
        abline(h=minExplainedVariance, lty=2, col="red")
        y4Label <- minExplainedVariance - 2*(max(explVar)-min(explVar))/100
        text(c(maxNumSignatures), c(y4Label),
             paste("min. explained variance: ",minExplainedVariance),
             cex=0.6, pos=2, col="red")

        # also highlight the signatures exceeding the threshold
        points(c(k[minExplVarIndex]), c(explVar[minExplVarIndex]), bg="red",
               col="black", pch=21)

        text(c(k[minExplVarIndex]),
             c(min(explVar)+(max(explVar)-min(explVar))/2),
             paste(sigList[[minExplVarIndex]], collapse="\n"), cex=0.6,
             pos=4, col="red")
        
        abline(v=k[minExplVarIndex], lty=3, col="red")
    }
}

