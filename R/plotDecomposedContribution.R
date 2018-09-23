#' Plot the decomposition (contributions/exposures) of a tumor genome.
#'
#' `plotDecomposedContribution()` plots the decomposition of a tumor genome,
#' i.e., the contributions/exposures obtained from \code{decomposeTumorGenomes}
#' for a set of signatures.
#'
#' @usage plotDecomposedContribution(decomposition, signatures=NULL,
#' removeNA=TRUE)
#' @param decomposition (Mandatory) A decomposition vector (exposure vector)
#' obtained for a single tumor genome.
#' @param signatures (Optional) A list object containing the signatures used
#' to compute the decomposition. If specified, the signature labels used in
#' the plot will be taken from the element names of the list; otherwise
#' signature names will be taken from the exposure object (decomposition) or
#' named from sign_1 to sign_N.
#' @param removeNA (Optional) If \code{TRUE} (default), signatures with
#' an NA as exposure will not be included on the x-axis of the the plot.
#' Exposures can be NA if they have been determined with a greedy search.
#' @return Returns (or draws) a plot of the decomposed tumor genome (i.e.,
#' contributions of the single signatures). 
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
#' ### plot signature composition of the first genome
#' plotDecomposedContribution(exposures[[1]], signatures=NULL)
#' 
#' @import ggplot2
#' @importFrom graphics plot text
#' @export plotDecomposedContribution
plotDecomposedContribution <- function(decomposition, signatures=NULL,
                                       removeNA=TRUE) {

    # dummy declarations to avoid NOTEs by R CMD check together with
    # ggplot2 syntax
    Signatures <- NULL
    Exposures <- NULL
    # end of dummy stuff
    
    if (isExposureSet(decomposition)) {
        if (length(decomposition) > 1) {
            warning(paste("object 'decomposition' passed to",
                          "plotDecomposedContribution is a list with multiple",
                          "elements. Taking only the first element and",
                          "trying to plot"))
        }
        
        decomposition <- decomposition[[1]]
    }

    if (!is.probability.vector(decomposition)) {
        stop("'decomposition' must be a single exposure vector!")
    }

    if (!is.null(signatures) & !isSignatureSet(signatures)) {
        stop("'signatures', if specified, must be a set (list) of signatures!")
    }
    
    # determine signature names

    if (!is.null(signatures) && !is.null(names(signatures))) {
        
        # first choice: explicitly specified signatures
        sigNames <- names(signatures)
        
    } else if (!is.null(names(decomposition))) {
        
        # second guess: directly from the exposure vector
        sigNames <- names(decomposition)
      
    } else {
        
        # last solution: just number them
        sigNames <- paste0("sign_", seq_along(decomposition))

    }

    
    if (length(sigNames) != length(decomposition)) {
        # might happen if signatures aren't associated with the exposure vector
        stop(paste("The number of exposures is different from the",
                   "number of signatures!"))
    }


    if(removeNA) {
        # ignore NAs; might be due to a greedy search
        sigNames <- sigNames[!is.na(decomposition)]
        decomposition <- decomposition[!is.na(decomposition)]
    } else {
        decomposition[is.na(decomposition)] <- 0
    }
    
    # construct data frame for exposures
    df <- data.frame(Signatures=factor(sigNames, levels=sigNames),
                     Exposures=decomposition)

    ggplot(data=df, aes(x=Signatures, y=Exposures)) + xlab(NULL) +
        ylab("Exposures (percent contribution)") +
        geom_bar(stat="identity", width=0.75, color="black", fill="steelblue") +
        geom_text(aes(label=round(Exposures,digits=2)), vjust=1.6,
                  color="white", size=2) +
        theme(panel.background = element_rect(fill="white", colour="black"),
              panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                  colour = "lightgray"),
              panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                  colour = "lightgray"),
              axis.text.x = element_text(face="bold", angle=90, vjust=0.5))

}

