#' Evaluate tumor decomposition quality.
#'
#' `evaluateDecompositionQuality()` evaluates the quality of the decomposition
#' into exposures of a single tumor. The function evaluates the quality of a
#' decomposition obtained from the function \code{decomposeTumorGenomes}
#' by comparing the re-composed (=re-constructed) tumor genome mutation
#' frequencies to those actually observed in the tumor genome. Tumor genome
#' mutation frequencies are reconstructed using
#' \code{composeGenomesFromExposures} and the results can optionally be plotted.
#'
#' @usage evaluateDecompositionQuality(exposure, signatures, genome,
#' plot=FALSE)
#' @param exposure (Mandatory) A single vector containing the estimated
#' signature contributions, or exposures, of a single tumor as provided by
#' \code{decomposeTumorGenomes}. The number of elements of the
#' exposure vector must correspond to the number of signatures (see below). 
#' @param signatures (Mandatory) The list of signatures (vectors, data
#' frames or matrices) for which the exposures were obtained. Each of the
#' list objects represents one mutational signature. Vectors are used for
#' Alexandrov signatures, data frames or matrices for Shiraishi signatures.
#' @param genome (Mandatory) A single tumor genome in form of mutation
#' frequencies specified either in the Alexandrov or the Shiraishi format
#' (must match the format used for \code{signatures}, see above). 
#' @param plot (Optional) If \code{FALSE} (default), the numerical results
#' (see below) will be returned. If \code{TRUE}, the reconstructed mutation
#' frequencies will be plotted against the original, observed mutation
#' frequencies and the numerical results will be integrated as text labels
#' in the plot.
#' @return A named list object containing measurements for the Pearson
#' correlation coefficient between the reconstructed and observed mutation
#' frequencies, and the explained variance; or alternatively, a plot with
#' these measurements (see option \code{plot} above).
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @author Rosario M. Piro, Freie Universitaet Berlin, \email{rmpiro@@gmail.com}
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @seealso \code{\link{decompTumor2Sig}}\cr
#' \code{\link{decomposeTumorGenomes}}\cr
#' \code{\link{composeGenomesFromExposures}}\cr
#' \code{\link{computeExplainedVariance}}
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
#' ### evaluate the decomposition by comparing to the original data
#' evaluateDecompositionQuality(exposures[[1]], signatures, genomes[[1]])
#' 
#' @importFrom stats cor scatter.smooth
#' @importFrom graphics text
#' @export evaluateDecompositionQuality
evaluateDecompositionQuality <- function(exposure, signatures, genome,
                                         plot=FALSE) {
    if (!is.probability.object(genome)) {
        # must be a single genome
        stop(paste("'genome' must be an the mutation frequencies of an",
                   "individual genome (in Alexandrov or Shiraishi format)."))
    }
    
    if (!is.probability.vector(exposure)) {
        stop(paste("'exposure' must be the exposure vector of a single",
                   "decomposition."))
    }
    if (!isSignatureSet(signatures)) {
        stop("'signatures' must be a set (list) of signatures.")
    }

    # check the genome format; same as signature format?
    if (!sameSignatureFormat(list(genome), signatures)) {
        stop("Formats of genome and signatures must match!")
    }
    
    if (!is.logical(plot)) {
        stop("plot must be logical (TRUE or FALSE).")
    }
    

    decQual <- list()
    
    # compute explained variance
    decQual$explainedVariance <- computeExplainedVariance(exposure,
                                                          signatures, genome)

    # re-compose genome from exposure and signatures
    predGenome <- composeGenomesFromExposures(exposure, signatures)[[1]]

    # compute correlation of re-composed and original genome
    decQual$pearsonCorr <- cor(as.numeric(as.matrix(predGenome)),
                               as.numeric(as.matrix(genome)), method="pearson")
    # (as.numeric(as.matrix(X)) works on vectors, matrices and data frames!)

    
    if (!plot) {
        return (decQual)
    }
    # else continue and plot

    scatter.smooth(x=genome, y=predGenome,
                   main="Quality of decomposition",
                   #main="Observed and reconstructed mutation frequencies",
                   sub="mutation frequencies",
                   xlab="observed", ylab="reconstructed",
                   pch=20, cex=0.7)
    text(x=max(genome),y=0, pos=2,
         labels=paste("Explained variance:",
                      round(decQual$explainedVariance, digits=4)))
    text(x=max(genome),y=max(predGenome)/25, pos=2,
         labels=paste("Pearson correlation:",
                      round(decQual$pearsonCorr, digits=4)))
}

