#' Read tumor genomes from an MPF file (Mutation Position Format).
#'
#' `readGenomesFromMPF()` reads somatic mutations of a single tumor genome
#' (sample) or a set of genomes from an MPF file (Mutation Position Format;
#' see details below) and determines the mutation frequencies according to
#' a specific model of mutational signatures (Alexandrov or Shiraishi).
#'
#' An MPF file has the following format (one line per mutation and
#' patient/sample):
#'
#' [sampleID]<tab>[chrom]<tab>[position]<tab>[ref_bases]<tab>[alt_bases]
#'
#' @usage readGenomesFromMPF(file, numBases=5, type="Shiraishi", trDir=TRUE,
#' refGenome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
#' transcriptAnno=
#' TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
#' verbose=TRUE)
#' @param file (Mandatory) The name of the MPF file (can be compressed with
#' \code{gzip}).
#' @param numBases (Mandatory) Total number of bases (mutated base and
#' flanking bases) to be used for sequence patterns. Must be odd. Default: 5
#' @param type (Mandatory) Signature model or type (\code{"Alexandrov"} or
#' \code{"Shiraishi"}). Default: \code{"Shiraishi"}
#' @param trDir (Mandatory) Specifies whether the transcription direction is
#' taken into account in the signature model. If so, only mutations within
#' genomic regions with a defined transcription direction can be considered.
#' Default: \code{TRUE}
#' @param refGenome (Mandatory) The reference genome (\code{BSgenome}) needed
#' to extract sequence patterns. Default: \code{BSgenome} object for hg19.
#' @param transcriptAnno (Optional) Transcript annotation (\code{TxDb} object)
#' used to determine the transcription direction. This is required only if
#' \code{trDir} is \code{TRUE}. Default: \code{TxDb} object for hg19.
#' @param verbose (Optional) Print information about reading and processing the
#' mutation data. Default: \code{TRUE}
#' @return A list containing the genomes in terms of frequencies of the mutated
#' sequence patterns. This list of genomes can be used for
#' \code{decomposeTumorGenomes}. 
#' @author Rosario M. Piro\cr Freie Universitaet Berlin\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <r.piro@@fu-berlin.de>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2017) Identification of Mutational Signatures Active in
#' Individual Tumors. NETTAB 2017 - Methods, Tools & Platforms for
#' Personalized Medicine in the Big Data Era, October 16-18, 2017, Palermo,
#' Italy. PeerJ Preprints 5:e3257v1, 2017.
#' @seealso \code{\link{decompTumor2Sig}}\cr
#' \code{\link{decomposeTumorGenomes}}\cr
#' \code{\link{readGenomesFromVCF}}\cr
#' \code{\link{getGenomesFromMutFeatData}}
#' @examples
#' 
#' ### load reference genome and transcript annotation (if direction is needed)
#' refGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' transcriptAnno <-
#'   TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#' 
#' ### read breast cancer genomes from Nik-Zainal et al (PMID: 22608084) 
#' gfile <- system.file("extdata", "Nik-Zainal_PMID_22608084-MPF.txt.gz", 
#'          package="decompTumor2Sig")
#' genomes <- readGenomesFromMPF(gfile, numBases=5, type="Shiraishi",
#'          trDir=TRUE, refGenome=refGenome, transcriptAnno=transcriptAnno,
#'          verbose=FALSE)
#' 
#' @importFrom plyr aaply
#' @importFrom utils read.table
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene
#' TxDb.Hsapiens.UCSC.hg19.knownGene
#' @export readGenomesFromMPF
readGenomesFromMPF <- function(file, numBases=5, type="Shiraishi", trDir=TRUE,
    refGenome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
    transcriptAnno=
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
    verbose=TRUE) {

    # read mutation data
    if(verbose) {
        cat("Reading mutations and patient/sample IDs from MPF file:\n")
    }
    mpf <- as.matrix(read.table(file, header=FALSE, row.names=NULL, sep="\t"))
    colnames(mpf) <- c("SAMPLE", "CHROM", "POS", "REF", "ALT")

    smpIndex <- which(colnames(mpf) == "SAMPLE")
    chrIndex <- which(colnames(mpf) == "CHROM")
    posIndex <- which(colnames(mpf) == "POS")
    refIndex <- which(colnames(mpf) == "REF")
    altIndex <- which(colnames(mpf) == "ALT")

    # remove all spaces (some files contain them, e.g., in the POS field ...
    mpf <- gsub(" ", "", mpf)
    
    # get only SNVs (one REF base and one ALT base)
    snvRows <- (nchar(mpf[,refIndex]) == 1) & (nchar(mpf[,altIndex]) == 1)
    mpf <- mpf[snvRows,]

    
    # we need to map mutations to samples, so that we can create
    # VCF-like, dummy genotype information

    if(verbose) {
        cat(paste("Collapsing variant information by mapping multiple",
                  "samples to unique variants.\n"))
    }
    
    sampleIDlist <- sort(unique(mpf[,smpIndex]))

    snvs <- plyr::aaply(mpf, 1, function(x) {
        gtVec <- c("GT", rep(NA, length(sampleIDlist)))
        names(gtVec) <- c("FORMAT", sampleIDlist)
        gtVec[x[smpIndex]] <- "1/1"
        
        c(x[c(chrIndex,posIndex,refIndex,altIndex)], gtVec)
    })

    
    ## we have now (example):
    #> snvs[1:10,]
    #X1  CHROM  POS      REF ALT FORMAT PD3851a PD3890a PD3904a PD3905a PD3945a
    #  1 "chr1" "809687" "G" "C" "GT"   "1/1"   NA      NA      NA      NA     
    #  2 "chr1" "819245" "G" "T" "GT"   "1/1"   NA      NA      NA      NA     

    genomes <- buildGenomesFromMutationData(snvs=snvs, numBases=numBases,
                                            type=type, trDir=trDir,
                                            refGenome=refGenome,
                                            transcriptAnno=transcriptAnno,
                                            verbose=verbose)

    if(verbose) {
        cat("Done reading genomes.\n")
    }

    return(genomes)

}


