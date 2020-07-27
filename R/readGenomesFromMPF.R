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
#' enforceUniqueTrDir=TRUE, 
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
#' @param enforceUniqueTrDir (Optional) Used only if \code{trDir} is
#' \code{TRUE}. If \code{enforceUniqueTrDir} is TRUE (default), then mutations
#' which map to a region with multiple overlapping genes with opposing
#' transcription directions will be excluded from the analysis. If \code{FALSE},
#' the transcript direction encountered first in the transcript database (see
#' \code{transcriptAnno}) is assigned to the mutation. The latter was the
#' behavior until version 1.3.5 of \code{decompTumor2Sig} and is also the
#' behavior of \code{pmsignature}. However, it is preferable to exclude
#' these mutations from the count (default) because from mutation data alone
#' it cannot be inferred which of the two genes has the higher transcriptional
#' activity which might potentially be linked to the occurrence of the mutation.
#' (If you are unsure, use the default setting; this option exists mostly for
#' backward compatibility with older versions.)
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
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
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
#'          trDir=TRUE, enforceUniqueTrDir=TRUE, refGenome=refGenome,
#'          transcriptAnno=transcriptAnno, verbose=FALSE)
#' 
#' @importFrom plyr aaply
#' @importFrom utils read.table
#' @importFrom data.table as.data.table
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene
#' TxDb.Hsapiens.UCSC.hg19.knownGene
#' @export readGenomesFromMPF
readGenomesFromMPF <- function(file, numBases=5, type="Shiraishi", trDir=TRUE,
    enforceUniqueTrDir=TRUE,
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

    # remove all spaces (some files contain them, e.g., in the POS field ...
    # (this is faster for a table than a data.table, so we do it first)
    mpf <- gsub(" ", "", mpf)
    
    # data.table objects can be modified much more efficiently
    mpf <- as.data.table(mpf)

    # get only SNVs (one REF base and one ALT base)
    # (indels can be specified as, e.g. deletion, "AG > A" or as "G > -")
    snvRows <- (nchar(mpf$REF) == 1) & (nchar(mpf$ALT) == 1) &
        (mpf$REF != "-") & (mpf$ALT != "-") &
        (mpf$REF != "N") & (mpf$ALT != "N")
    mpf <- mpf[snvRows]

    # we need to map mutations to samples, so that we can create
    # VCF-like, dummy genotype information

    sampleIDlist <- sort(unique(mpf$SAMPLE))


    ## add genotype information
    #mpf$FORMAT <- rep("GT", nrow(mpf))
    #
    ## construct genotype information for individual samples
    #for (sid in sampleIDlist) {
    #    gt <- rep(NA, nrow(mpf))
    #    gt[mpf$SAMPLE == sid] <- "1/1"
    #
    #    mpf[,sid] <- gt
    #}
    #
    ## now remove the sample ID (we have it in the genotype information)
    #mpf$SAMPLE = NULL
    #
    ## we have now (example):
    ##> mpf[1:10,]
    ##X1  CHROM  POS      REF ALT FORMAT PD3851a PD3890a PD3904a PD3905a 
    ##  1 "chr1" "809687" "G" "C" "GT"   "1/1"   NA      NA      NA      
    ##  2 "chr1" "819245" "G" "T" "GT"   "1/1"   NA      NA      NA 
    #
    #genomes <- buildGenomesFromMutationData(snvs=as.matrix(mpf),
    #                                        numBases=numBases,
    #                                        type=type, trDir=trDir,
    #                                        uniqueTrDir=enforceUniqueTrDir,
    #                                        refGenome=refGenome,
    #                                        transcriptAnno=transcriptAnno,
    #                                        verbose=verbose)

    maxN <- 20        # don't process more than 20 genomes at once,
                      # otherwise the memory consumption explodes
    
    if (verbose & (length(sampleIDlist) > maxN)) {
        cat(paste0("Too many tumors genomes; processing in batches ",
                   "(max ",maxN,") to limit memory usage!\n"))
    }

    genomes <- list()
        
    while(length(sampleIDlist)) {

        if (length(sampleIDlist) < maxN) {
            # remainder is smaller than batch size
            maxN <- length(sampleIDlist)
        }

        # get subset of mutations for this sample
        mpfS <- mpf[mpf$SAMPLE %in% sampleIDlist[seq(maxN)],]
        mpfS$FORMAT <- rep("GT", nrow(mpfS))

        for (sid in sampleIDlist[seq(maxN)]) {
            gt <- rep(NA, nrow(mpfS))
            gt[mpfS$SAMPLE == sid] <- "1/1"

            mpfS[,sid] <- gt
        }

        mpfS$SAMPLE = NULL

        genomes <-
            c(genomes, 
              buildGenomesFromMutationData(snvs=as.matrix(mpfS),
                                           numBases=numBases,
                                           type=type, trDir=trDir,
                                           uniqueTrDir=enforceUniqueTrDir,
                                           refGenome=refGenome,
                                           transcriptAnno=transcriptAnno,
                                           verbose=verbose)
              )

        # remove this batch of samples from the list
        sampleIDlist <- sampleIDlist[-seq(maxN)]

    } # end while

    if(verbose && length(genomes) > 0) {
        cat("Done reading genomes.\n")
    }

    return(genomes)

}


