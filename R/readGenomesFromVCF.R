#' Read tumor genomes from a VCF file (Variant Call Format).
#'
#' `readGenomesFromVCF()` reads somatic mutations of a single tumor genome
#' (sample) or a set of genomes from a VCF file (Variant Call Format) and
#' determines the mutation frequencies according to a specific model of
#' mutational signatures (Alexandrov or Shiraishi).
#'
#' @usage readGenomesFromVCF(file, numBases=5, type="Shiraishi", trDir=TRUE,
#' enforceUniqueTrDir=TRUE,
#' refGenome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
#' transcriptAnno=
#' TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
#' verbose=TRUE)
#' @param file (Mandatory) The name of the VCF file (can be compressed with
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
#' \code{\link{decomposeTumorGenomes}}\cr
#' \code{\link{readGenomesFromMPF}}\cr
#' \code{\link{getGenomesFromMutFeatData}}
#' @examples
#' 
#' ### load reference genome and transcript annotation (if direction is needed)
#' refGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' transcriptAnno <-
#'   TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#' 
#' ### read breast cancer genomes from Nik-Zainal et al (PMID: 22608084) 
#' gfile <- system.file("extdata",
#'          "Nik-Zainal_PMID_22608084-VCF-convertedfromMPF.vcf.gz", 
#'          package="decompTumor2Sig")
#' genomes <- readGenomesFromVCF(gfile, numBases=5, type="Shiraishi",
#'          trDir=TRUE, enforceUniqueTrDir=TRUE, refGenome=refGenome,
#'          transcriptAnno=transcriptAnno, verbose=FALSE)
#' 
#' @importFrom vcfR read.vcfR extract.gt getFIX
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene
#' TxDb.Hsapiens.UCSC.hg19.knownGene
#' @export readGenomesFromVCF
readGenomesFromVCF <- function(file, numBases=5, type="Shiraishi", trDir=TRUE,
    enforceUniqueTrDir=TRUE,
    refGenome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
    transcriptAnno=
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
    verbose=TRUE) {

    # read mutation data
    if(verbose) {
        cat("Reading mutations and genotype information from VCF file:\n")
    }
    vcf <- read.vcfR(file, verbose=verbose, checkFile=TRUE, convertNA=TRUE)

    # get only SNVs (one REF base and one ALT base)
    # (indels can be specified as, e.g. deletion, "AG > A" or as "G > -")
    snvRows <-
        (nchar(getFIX(vcf)[,"REF"]) == 1) & (nchar(getFIX(vcf)[,"ALT"]) == 1) &
            (getFIX(vcf)[,"REF"] != "-") & (getFIX(vcf)[,"ALT"] != "-")

    # basic variant information (chr, pos, ref, alt)
    snvs <- getFIX(vcf)[snvRows, c("CHROM","POS","REF","ALT")]

    # get genotypes (if FORMAT and genotypes are present in the VCF)
    gtypes <- tryCatch(extract.gt(vcf)[snvRows,],
                       error=function(e) { return(NULL) })

    # if this was only one sample: have a vector, need keep sample name
    gtypes <- as.matrix(gtypes)
    colnames(gtypes) <- colnames(extract.gt(vcf))
 
    # genotype information
    if (!is.null(gtypes)) {
        # if we do have genotype information, add it!
        snvs <- cbind(snvs,                        # CHR, POS, REF, ALT
                      "FORMAT"=rep("GT", sum(snvRows)),  # FORMAT
                      gtypes)                      # genotype info for samples
        
    } else {
        # we don't have genotype information from the VCF file
        # take all variants as originating from the same sample/genome,
        # create dummy genotype
        snvs <- cbind(snvs, matrix(c("GT", "1/1"),
                                   nrow=length(which(snvRows)),
                                   ncol=2, byrow=TRUE)
                      )
        colnames(snvs)[seq((ncol(snvs)-1),ncol(snvs))] =
            c("FORMAT", "variants_without_genotype_info")
    }

    ## we have now (example):
    #> snvs[1:10,]
    #     CHROM POS  REF ALT FORMAT           sample1                 sample2
    #[1,] "2"  "947" "C" "T" "GT:PL:GQ:AD:DP" "1/1:84,6,0:6:0,2:2"    NA    
    #[2,] "2"  "992" "G" "A" "GT:PL:GQ:AD:DP" "0/1:123,0,33:33:1,3:4" "0/0:..."

    genomes <- buildGenomesFromMutationData(snvs=snvs, numBases=numBases,
                                            type=type, trDir=trDir,
                                            uniqueTrDir=enforceUniqueTrDir,
                                            refGenome=refGenome,
                                            transcriptAnno=transcriptAnno,
                                            verbose=verbose)

    if(verbose) {
        cat("Done reading genomes.\n")
    }

    return(genomes)

}
