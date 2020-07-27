#' Convert genomes from a \code{VRanges} object
#'
#' `convertGenomesFromVRanges()` converts the SNVs of a single tumor genome
#' (sample) or a set of genomes from a \code{VRanges} object (package
#' \code{VariantAnnotation}) and determines the mutation frequencies according
#' to a specific model of mutational signatures (Alexandrov or Shiraishi),
#' such that the resulting format can be used as genomes input for
#' \code{decomposeTumorGenomes}.
#'
#' @usage convertGenomesFromVRanges(vranges, numBases=5, type="Shiraishi",
#' trDir=TRUE, enforceUniqueTrDir=TRUE, 
#' refGenome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
#' transcriptAnno=
#' TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
#' verbose=TRUE)
#' @param vranges (Mandatory) The \code{VRanges} object which specifies the
#' mutations.
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
#' @param verbose (Optional) Print information about reading and processing
#' the mutation data. Default: \code{TRUE}
#' @return A list containing the genomes in terms of frequencies of the
#' mutated sequence patterns. This list of genomes can be used for
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
#' \code{\link{readGenomesFromMPF}}\cr
#' \code{\link{getGenomesFromMutFeatData}}
#' @examples
#' 
#' ### load the reference genome and the transcript annotation database
#' refGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' transcriptAnno <-
#'   TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#' 
#' ### take the breast cancer genomes from Nik-Zainal et al (PMID: 22608084) 
#' gfile <- system.file("extdata",
#'          "Nik-Zainal_PMID_22608084-VCF-convertedfromMPF.vcf.gz", 
#'          package="decompTumor2Sig")
#' 
#' ### get the corresponding VRanges object (using the VariantAnnotation
#' ### package)
#' library(VariantAnnotation)
#' vr <- readVcfAsVRanges(gfile, genome="hg19")
#' 
#' ### convert the VRanges object to the decompTumor2Sig format
#' genomes <- convertGenomesFromVRanges(vr, numBases=5, type="Shiraishi",
#'          trDir=TRUE, enforceUniqueTrDir=TRUE, refGenome=refGenome, 
#'          transcriptAnno=transcriptAnno, verbose=FALSE)
#' 
#' @importFrom VariantAnnotation asVCF isSNV readVcfAsVRanges header geno ref
#' alt
#' @importFrom GenomicRanges ranges seqnames start
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom methods is
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene
#' TxDb.Hsapiens.UCSC.hg19.knownGene
#' @export convertGenomesFromVRanges
convertGenomesFromVRanges <- function(vranges,
    numBases=5, type="Shiraishi", trDir=TRUE, enforceUniqueTrDir=TRUE,
    refGenome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
    transcriptAnno=
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
    verbose=TRUE) {

    if (!is(vranges, "VRanges")) {
        stop("vranges must be an object of type VRanges!")
    }
    
    # read mutation data
    if(verbose) {
        cat(paste("Extracting mutations and sample/genotype information",
                  "from VRanges:\n"))
    }

    # convert VRanges to VCF
    vcf <- asVCF(vranges)

    # reduce to contain only SNVs
    vcf <- vcf[isSNV(vcf, singleAltOnly=FALSE)]
    
    samples <- samples(header(vcf))


    # construct table with SNVs
    sampleGT <- geno(vcf)$GT
    sampleGT[sampleGT == "."] <- NA

    snvs <- cbind(as.character(seqnames(rowRanges(vcf))),  # CHROM
                  start(ranges(rowRanges(vcf))),           # POS
                  as.character(ref(vcf)),                  # REF
                  as.character(alt(vcf)),                  # ALT
                  rep("GT", nrow(sampleGT)))               # FORMAT
    colnames(snvs) <- c("CHROM", "POS", "REF", "ALT", "FORMAT")

    snvs <- cbind(snvs, sampleGT)
    
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
        cat("Done converting genomes.\n")
    }

    return(genomes)

}

