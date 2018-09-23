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
#' trDir=TRUE,
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
#'          trDir=TRUE, refGenome=refGenome, transcriptAnno=transcriptAnno,
#'          verbose=FALSE)
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
    numBases=5, type="Shiraishi", trDir=TRUE,
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
                                            refGenome=refGenome,
                                            transcriptAnno=transcriptAnno,
                                            verbose=verbose)

    if(verbose) {
        cat("Done converting genomes.\n")
    }

    return(genomes)

}

