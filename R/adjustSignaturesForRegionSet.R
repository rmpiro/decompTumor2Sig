#' Adjust (normalize) signatures for a set of genomic regions.
#'
#' `adjustSignaturesForRegionSet()` takes a set of signatures that have
#' been orginally defined with respect to the nucleotide frequencies within
#' a specific reference genome or region (e.g., by deriving them from whole
#' genome mutation data) and adjusts or normalizes them to the often different
#' nucleotide frequencies of another specific subset of genomic regions.
#'
#' This may be useful, for example, to perform signature refitting (using
#' \code{\link{decomposeTumorGenomes}}) for mutation data from targetted
#' sequencing (e.g., only a subset of genes), whole exome sequencing (only
#' exonic regions), or other limited subsets of the genome with particular
#' nucleotide frequencies.
#'
#' For Alexandrov-type signatures, the important frequencies are those of the
#' whole sequence patterns (e.g., trinucleotides) whose central base can be
#' mutated. Therefore, adjustment factors for individual mutation types
#' (e.g., A[C>T]G) are computed by comparing the corresponding sequence
#' pattern frequencies (e.g., ACG) between the original reference
#' regions (e.g., whole genome) and the target regions (e.g., target
#' regions of whole exome sequencing).
#' 
#' In the Shiraishi-type signature model, individual bases of the sequence
#' patterns are considered as independent features. Thus, to compute nucleotide
#' frequencies for such signatures, the frequencies of the sequence patterns
#' (computed as for Alexandrov-type signatures) are broken down to
#' single-nucleotide frequencies for the individual positions of the patterns.
#'
#' In both cases, after the appropriate adjustment of individual features,
#' signatures are re-normalized such that overall probabilites sum up to 1.
#'
#' @usage adjustSignaturesForRegionSet(signatures,
#' regionsTarget, regionsOriginal=NULL, 
#' refGenome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
#' @param signatures (Mandatory) Signatures to be adjusted to the nucleotide
#' frequencies of the genomic regions defined by the parameter \code{regions}.
#' @param regionsTarget (Mandatory) \code{GRanges} object defining a subset
#' of the genome (i.e., a set of genomic regions) for which the signatures
#' need to be adjusted (can be set to \code{NULL} for the whole genome).
#' @param regionsOriginal (Optional) \code{GRanges} object defining the subset
#' of the genome (i.e., set of genomic regions) from which the signatures
#' where originally derived. Default: \code{NULL} (whole genome).
#' @param refGenome (Optional) Reference genome sequence from which to
#' compute the nucleotide frequencies. Default:\cr
#' \code{BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19}.
#' @return A set of adjusted mutational signatures in the same format as those
#' specified for \code{signatures}.
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @seealso \code{\link{decompTumor2Sig}}\cr
#' \code{\link{decomposeTumorGenomes}}
#' @examples
#' 
#' ### get Alexandrov signatures from COSMIC
#' signatures <- readAlexandrovSignatures()
#' 
#' ### get gene annotation for the default reference genome (hg19)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' ### get a GRanges object for gene promoters (-2000 to +200 bases from TSS)
#' ### [taking only the first 1000 for testing purpose]
#' library(GenomicRanges)
#' regionsTarget <- promoters(txdb, upstream=2000, downstream=200)[seq(1000)]
#'
#' ### assume these signatures were derived only from mutation data from
#' ### exons on chromosome X [not true; just for illustrative purpose]
#' filter <- list(tx_chrom = c("chrX"))
#' regionsOriginal <- exons(txdb, filter=filter)
#' 
#' ### adjust signatures according to nucleotide frequencies in the target
#' ### subset of the genome
#' sign_adj <- adjustSignaturesForRegionSet(signatures, regionsTarget,
#'                                                      regionsOriginal)
#'
#' @importFrom GenomicRanges start end trim
#' @importFrom Biostrings getSeq
#' @export adjustSignaturesForRegionSet
adjustSignaturesForRegionSet <- function(signatures,
                                         regionsTarget, regionsOriginal=NULL,
        refGenome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19) {
    
    if (!isSignatureSet(signatures)) {
        stop("Parameter 'signatures' must be a set (list) of signatures!")
    }

    if (!is.null(regionsTarget) && !is(regionsTarget, "GRanges")) {
        stop("Parameter 'regionsTarget' must be NULL or of type GRanges!")
    }

    if (!is.null(regionsOriginal) && !is(regionsOriginal, "GRanges")) {
        stop("Parameter 'regionsOriginal' must be NULL or of type GRanges!")
    }

    if ( (is.null(regionsTarget) && is.null(regionsOriginal) ) ||
         (!is.null(regionsTarget) && !is.null(regionsOriginal) &&
          length(regionsTarget) == length(regionsOriginal) &&
          all(regionsOriginal == regionsTarget)
         )
        ) {
        stop("Parameters 'regionsOriginal' and 'regionsTarget' must differ!")
    }

    # Size (and type) of sequence patterns of the signatures
    sigType <- determineTypeNumBasesAndTrDir(signatures[[1]])
    numBases <- sigType$numBases

    # If necessary, extend the regions
    if (numBases > 1) {
        # With sequence patterns of trinculeotides (or more) each nucleotide in
        # a region set needs to be the _center_ of a sequence pattern, 
        # so we need to extend each region by the number of flanking bases!

        if (!is.null(regionsTarget)) {
            tryCatch( {
                GenomicRanges::start(regionsTarget) <-
                    GenomicRanges::start(regionsTarget) - (numBases-1)/2
                GenomicRanges::end(regionsTarget) <-
                    GenomicRanges::end(regionsTarget) + (numBases-1)/2
                # Got a warning here? we're out of bounds and need to trim
            }, warning=function(w) {regionsTarget <-
                                        GenomicRanges::trim(regionsTarget)} )
            
        }

        if (!is.null(regionsOriginal)) {
            tryCatch( {
                GenomicRanges::start(regionsOriginal) <-
                    GenomicRanges::start(regionsOriginal) - (numBases-1)/2
                GenomicRanges::end(regionsOriginal) <-
                    GenomicRanges::end(regionsOriginal) + (numBases-1)/2
                # Got a warning here? we're out of bounds and need to trim
            }, warning=function(w) {regionsOriginal <-
                                        GenomicRanges::trim(regionsOriginal)} )
        }
    }


    # Compute nucleotide frequencies for original reference sequences
    # (merging counts from reverse complement sequence patterns)
    freqO <- compNucFreq(refGenome, regionsOriginal,
                         numBases=numBases, mergeByRevComp=TRUE)

    # Compute the same for target sequences
    freqT <- compNucFreq(refGenome, regionsTarget,
                         numBases=numBases, mergeByRevComp=TRUE)


    # compute adjustment factors: 
    signaturesAdj <- NULL

    if (isAlexandrovSet(signatures)) {
        # These are Alexandrov signatures with a full dependency between all
        # bases of the basic sequence pattern (usually trinucleotides)
        
        # Divide by original (genome-wide?) frequency, then multiply 
        # by frequency in the target region set, i.e., if the limited
        # region has double frequency of a trinucleotide (e.g., ACG) with 
        # respect to the whole genome, the mutation probabilities of the 
        # corresponding mutation types (A[C>A]G, A[C>G]G, A[C>T]G) in the 
        # signatures will be multiplied by 2. Reasoning: the signatures
        # has more opportunities to mutate that trinucleotide in the limited
        # region set then in the whole genome.
        factors <- freqT / freqO

        # now adjust all signatures according to these factors
        signaturesAdj <- base::lapply(signatures,
                                      adjustAlexandrovSignature, factors)

    } else if (isShiraishiSet(signatures)) {
        # these are Shiraishi signatures; need to compute frequencies for
        # the individual bases of the sequence pattern (e.g., frequencies of
        # A,C,G,T for the first base, etc., and frequencies of C and T for
        # the central base; these are determined from the frequencies of
        # complete sequence patterns (freqOMerged, freqTMerged)

        # convert pattern frequency vectors to base frequency matrices 
        freqOMatrix <- convertSeqFreqToBaseFreq(freqO, sigType)
        freqTMatrix <- convertSeqFreqToBaseFreq(freqT, sigType)
        
                
        # adjustment factors 
        factorsMatrix <- freqTMatrix / freqOMatrix

        # now adjust all signatures according to these factors
        signaturesAdj <- base::lapply(signatures,
                                      adjustShiraishiSignature, factorsMatrix)
    } else {
        stop("Unknown signature type!")
    }

    return(signaturesAdj)
}

