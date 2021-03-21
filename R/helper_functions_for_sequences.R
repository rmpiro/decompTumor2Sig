

###########################################
# internal helper functions for sequences #
###########################################

#' Compute nucleotide frequencies (internal function)
#'
#' Compute nucleotide frequencies for a reference genome and a region set.
#'
#' @usage compNucFreq(refGenome, regions=NULL, numBases=1, mergeByRevComp=TRUE)
#' @param refGenome Reference genome (\code{BSgenome})
#' @param regions (Optional) Regions (\code{GRanges}); default:
#' \code{NULL} (whole genome).
#' @param numBases (Optional) Sequence pattern length for which to compute
#' frequencies (e.g., '3' for trinucleotides). Default: 1 (single nucleotides).
#' @param mergeByRevComp (Optional) Reduce redundancy by merging counts
#' for reverse complement sequence patterns. This function allows this only
#' for sequence patterns with odd length (e.g., single nucleotides,
#' trinucleotides). Default: \code{TRUE}.
#' @return Vector of computed nucleotide frequencies
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @importFrom Biostrings oligonucleotideFrequency DNAStringSet
#' reverseComplement
#' @keywords internal
compNucFreq <- function(refGenome, regions=NULL, numBases=1,
                        mergeByRevComp=TRUE) { 

    if (is.null(regions)) {
        sequences <- Biostrings::getSeq(refGenome)
    } else {
        sequences <- Biostrings::getSeq(refGenome, regions)
    }

    # compute nucleotide pattern counts for individual sequences
    counts <- oligonucleotideFrequency(sequences, numBases)
    # sum counts over all sequences
    if(is.matrix(counts)) {
        counts <- apply(counts, 2, sum)
    }

    # merge counts from reverse complement patterns?
    if (mergeByRevComp) {
        if ((numBases %% 2) == 0) {
            stop(paste0("Nucleotide patterns must be of odd length for ",
                        "merging reverse complements!"))
        }
    
        # nucloetide patterns are equivalent to their reverseComplement!
        revNuc <- as.character(reverseComplement(DNAStringSet(names(counts))))
        counts <- counts + counts[revNuc]

        # take only those with C and T at center (convention of signatures)
        pattern <- paste0(strrep(".", (numBases-1)/2),
                          "[CT]", strrep(".", (numBases-1)/2))
    
        counts <- counts[grep(pattern, names(counts))]
    }

    # convert to frequencies
    freqs <- counts / sum(counts)

    return(freqs)
}


#' Convert sequence pattern frequencies to base frequencies (internal function)
#'
#' Convert sequence pattern frequencies (Alexandrov-like) to frequencies for
#' individual bases (Shiraishi-like).
#'
#' @usage convertSeqFreqToBaseFreq(frequencies, sigType)
#' @param frequencies Vector for frequencies of individual sequence patterns
#' @param sigType Signature type information as returned by
#' \code{\link{determineTypeNumBasesAndTrDir}}
#' @return Shiraishi signature-like matrix of per-base nucleotide frequencies
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @keywords internal
convertSeqFreqToBaseFreq <- function(frequencies, sigType) { 

    numBases <- sigType$numBases    # number of bases in sequence pattern
    center <- (numBases+1)/2        # position of central base
    
    # initialize a Shiraishi signature-like matrix for frequencies
    # add a row for transcription direction info, if necessary;
    # initialize everything to 1 because we later divide one matrix by another;
    # hence the initial base factor for normalization is 1 for each element
    # in the matrix
    freqMatrix <- matrix(1, ncol=6, nrow=(numBases+sigType$trDir))

    # compute nucleotide frequencies for individual bases of the pattern
    for (pos in seq(numBases)) {
        # for this position, get the corresponding base for each element in
        # the frequency vector
        posBases <- substr(names(frequencies), start=pos, stop=pos)

        # split the frequency vector according to the bases
        freqSplit <- split(frequencies, posBases)

        # convert split list to overall base frequencies (by summing)
        baseFreq <- vapply(freqSplit, sum, FUN.VALUE=numeric(1))

        # decide where to put this in the matrix
        if (pos == center) {
            # central (mutated) base; row 1 of the matrix
            # each frequency 3 times (row is: C>A, C>G, C>T, T>A, T>C, T>G)
            freqMatrix[1,] <- rep(baseFreq, each=3)
        } else {
            # flanking bases
            # if (pos<center): upstream base; starting from row 2 in the matrix
            #                  row=pos+1 (because 1st row is the center!)
            # if (pos>center): downstream base; row=pos
            # just copy to first 4 row elements (row is: A, C, G, T, 0, 0)
            freqMatrix[pos+(pos<center),] <- c(baseFreq, rep(1,2))
            # Note: we use dummy score 1 for two remaining row elements
            # because we later divide one matrix by another, so that the
            # respective factors becomes 1
        }
    }

    return(freqMatrix)
}


#' Adjust an Alexandrov signature (internal function)
#'
#' Adjust an Alexandrov signature by a set of pattern-specific factors.
#'
#' @usage adjustAlexandrovSignature(signature, factors)
#' @param signature Single Alexandrov signature to be adjusted
#' @param factors Factors for trinucleotides 
#' @return Adjusted Alexandrov signature
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signature active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @keywords internal
adjustAlexandrovSignature <- function(signature, factors) { 

    # keep names of the mutation types
    muttypes <- names(signature)

    # get trinucleotides in the same order
    tri <- gsub("\\[([CT])>.\\]", "\\1", muttypes)

    # adjust mutation probabilities according to trinucleotide factors
    sign <- signature * factors[tri]

    # finally, re-normalize signature to sum up to 1
    sign <- sign / sum(sign)

    return(sign)
}



#' Adjust an Shiraishi signature (internal function)
#'
#' Adjust an Shiraishi signature by a set of base-specific factors.
#'
#' @usage adjustShiraishiSignature(signature, factorsMatrix)
#' @param signature Single Shiraishi signature to be adjusted
#' @param factorsMatrix Factors for single nucleotides (Shiraishi-like matrix) 
#' @return Adjusted Shiraishi signature
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signature active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @keywords internal
adjustShiraishiSignature <- function(signature, factorsMatrix) { 

    if ( !(is.matrix(signature) && is.matrix(factorsMatrix)) ||
         !(all(dim(signature) == dim(factorsMatrix))) ) {
        stop("signature and factorsMatrix need to be matrices of same format!")
    }

    # multiply signature elementwise by adjustment factors
    sign <- signature * factorsMatrix

    # renormalize each row (=base of the pattern) to sum up to 1
    sign <- t(apply(sign, 1, function(x){ x/sum(x)}))
    
    return(sign)
}


