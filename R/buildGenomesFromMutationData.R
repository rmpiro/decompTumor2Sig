

#####################
# internal function #
#####################

#' buildGenomesFromMutationData (internal function)
#'
#' Build genome data structures (same as signatures) and fill then with
#' mutation data.
#'
#' SNVs are specified as a matrix of the following format (adapted from VCF):\cr
#' #> snvs[1:2,]\cr
#' #     CHROM POS  REF ALT FORMAT           sample1                 sample2\cr
#' #[1,] "2"  "947" "C" "T" "GT:PL:GQ:AD:DP" "1/1:84,6,0:6:0,2:2"    NA\cr
#' #[2,] "2"  "992" "G" "A" "GT:PL:GQ:AD:DP" "0/1:123,0,33:33:1,3:4" "0/0:..."
#' 
#' @usage buildGenomesFromMutationData(snvs, numBases, type, trDir,
#' uniqueTrDir=TRUE, refGenome, transcriptAnno, verbose)
#' @param snvs SNV matrix (see description above).
#' @param numBases Number of bases for the sequence pattern (odd integer).
#' @param type Type of signature to be used ("Alexandrov", "Shiraishi").
#' @param trDir Logical: use transcription-strand information?
#' @param uniqueTrDir Logical; used only if trDir is also \code{TRUE}: if
#' \code{uniqueTrDir} is \code{TRUE} (default), then only mutations with only
#' one defined transcription strand will be used, mutations for which both 
#' strands are valid are ignored. If \code{FALSE}, these mutations are accepted 
#' and one of the two transcription strands will be arbitrarily taken (the  
#' first one encountered in the databse specified for \code{transcriptAnno}). 
#' The latter was the behavior until version 1.3.5 of \code{decompTumor2Sig} 
#' and is also the behavior of \code{pmsignature}.
#' @param refGenome Reference genome (\code{BSgenome} object).
#' @param transcriptAnno Transcription information (\code{TxDb} object).
#' @param verbose Logical. Print additional information?
#' @return A list of genomes: each genome is represented by the observed
#' frequencies of mutation patterns according to the selected signature type.
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
#' @importFrom GenomicRanges seqnames makeGRangesFromDataFrame resize
#' findOverlaps
#' @importFrom Biostrings getSeq reverseComplement
#' @importFrom GenomicFeatures transcripts
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom BiocGenerics strand
#' @importFrom plyr alply
#' @keywords internal
buildGenomesFromMutationData <- function(snvs, numBases, type, trDir,
                                         uniqueTrDir = TRUE,
                                         refGenome, transcriptAnno, verbose) {

    if (type == "independent") {
        type <- "Shiraishi"
    } else if (type == "full") {
        type <- "Alexandrov"
    } else if (type != "Shiraishi" && type != "Alexandrov") {
        stop("Unkown signature 'type'!")
    }

    if (!is.logical(trDir)) {
        stop("The parameter 'trDir' must be logical/boolean!")
    }

    if (numBases < 1 || (numBases %% 2) != 1) {
        stop("The parameter 'numBases' must be positive and odd!")
    }

    
    # check for MT versus M for mitochondrion
    # first: ref. genome
    if(verbose) {
        cat(paste("Comparing specification of mitochondrion in mutation data",
                  "(M or MT?) to reference genome and adjusting if",
                  "necessary.\n"))
    }

    # check whether the mitochondrion is specified as chromosome "MT" or "M"

    refMT <- ""
    if (length(grep("MT", GenomicRanges::seqnames(refGenome))) == 1) { 
        refMT <- "MT"
    } else if (length(grep("M", GenomicRanges::seqnames(refGenome)))
               == 1) {
        refMT <- "M"
    } # else empty, no clear mitochondrion found

    #second: mutation file
    fileMT <- ""
    if (length(grep("MT", unique(snvs[,"CHROM"]))) == 1) { # having MT
        fileMT <- "MT"
    } else if (length(grep("M", unique(snvs[,"CHROM"]))) == 1) {
        fileMT <- "M"
    } # else empty, no clear mitochondrion found

    # third: check and adjust mutation data, if necessary
    if (fileMT != "" && fileMT != refMT) {
        if (refMT != "") {
            # invert
            snvs[,"CHROM"] <- gsub(fileMT, refMT, snvs[,"CHROM"])
        }
    }


    
    # check if chromosome IDs match what we have in the reference genome
    # (with or without prefix "chr")
    if(verbose) {
        cat(paste("Checking chromosome IDs of mutation data (1 or chr1?)",
                  "and adjusting if necessary.\n"))
    }

    setF = unique(snvs[,"CHROM"])
    setR = GenomicRanges::seqnames(refGenome)

    # check what "chr"-prefix we have in the ref. genome: chr? Chr? CHR? none?
    refPrefix <- unique(substr(grep("^chr", setR,
                                    ignore.case=TRUE, value=TRUE), 1, 3))
    if (!length(refPrefix)) {
        refPrefix <- ""
    }
    
    # check what "chr"-prefix we have in the mut. file: chr? Chr? CHR? none?
    filePrefix <- unique(substr(grep("^chr", setF,
                                    ignore.case=TRUE, value=TRUE), 1, 3))
    if (!length(filePrefix)) {
        filePrefix <- ""
    }

    if (refPrefix != filePrefix) {
        # not the same "chr"-prefix
        if (verbose) {
            cat(paste0("Having prefix for chromosome names: '", filePrefix,
                       "'; requiring prefix for reference genome: '",
                       refPrefix,"' ... replacing.\n"))
        }
        
        # first, remove whatever different prefix we have in the mutation file
        if (filePrefix != "") {
            snvs[,"CHROM"] <- gsub(paste0("^",filePrefix), "", snvs[,"CHROM"])
        }

        # then, add whatever prefix we need for the reference genome
        if (refPrefix != "") {
            snvs[,"CHROM"] <- paste(refPrefix,snvs[,"CHROM"], sep="")
        }
    }

    # now check, which chromosomes we don't find in the reference
    # genome and exclude them
    setF = unique(snvs[,"CHROM"])
    # setR has in any case remained the same

    excludeSeq <- setF[!(setF %in% setR)]

    if (length(excludeSeq)) {
        cat(paste("Warning: cannot find the following chromosomes in the",
                  "reference genome (will be ignored):\n"))
        cat(excludeSeq)
        cat("\n")

        snvs <- snvs[!(snvs[,"CHROM"] %in% excludeSeq),]
    }

    if (verbose) {
        cat(paste("Processing a total of", nrow(snvs), "mutations.\n"))
    }

    
    # get sequences from the reference genome
    if (verbose) {
        cat(paste("Extracting sequences with flanking bases of mutations",
                  "from the reference genome.\n"))
    }
    
    gr <- GenomicRanges::makeGRangesFromDataFrame(
                      data.frame(chr = snvs[,"CHROM"], 
                                 start = as.numeric(snvs[,"POS"]), 
                                 end = as.numeric(snvs[,"POS"])),
                      ignore.strand = TRUE)
    seqs <- Biostrings::getSeq(refGenome, GenomicRanges::resize(gr, numBases,
                                                                fix = "center"))


    # verify that the obtained center base is the REF base specified in the
    # mutation data (i.e., check that this is the correct reference genome 
    # for the mutation data)

    if( ! all( snvs[,"REF"] ==
                  substr(as.character(seqs), (numBases%/%2)+1, (numBases%/%2)+1)
              )
       ) {
        stop(paste("Reference (REF) bases in mutation data do not match the",
                   "specified reference genome!"))
    }

    
    # set strand info to NA
    strands <- rep(NA, length(seqs))
    
    # get strand info if we need it
    if(trDir) {
        if(verbose) {
            cat(paste("Getting information on transcription directions for",
                      "mutations within genes.\n"))
        }

        tr <- GenomicFeatures::transcripts(transcriptAnno)
        trHit <- GenomicRanges::findOverlaps(gr, tr, ignore.strand = FALSE)
        trStr <- cbind(S4Vectors::queryHits(trHit),
                       as.character(BiocGenerics::strand(
                           tr[S4Vectors::subjectHits(trHit)])
                                    )
                       )
        trStrUnique <- unique(trStr[trStr[,2] != "*",],MARGIN=1)

        if (!uniqueTrDir) {
            ## the following was the approach also used by pmsignature:
            ## for mutations in regions with transcripts in both directions,
            ## use the direction of the first transcript encountered in the
            ## transcript database; we did this also until version 1.3.5
            ## we now exclude these mutations! See below.
            trStrUnique <- trStrUnique[!duplicated(trStrUnique[,1]),]

        } else {
            ## uniqueTrDir == TRUE (default)

            # exclude cases where both transcription directions are valid
            # (e.g., ambiguous due to overlapping transcripts)
            trStrUnique <-
                trStrUnique[!(trStrUnique[,1] %in%
                              trStrUnique[duplicated(trStrUnique[,1])]),]
        }
        
        trPlus <- as.integer(trStrUnique[trStrUnique[,2]=="+",1])
        trMinus <- as.integer(trStrUnique[trStrUnique[,2]=="-",1])

        strands[trPlus] <- "+"
        strands[trMinus] <- "-"

    }


    # take reverse complement of all sequences without pyrimidine (C or T)
    # at center
    if (verbose) {
        cat(paste("Building reverse complement of sequences for which the",
                  "REF base is a purine (A or G).\n"))
    }

    revCompIndices <- grep("[^CT]", substr(as.character(seqs), (numBases%/%2)+1,
                                           (numBases%/%2)+1))

    seqs[revCompIndices] <- Biostrings::reverseComplement(seqs[revCompIndices])
    
    # also invert the corresponding strands; thus, "strands" no longer
    # contains the _genomic strand_ but becomes the "location" of the
    # pyrimidine (C,T) with respect to the _transcription direction_!
    # rationale:
    # a T in a gene on the "+" strand -> remains T and is in transcription
    #     direction ("+")
    # a T in a gene on the "-" strand -> remains T but is not on the
    #     transcribed strand (not in tr. dir., "-")
    # an A in a gene on the "+" strand -> translated to T but the T is not
    #     on the transcribed strand ("-")
    # an A in a gene on the "-" strand -> translated to T now this T is on
    #     the transcribed strand ("+")
    # so whenever we build the reverse compliment, we can also invert the
    # "strand" ...
    strands[revCompIndices] <- chartr("+-", "-+", strands[revCompIndices])

    # add to mutation table
    snvs <- cbind(snvs, as.character(seqs), strands)
    colnames(snvs)[seq((ncol(snvs)-1),ncol(snvs))] <- c("SEQ", "STRAND")

    # now, we have everything we need; start counting occurrences ...

    # get location of genotype info in sample specifications ("GT" in format)
    # [assume it's the same for all variants ...]

    gtIndex <- as.numeric(which(unlist(strsplit(snvs[1,"FORMAT"], ":",
                                                fixed=TRUE)) == "GT"))


    #
    # first: define empty data structure for this mutational load (=signature)
    # model and size
    #

    sigModel <- NULL
    if (type == "Alexandrov" || type == "full") {
        sigModel <- rep(0, 4^(numBases-1)*6*(1+as.numeric(trDir)))

        names(sigModel) <-
            buildSortedAlexandrovSignaturePatternList(numBases=numBases,
                                                      trDir=trDir)
        
    } else { # Shiraishi (table)
        sigModel <- matrix(0, ncol=6, nrow=(numBases+as.numeric(trDir)))

        sigModel <- setNames4ShiraishiTable(sigModel)

        # we need additional mapping to the colnames from single flanking
        # bases and transcription directions
        shColMapping <-
            c("[C>A]", "[C>G]", "[C>T]", "[T>A]", "[C>A]", "[C>G]")
        names(shColMapping) <-
            c(  "A",     "C",     "G",     "T",     "+",     "-"  )
    }

    #
    # now, iterate over samples/genomes!
    #

    sampleCols <- colnames(snvs)[!(colnames(snvs) %in%
                                   c("CHROM", "REF", "FORMAT",
                                     "POS", "ALT", "SEQ", "STRAND"))]

    if(verbose) {
        cat("Samples/genomes to be processed:\n")
        cat(paste0("  ",sampleCols,"\n"))
    }

    genomes <- plyr::alply(as.matrix(snvs[,sampleCols]), 2, function(sample) {
        ## apply to all samples (genotype column)

        if(verbose) {
            cat("Processing new genome:\n")
        }
        
        # get indices of mutations we need to process (genotype not NA and at
        # least heterozygous)
        if(verbose) {
            cat(paste("  Selecting mutations present in this genome",
                      "(genotype information), with defined ALT base and",
                      "without flanking N: "))
        }
        
        mutIndices <- grep("[^0/|]",
                           rapply(strsplit(sample,  ":", fixed=TRUE),
                                  function(x) { x[gtIndex] }))

        # but: exclude mutations with REF or ALT base different from
        # A, C, G or T (e.g. "N")
        mutIndices <- mutIndices[which(snvs[mutIndices, "REF"] %in%
                                       c("A","C","G","T"))]
        mutIndices <- mutIndices[which(snvs[mutIndices, "ALT"] %in%
                                       c("A","C","G","T"))]
        
        # exclude also those which have an N in their sequence pattern
        mutIndices <- mutIndices[grep("[^ACGT]", snvs[mutIndices, "SEQ"],
                                      invert=TRUE)]

        if(verbose) {
            cat(paste0(length(mutIndices)," variants left.\n"))
        }
        
        # but: in case we need transcriptional direction, remove all mutation
        # that don't have it
        # (i.e., that don't lie within genes)
        if (trDir) {
            if (verbose) {
                cat(paste("  Selecting mutations with information on",
                          "transcription direction (within genes): "))
            }

            mutIndices <- mutIndices[which(!is.na(snvs[mutIndices, "STRAND"]))]

            if(verbose) {
                cat(paste0(length(mutIndices)," variants left.\n"))
            }
        }


        # now build count occurrences for this sample/genome
        if (verbose) {
            cat(paste0("  Processing ", length(mutIndices),
                       " mutations to obtain frequencies.\n"))
        }

        genome <- sigModel
        
        for (mutIndex in mutIndices) {

            snv <- snvs[mutIndex,]
         
            ref <- snv["REF"]
            alt <- snv["ALT"]

            if (!ref %in% c("C","T")) {
                # need the complement also of the ALT base
                ref <- chartr("AG", "TC", ref)
                alt <- chartr("ACGT", "TGCA", alt)
            }
            basechange <- paste0("[",ref,">",alt,"]")
            
            if (type == "Shiraishi") {
                # increment base change
                genome["mut", basechange] <- genome["mut", basechange] + 1

                # increment counts for flanking bases
                flanking <-
                    unlist(strsplit(snv["SEQ"], ""))[-((numBases %/% 2)+1)]

                for (ii in seq_along(flanking)) {
                    genome[1+ii, shColMapping[flanking[ii]]] =
                        genome[1+ii, shColMapping[flanking[ii]]] + 1
                }

                # increment count for transcription direction, if required
                if (trDir) {
                    genome["tr", shColMapping[snv["STRAND"]]] =
                        genome["tr", shColMapping[snv["STRAND"]]] + 1
                }

            } else { # Alexandrov
                # construct element name (index of signature vector)
                index <- paste0(substr(snv["SEQ"], 1, (numBases %/% 2)),
                                basechange,
                                substr(snv["SEQ"], (numBases %/% 2)+2,
                                       numBases))
                if (trDir) {
                    index <- paste0(index, snv["STRAND"])
                }

                # increment this single element
                genome[index] <- genome[index] + 1
            }

        }

        # normalize the counts to get fractions
        genome <- genome/length(mutIndices)
        
        genome

    }) # genomes <- apply over samples ...

    names(genomes) <- sampleCols

    return(genomes)
    
}
