

######################
# internal functions #
######################

computeFrobeniusNorm <- function(A) { # compute the Frobenius norm of matrix A
    A = as.matrix(A) # in case A was a data.frame

    trace = sum(diag(A %*% t(A)))       # need the conjugate transpose here, 
                                        # but for real-valued matrices this is 
                                        # equivalent to the transpose!

    return(sqrt(trace))
}


computeRSS <- function(x, y) {
    if (length(x) != length(y)) {
        stop("In computeRSS: x and y need to have the same length!")
    }

    return(sum((x-y)^2))
}

determineTypeNumBasesAndTrDir <- function(mutData) {

    haveTrDir = NULL
    haveNumBases = NULL
    haveType = NULL
    
    if(class(mutData) == "matrix" && ncol(mutData) == 6) {
        # Shiraishi-type

        # number of bases of sequence pattern and transcription dir.
        haveNumBases = nrow(mutData)
        haveTrDir = FALSE
        if ( (nrow(mutData) %% 2) == 0 ) { # even number of rows,
                                           # must have transcript direction!
            haveNumBases = haveNumBases -1
            haveTrDir = TRUE
        }
    
        haveType = "Shiraishi"
        
    } else if (class(mutData) == "numeric") {
        # Alexandrov-type
 
        if ( ( (log10(length(mutData)/6)/log10(4)) %% 1 ) == 0 ) {
            
            # this must be without transcription direction
            haveNumBases = log10(length(mutData)/6)/log10(4) + 1
            haveTrDir = FALSE
            haveType = "Alexandrov"
            
        } else if ( ( (log10(length(mutData)/(6*2))/log10(4)) %% 1 ) == 0 ) {
            
            # this must be with transcription direction
            haveNumBases = log10(length(mutData)/(6*2))/log10(4) + 1
            haveTrDir = TRUE
            haveType = "Alexandrov"
        }
        
    } else{
        cat("Warning: cannot determine the number of bases and presence of transcription direction for signature/genome!\n")
    }

    return(list("type"=haveType, "numBases"=haveNumBases, "trDir"=haveTrDir))
}


buildSortedAlexandrovSignaturePatternList <- function(numBases, trDir) {
    seqpatterns <- NULL
    for (ii in (numBases-1):1) {
        seqpatterns <- paste0(seqpatterns, gl(4, 4^(ii-1), 4^(numBases-1),
                                              labels=c("A","C","G","T")))
    }
    
    # add the central base
    seqpatterns <-
        paste0(substr(seqpatterns,1,(numBases-1)/2),
               sort(rep(c("[C>A]","[C>G]","[C>T]","[T>A]","[T>C]","[T>G]"),
                        length(seqpatterns))),
               substr(seqpatterns,(numBases-1)/2+1,numBases-1))

    # double if we consider the transcription direction ...
    if (trDir) {
        seqpatterns <- paste0(c(seqpatterns,seqpatterns),
                              c(rep("+", length(seqpatterns)),
                                rep("-", length(seqpatterns))))
    }

    return(seqpatterns)
}


QPforSig <- function(counts, signatures, constrainToMaxContribution=FALSE,
                     tolerance=0.1) {

    if (!is.list(signatures)) {
        stop("Signatures must be a list for function QPforSig")
    }

    if (constrainToMaxContribution && (tolerance < 0 ||  tolerance > 1)) {
        stop("Tolerance must be between 0 and 1 when constraining the maximum contribution of signatures in function QPforSig")
    }

    # if necessary, convert signatures from matrices to vectors
    if (class(signatures[[1]]) == "matrix"
        || class(signatures[[1]]) == "data.frame") {
        
        signatures <-
            lapply(signatures, function(x){as.vector(t(as.matrix(x)))})
    }

    # if necessary, convert counts from matrix to vector
    if(class(counts) == "matrix" || class(counts) == "data.frame") {
        counts <- as.vector(t(as.matrix(counts)))
    }

    # the signatures should be stored in a matrix, one row is one signature
    sigMa <- do.call(rbind,signatures) 

    ### needed to make sure that the matrix (sigMa %*% t(sigMa)) is
    ### positive definite!
    sigMaSquared.nearPD =
        nearPD((sigMa %*% t(sigMa)), corr=FALSE, keepDiag=TRUE,
               ensureSymmetry=TRUE)$mat

    #Rinverse <- backsolve(chol(sigMa %*% t(sigMa)),diag(dim(sigMa)[1]))
    Rinverse <- backsolve(chol(sigMaSquared.nearPD),diag(dim(sigMa)[1]))

    if (constrainToMaxContribution) {
        # We want to contrain the maximum contribution each signature can
        # have (given some tolerance):
        # For each row in sigMa (each signature), we divide the counts from
        # the genome by the correpsonding fractions in the signature and take
        # the minimum over the signature, the signature can't have a higher
        # contribution than this ratio, but we add a tolerance value (given
        # that the data is noisy!)
        # Example: if in the genome 30% of variants have a specific feature
        # (e.g. a specific flanking base), i.e., the counts score is 0.3, and
        # 60% (0.6) of the variants produced by a process/signature have this 
        # feature, then the signature can have contributed up to 0.3/0.6 = 0.5
        # (50%) of the genome's variants. 
        # With the default tolerance of 0.1, we would allow the quadratic
        # programming approach to assign up to 0.5+0.1 (60%) of the variants
        # to this signature

        maxContributions <-
            apply(sigMa, 1, function(x) {
                min((counts/x)[!is.na(counts/x)]) + tolerance })

        # limit maximum contribution to 100%
        maxContributions[maxContributions>1] = 1


        # for the additional contraints, we use negative values because we 
        # require contribution/weight w <= w_max, which we can implement 
        # as -w >= -w_max (because quadprog uses Ax >= b, so we have to
        # use >= instead of <=!)
        Amat <- cbind(rep(1,dim(Rinverse)[1]),
                      diag(dim(Rinverse)[1]), -diag(length(maxContributions)))
        bvec <- c(1, rep(0,dim(Rinverse)[1]), -maxContributions)
    }
    else { # !constrainToMaxContribution
        Amat <- cbind(rep(1,dim(Rinverse)[1]), diag(dim(Rinverse)[1]))
        bvec <- c(1, rep(0,dim(Rinverse)[1]))
    }

    dvec <- (counts) %*% t(sigMa)
    resQP <- quadprog::solve.QP(Dmat=Rinverse, dvec=dvec, Amat=Amat,
                                bvec=bvec, meq=1, factorized=TRUE)
    resQP$solution[resQP$solution<0]=0  # there can be some values below 0
                                        # which should actually be 0 (very 
                                        # close to 0 in any case, e.g., -10^-15)

    return(resQP$solution)
}



#
# convert a single Alexandrov signature to a Shiraishi signature
#

convAlx2Shi <- function(x) {

    strLen <- unique(nchar(names(x)))   # example: A[T>G]C

    haveTrDir = FALSE
    if(substr(names(x)[1], strLen, strLen) %in% c("+","-")) {
        # do we have info on transcription direction?
        haveTrDir = TRUE
    }

    basesAdj <- (strLen-5-as.numeric(haveTrDir))/2

    # first row: base change frequencies
    shSig <- as.vector(tapply(x, substr(names(x),basesAdj+2,basesAdj+4), sum))

    # upstream base frequencies
    for (ii in 0:(basesAdj-1)) {
        shSig <- rbind(shSig,
            c(as.vector(tapply(x, substr(names(x),1+ii,1+ii), sum)), rep(0,2))
                       )
    }

    # downstream base frequencies
    for (ii in 0:(basesAdj-1)) {
        shSig <- rbind(shSig,
                       c(as.vector(tapply(x, substr(names(x),
                                                    strLen-ii
                                                        -as.numeric(haveTrDir),
                                                    strLen-ii
                                                        -as.numeric(haveTrDir)),
                                          sum)), rep(0,2))
                       )
    }

    if(haveTrDir) {
        shSig <- rbind(shSig,
                       c(as.vector(tapply(x, substr(names(x),strLen,strLen),
                                          sum)[c("+","-")]), rep(0,4))
                       )
    }

    rownames(shSig) <- NULL

    return(shSig)
}



#
# build genome data structures (same as signatures) and fill with mutation data
#

buildGenomesFromMutationData <- function(snvs, numBases, type, trDir,
                                         refGenome, transcriptAnno, verbose) {

    if (type == "independent") {
        type = "Shiraishi"
    } else if (type == "full") {
        type = "Alexandrov"
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
        cat("Comparing specification of mitochondrion in mutation data (M or MT?) to reference genome and adjusting if necessary.\n")
    }

    # check whether the mitochondrion is specified as chromosome "MT" or "M"

    refMT = ""
    if (length(grep("MT", GenomicRanges::seqnames(refGenome@seqinfo))) == 1) { 
        refMT = "MT"
    } else if (length(grep("M", GenomicRanges::seqnames(refGenome@seqinfo)))
               == 1) {
        refMT = "M"
    } # else empty, no clear mitochondrion found

    #second: VCF file
    vcfMT = ""
    if (length(grep("MT", unique(snvs[,"CHROM"]))) == 1) { # having MT
        vcfMT = "MT"
    } else if (length(grep("M", unique(snvs[,"CHROM"]))) == 1) {
        vcfMT = "M"
    } # else empty, no clear mitochondrion found

    # third: check and adjust VCF data, if necessary
    if (vcfMT != "" && vcfMT != refMT) {
        if (refMT != "") {
            # invert
            snvs[,"CHROM"] = gsub(vcfMT, refMT, snvs[,"CHROM"])
        }
    }


    
    # check if chromosome IDs match what we have in the reference genome
    # (with or withour prefix "chr")
    if(verbose) {
        cat("Checking chromosome IDs of mutation data (1 or chr1?) and ajusting if necessary.\n")
    }

    if (!all(unique(snvs[,"CHROM"]) %in%
             GenomicRanges::seqnames(refGenome@seqinfo))) {
        
        # some problem, first try adding "chr"
        if(all(paste("chr",unique(snvs[,"CHROM"]),sep="") %in%
               GenomicRanges::seqnames(refGenome@seqinfo))) {
            # add "chr"
            snvs[,"CHROM"] <- paste("chr",snvs[,"CHROM"], sep="")

        } else if (all(gsub("^chr", "", snvs[,"CHROM"]) %in%
                       GenomicRanges::seqnames(refGenome@seqinfo))) {
            # try removing "chr"
            snvs[,"CHROM"] <- gsub("^chr", "", snvs[,"CHROM"])
        }
        
    }

    
    # get sequences from the reference genome
    if (verbose) {
        cat("Extracting sequences with flanking bases of mutations from the reference genome.\n")
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
        stop("Reference (REF) bases in VCF do not match the specified reference genome!")
    }

    
    # set strand info to NA
    strands <- rep(NA, length(seqs))
    
    # get strand info if we need it
    if(trDir) {
        if(verbose) {
            cat("Getting information on transcription directions for mutations within genes.\n")
        }

        tr <- GenomicFeatures::transcripts(transcriptAnno)
        trHit <- GenomicRanges::findOverlaps(gr, tr, ignore.strand = FALSE)
        trStr <- cbind(S4Vectors::queryHits(trHit),
                       as.character(BiocGenerics::strand(
                           tr[S4Vectors::subjectHits(trHit)])
                                    )
                       )
        trStrUnique <- unique(trStr[trStr[,2] != "*",],MARGIN=1)
        trStrUnique <- trStrUnique[!duplicated(trStrUnique[,1]),]

        trPlus <- as.integer(trStrUnique[trStrUnique[,2]=="+",1])
        trMinus <- as.integer(trStrUnique[trStrUnique[,2]=="-",1])

        strands[trPlus] <- "+"
        strands[trMinus] <- "-"
    }


    # take reverse complement of all sequences without pyrimidine (C or T)
    # at center
    if (verbose) {
        cat("Building reverse complement of sequences for which the REF base is a purine (A or G).\n")
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
    colnames(snvs)[(ncol(snvs)-1):ncol(snvs)] <- c("SEQ", "STRAND")

    # now, we have everything we need; start counting occurrences ...

    # get location of gentype info in the sample specifications ("GT" in format)
    # [assume it's the same for all variants ...]

    gtIndex <- as.numeric(which(unlist(strsplit(snvs[1,"FORMAT"], ":",
                                                fixed=TRUE)) == "GT"))


    #
    # first: define empty data structure for this mutational load (=signature)
    # model and size
    #

    sigModel <- NULL
    if (type == "Alexandrov") {
        sigModel <- rep(0, 4^(numBases-1)*6*(1+as.numeric(trDir)))

        names(sigModel) <-
            buildSortedAlexandrovSignaturePatternList(numBases=numBases,
                                                      trDir=trDir)
        
    } else { # Shiraishi (table)
        sigModel <- matrix(0, ncol=6, nrow=(numBases+as.numeric(trDir)))

        colnames(sigModel) <- c("[C>A]","[C>G]","[C>T]","[T>A]","[T>C]","[T>G]")

        rnames <- c("mut", c(-(numBases%/%2):-1), c(1:(numBases%/%2)))
        if(trDir) {
            rnames <- c(rnames, "tr")
        }

        rownames(sigModel) <- rnames

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
            cat("  Selecting mutations present in this genome (genotype information), with defined ALT base and without flanking N: ")
        }
        
        mutIndices <- grep("[^0/|]",
                           rapply(strsplit(sample,  ":", fixed=TRUE),
                                  function(x) { x[gtIndex] }))

        # but: exclude mutations with REF or ALT base different from
        # A, C, G ot T (e.g. "N")
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
                cat("  Selecting mutations with information on transcription direction (within genes): ")
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

            snv = snvs[mutIndex,]
         
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
                genome["mut", basechange] = genome["mut", basechange] + 1

                # increment counts for flanking bases
                flanking <-
                    unlist(strsplit(snv["SEQ"], ""))[-((numBases %/% 2)+1)]

                for (ii in 1:length(flanking)) {
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
                genome[index] = genome[index] + 1
            }

        }

        # normalize the counts to get fractions
        genome <- genome/length(mutIndices)
        
        genome

    }) # genomes <- apply over samples ...

    names(genomes) <- sampleCols

    return(genomes)
    
}



getBestDecomp4Ksignatures <- function(genome, signatures, k,
                                      constrainToMaxContribution=FALSE,
                                      tolerance=0.1) {

    # get list of all possible combinations of k signatures
    sigCombn <- combn(1:length(signatures), k)
    sigCombnNames <-
        apply(sigCombn, 2, function(x) {
            paste(names(signatures)[x], collapse="|") } )

    sigCombn <- split(t(sigCombn), seq(ncol(sigCombn)))
    names(sigCombn) <- sigCombnNames
            
    decompTmp <- lapply(sigCombn, function(sigIndices) {

        # signatures to be used
        sigs <- signatures[sigIndices]
                
        # decompose for these signatures
        QPforSig(genome, sigs,
                 constrainToMaxContribution=constrainToMaxContribution,
                 tolerance=tolerance)
    } )

    #explVarTmp <- plyr::laply(decompTmp, function(decomp) {
    explVarTmp <- sapply(names(decompTmp), function(sigsNames) {

        # signatures to be used
        sigs <- signatures[unlist(strsplit(sigsNames, "|", fixed=TRUE))]
                
        # determine the explained variance for this decomposition
        computeExplainedVariance(decompTmp[sigsNames], sigs, list(genome))
    } )            
    names(explVarTmp) <- names(decompTmp)

    # which was the best explained variance for this number of signatures?
    explVar <- max(explVarTmp)

    # which decomposition achieved this explained variance? If multiple do,
    # select one at random
    if (length(which(explVarTmp == explVar)) > 1) {
        explVarIndex <- sample(which(explVarTmp == explVar),1)
    } else {
        # only on, using sample would lead to wrong behavior!
        explVarIndex <- which(explVarTmp == explVar)
    }
    
    sigList <- unlist(strsplit(names(explVarTmp)[explVarIndex], "|",
                               fixed=TRUE))
    decomposition <- decompTmp[[explVarIndex]]
    names(decomposition) <- sigList
    
    list(k = k,
         explVar = explVar,
         sigList = sigList,
         decomposition = decomposition)
}



addBestSignatureToSubset <- function(genome, signatures, subset,
                                     constrainToMaxContribution=FALSE,
                                     tolerance=0.1) {

    subsetIndices <- which(names(signatures) %in% subset)
    restIndices <- which(!(names(signatures) %in% subset))

    # get all possible combinations of subset + one assitional signature
    sigCombn <-
        lapply(restIndices, function(ix) { sort(c(subsetIndices, ix)) } )
    
    names(sigCombn) <-
        unlist(lapply(sigCombn, function(set) {
            paste(names(signatures)[set],collapse="|") } ))

    #print(sigCombn)
    
    decompTmp <- lapply(sigCombn, function(sigIndices) {

        # signatures to be used
        sigs <- signatures[sigIndices]
                
        # decompose for these signatures
        QPforSig(genome, sigs,
                 constrainToMaxContribution=constrainToMaxContribution,
                 tolerance=tolerance)
    } )

    #print(decompTmp)
    
    explVarTmp <- sapply(names(decompTmp), function(sigsNames) {

        # signatures to be used
        sigs <- signatures[unlist(strsplit(sigsNames, "|", fixed=TRUE))]
                
        # determine the explained variance for this decomposition
        computeExplainedVariance(decompTmp[sigsNames], sigs, list(genome))
    } )            
    names(explVarTmp) <- names(decompTmp)

    #print(explVarTmp)
    
    # which was the best explained variance for this number of signatures?
    explVar <- max(explVarTmp)

    # which decomposition achieved this explained variance? If multiple do,
    # select one at random
    if (length(which(explVarTmp == explVar)) > 1) {
        explVarIndex <- sample(which(explVarTmp == explVar),1)
    } else {
        # only on, using sample would lead to wrong behavior!
        explVarIndex <- which(explVarTmp == explVar)
    }
    
    sigList <-
        unlist(strsplit(names(explVarTmp)[explVarIndex], "|", fixed=TRUE))

    decomposition <- decompTmp[[explVarIndex]]
    names(decomposition) <- sigList
    
    list(k = length(subset)+1,
         explVar = explVar,
         sigList = sigList,
         decomposition = decomposition)
}
