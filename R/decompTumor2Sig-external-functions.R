

####################
# public functions #
####################


### given a tumor genome (or list of genomes) and a set of signatures,
### decompose each individual tumor

decomposeTumorGenomes <- function(genomes, signatures,
                                  minExplainedVariance=NULL,
                                  minNumSignatures=2, maxNumSignatures=NULL,
                                  greedySearch=FALSE, 
                                  constrainToMaxContribution=FALSE,
                                  tolerance=0.1, verbose=FALSE) {
    
    # input: gnomes = either a list of genome mutation counts as matrices or
    #                  vectors, or a single genome (matrix or vector)
    #        signatures = a list of signatures (matrices or vectors)
    #        [need to have the same format]

    # constrainToMaxContribution and tolerance:
    # We want to contrain the maximum contribution each signature can have
    # (given some tolerance): For each row in sigMa (each signature), we
    # divide the counts from the genome by the correpsonding fractions in
    # the signature and take the minimum over the signature, the signature
    # can't have a higher contribution than this ratio, but we add a tolerance
    # value (given that the data is noisy!)
    # Example: if in the genome 30% of variants have a specific feature
    # (e.g., a specific flanking base), i.e., the counts score is 0.3, and
    # 60% (0.6) of the variants produced by a process/signature have this 
    # feature, then the signature can have contributed up to 0.3/0.6 = 0.5 (50%)
    # of the genome's variants. 
    # With the default tolerance of 0.1, we would allow the quadratic
    # programming approach to assign up to 0.5+0.1 (60%) of the variants to
    # this signature

    # minExplainedVariance and minNumSignatures and maxNumSignatures:
    # If minExplainedVariance is specified, the minimum number of signatures
    # required will be determined such that the minimum threshold for
    # explained variance is satsified. The tested numbers of signatures range
    # from minNumSignatures (default: 2) to maxNumSignatures (default: all)

    # if greedySearch=TRUE then _not_ all possible combinations of
    # minNumSignatures to maxNumSignatures signatures will be checked. Instead,
    # first all possible combinations for exactly minNumSignatures will be
    # checked, then the next best signature will be added (maximum increase
    # in explained variability) until minExplainedVariance is reached
    # (or maxNumSignatures is exceeded).
    
    # If minExplainedVariance is NULL, then exactly maxNumSignatures
    # signatures will be taken (default: all)

    # Returns: exposure(s), but only if they satisfy minExplainedVariance or
    # if not minimum explained variance is requested.

    
    if (!is.list(signatures)) {
        stop("Parameter 'signatures' must be a list of signature objects!")
    }
    if (!is.data.frame(signatures[[1]]) & !is.matrix(signatures[[1]])
        & !(is.vector(signatures[[1]]) & is.numeric(signatures[[1]]))) {

        stop("Signatures must be data.frames, matrices or numeric vectors!")
    }

    # if maximum number of signatures is not defined, set it to the total
    # number of signatures
    if(is.null(maxNumSignatures)) {
        maxNumSignatures = length(signatures)
    }

    # if we don't require to find the minimum number of signatures for which
    # we exceed the minimum explained variance, we simply take the maximum
    # number of signatures (usually the default: all)
    # This way we make sure that we use exaclty k = maxNumSignatures
    if (is.null(minExplainedVariance)) {
        minNumSignatures = maxNumSignatures
        
    } else if(class(minExplainedVariance)!="numeric"
              || minExplainedVariance<0 || minExplainedVariance>1) {

        stop("minExplainedVariance must be NULL or between 0 and 1!")
    }

    
    # is the signatures are unnamed, name them by enumerating them
    if(is.null(names(signatures))) {
        names(signatures) <- paste0("sign_",1:length(signatures))
    }

    
    if (!is.logical(greedySearch)) {
         stop("greedySearch must be logical (TRUE or FALSE)!")
    }

    if (!is.logical(constrainToMaxContribution)) {
         stop("constrainToMaxContribution must be logical (TRUE or FALSE)!")
    }

    if (constrainToMaxContribution && (tolerance < 0 ||  tolerance > 1)) {
        stop("tolerance must be between 0 and 1 when constraining the maximum contribution of signatures (default: 0.1)!")
    }

    if (!is.list(genomes)) { # a single genome, make it a list of one genome
        genomes = list(genomes)
        #name(genomes) = "tumor_genome"
    }

    if (is.null(names(genomes))) { # numbering of genomes if they are not named!
                                   # (we need the names for the results)
        names(genomes) = paste0("genome_", as.character(1:length(genomes)))
    }

    # initialize an empty list with all elements set to NULL
    decompositions  = vector("list", length(genomes))  
    
    for (g in 1:length(genomes)) {
        # evaluate each genome individually

        if(verbose) {
            cat(paste0("Decomposing genome ",g," (",names(genomes)[g],")"))
        }
        
        if (!is.data.frame(genomes[[g]]) & !is.matrix(genomes[[g]])
            & !(is.vector(genomes[[g]]) & is.numeric(genomes[[g]]))) {
            
            stop("genomes must be data.frames, matrices or numeric vectors!")
        }
        
        counts = genomes[[g]]

        # first, check its format
        if (length(counts) != length(as.vector(as.matrix(signatures[[1]])))) {
            stop("Formats of genomes and signatures must match!")
        }

        
        # iterating for all possible numbers of signatures from
        # minNumSignatures to maxNumSignatures
        decompTmpList = list()

        if(!greedySearch) {
            # default: full search; all possible combinations!
            for (k in minNumSignatures:maxNumSignatures) {

                if(verbose) {
                    cat(paste(" with", k, "signatures ...\n"))
                }
            
                decompTmpList[[length(decompTmpList)+1]] <-
                    getBestDecomp4Ksignatures(counts, signatures, k,
                                              constrainToMaxContribution,
                                              tolerance)

                if (!is.null(minExplainedVariance)) {
                    # test the decomposition to check if we found one that
                    # explains at least minExplainedVariance of the variance
                    # of the observed mutation data of the genome

                    if (decompTmpList[[length(decompTmpList)]]$explVar
                        >= minExplainedVariance) {
                        
                        break
                    }
                }
            }
        } else {
            # greedy search! start with best combination of minNumSignatures;
            # then add one at a time
            k = minNumSignatures

            if(verbose) {
                cat(paste(" with", k, "signatures ...\n"))
            }

            decompTmpList[[length(decompTmpList)+1]] <-
                getBestDecomp4Ksignatures(counts, signatures, k,
                                          constrainToMaxContribution,
                                          tolerance)

            while(decompTmpList[[length(decompTmpList)]]$explVar
                       < minExplainedVariance
                  && k < maxNumSignatures) {

                if(verbose) {
                    cat(paste(" adding signature", k+1, "...\n"))
                }

                haveSubset <- decompTmpList[[length(decompTmpList)]]$sigList
                decompTmpList[[length(decompTmpList)+1]] <-
                    addBestSignatureToSubset(counts, signatures, haveSubset,
                                             constrainToMaxContribution,
                                             tolerance)
                k = k + 1
            }
        }

        haveExplVar <- decompTmpList[[length(decompTmpList)]]$explVar
        
        if(verbose) {
            cat(paste(" explained variance:", haveExplVar))
        }
        
        # The decomposition we have last, should be the correct one, either
        # because the k signatures exceed the threshold minExplainedVariance,
        # or because no such threshold was required
        # However, if no combination of signatures reached the specified
        # threshold, we do not report the result!

        if (is.null(minExplainedVariance)
            || haveExplVar >= minExplainedVariance) {
            # threshold not specified or satisfied

            ##previous version: return only used signatures
            #decompositions[[g]] =
            #    decompTmpList[[length(decompTmpList)]]$decomposition

            ##new behavior: return _all_ signatures but set unused signatures to NA

            # empty vector of NAs
            decompositions[[g]] = rep(NA, length(signatures))
            names(decompositions[[g]]) = names(signatures)

            # set used signatures to their exposure/contribution
            decompositions[[g]][decompTmpList[[length(decompTmpList)]]$sigList] = decompTmpList[[length(decompTmpList)]]$decomposition

            if(verbose) {
                cat("\n")
            }
        } else {
            # threshold specified not satisfied; do not report the results
            # (set to NULL instead)

            if(verbose) {
                cat(paste(" <", minExplainedVariance,"(rejected)\n"))
            }
        }

        names(decompositions)[g] = names(genomes)[g]

    }
    
    return(decompositions)
}



### convert tumor genomes that have been read with pmsignature

getGenomesFromMutationFeatureData <- function(countData, normalize=TRUE) {
    # convert MutationFeatureData (e.g. loaded with pmsignature's function
    # readMPFile) to a list of genome matrices

    # IMPORTANT: set normalize=FALSE only if you want to see full counts,
    # but USE ONLY NORMALIZED count tables for predicting/computing exposures,
    # i.e., the contributions of mutational signatures (because these are
    # normalized, too)
    
    if (class(countData) != "MutationFeatureData") {
        stop("countData must be an object of type MutationFeatureData (as produced by pmsignature's readMPFile)")
    }
    
    # if the user set trDir = TRUE, we need flankinBasesNum+1 rows,
    # otherwise only the flankingBasesNum
    if ((getElement(countData, "transcriptionDirection"))==TRUE) {
        rows <- as.numeric(getElement(countData, "flankingBasesNum")+1)
    } else {
        rows <- as.numeric(getElement(countData, "flankingBasesNum")) 
    }

    tableList <- list()

    # get list of samples (numbers)
    samples <- unique(getElement(countData,"countData")[2,])   

    for (l in 1:length(samples)){  # for each sample l

        # initialize the new count table/matrix with 0
        countTable <- matrix(rep(0, 6*rows), nrow=rows, ncol=6) 

        # in oneSample we extract the count data relative to sample l
        oneSample <-
            getElement(countData,
                       "countData")[,which(getElement(countData,
                                                      "countData")[2,]==l)]

        for (i in 1:ncol(oneSample)) { # for each i in the oneSample vector

            # index of the mutation type in featureVectorList
            mutTypeIndex <- oneSample[1,i] 

            # get the features for this mutation type
            # (corresponds to the column indices to be incremented for each
            # row in the count table)
            v <- getElement(countData,"featureVectorList")[,mutTypeIndex] 

            
            for (k in 1:rows){ #for each k in the length of our rows
                
                # in our countTable at position k and v[k] we store the
                # number of occurrences which is stored in the third row
                # and i-th column of oneSample
                countTable[k,v[k]]= countTable[k,v[k]] + oneSample[3,i] 
            }
        }

        # normalize the count table?
        if (normalize) {
            countTable <- countTable/(sum(countTable)/nrow(countTable))
        }
        
        tableList[[l]] <- countTable  # put each countTable in the tableList

    }

    # keep the samples' names/IDs!
    names(tableList) <- getElement(countData, "sampleList")

    return(tableList)
}


### load mutations from a VCF file

loadGenomesFromVCF <- function(file, numBases=5, type="Shiraishi", trDir=TRUE,
    refGenome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
    transcriptAnno=
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
    verbose=TRUE) {

    # read mutation data
    if(verbose) {
        cat("Loading mutations and genotype information from VCF file:\n")
    }
    vcf <- read.vcfR(file, verbose=verbose, checkFile=TRUE, convertNA=TRUE)

    chrIndex <- which(colnames(vcf@fix) == "CHROM")
    posIndex <- which(colnames(vcf@fix) == "POS")
    refIndex <- which(colnames(vcf@fix) == "REF")
    altIndex <- which(colnames(vcf@fix) == "ALT")

    # get only SNVs (one REF base and one ALT base)
    snvRows <-
        (nchar(vcf@fix[,refIndex]) == 1) & (nchar(vcf@fix[,altIndex]) == 1)

    # basic variant information (chr, pos, ref, alt)
    snvs <- vcf@fix[snvRows, c(chrIndex,posIndex,refIndex,altIndex)]

    # genotype information
    if (dim(vcf@gt)[1] > 0) {
        # if we do have genotype information, add it!
        snvs <- cbind(snvs, vcf@gt[snvRows,])   # 1st is "FORMAT",
                                                # following are samples/genomes
    } else {
        # take all variants as originating from the same sample/genome,
        # create dummy genotype
        snvs <- cbind(snvs, matrix(c("GT", "1/1"),
                                   nrow=length(which(snvRows)),
                                   ncol=2, byrow=TRUE)
                      )
        colnames(snvs)[(ncol(snvs)-1):ncol(snvs)] =
            c("FORMAT", "variants_without_genotype_info")
    }

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
        cat("Done loading genomes.\n")
    }

    return(genomes)

}



### load mutations from a MPF file (Mutation Position Format)

loadGenomesFromMPF <- function(file, numBases=5, type="Shiraishi", trDir=TRUE,
    refGenome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
    transcriptAnno=
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
    verbose=TRUE) {

    # read mutation data
    if(verbose) {
        cat("Loading mutations and patient/sample IDs from MPF file:\n")
    }
    mpf <- as.matrix(read.table(file, header=FALSE, row.names=NULL, sep="\t"))
    colnames(mpf) = c("SAMPLE", "CHROM", "POS", "REF", "ALT")

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
        cat("Collapsing variant information by mapping multiple samples to unique variants.\n")
    }
    
    sampleIDlist <- sort(unique(mpf[,smpIndex]))

    snvs <- plyr::aaply(mpf, 1, function(x) {
        gtVec = c("GT", rep(NA, length(sampleIDlist)))
        names(gtVec) = c("FORMAT", sampleIDlist)
        gtVec[x[smpIndex]] = "1/1"
        
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
        cat("Done loading genomes.\n")
    }

    return(genomes)

}



### load mutations from a VCF file

convertGenomesFromVRanges <- function(vranges,
    numBases=5, type="Shiraishi", trDir=TRUE,
    refGenome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
    transcriptAnno=
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
    verbose=TRUE) {

    if (class(vranges) != "VRanges") {
        stop("vranges must be an object of type VRanges!")
    }
    
    # read mutation data
    if(verbose) {
        cat("Extracting mutations and sample/genotype information from VRanges:\n")
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




### load a set of Shiraishi signatures from flat files

loadShiraishiSignatures <- function(files) {
    # load a set of Shiraishi signatures from flat files
    # files can be a single file name of vector of file names as returned
    #   by list.files(path="", pattern="" ....)

    if (class(files) == "list") {
        files = unlist(files)
    }

    if (class(files) != "character") {
        stop("Parameter 'files' must be a filename or a list/vector of file names!")
    }
    
    fnames = as.vector(files)
    
    sigList = list()

    for (fn in fnames) {
        sigmatrix <- as.matrix(read.table(fn, header=FALSE, row.names=NULL))
        colnames(sigmatrix) = NULL

        sigList[[length(sigList)+1]] <- sigmatrix
    }

    if (length(unique(basename(fnames))) == length(fnames)) {
        # base file names (without paths) are unique, take as signature names
        names(sigList) = basename(fnames)
    } else {
        # take the unique path/file names as signature names
        names(sigList) = fnames
    }
    
    return(sigList)
}


### convert a set of signatures from the result of pmsignature

getSignatureListFromEstimatedParameters <- function(Param) {
    # convert the signatures from an "EstimatedParameters" (Param),
    # as obtained from pmsignature, to a list of signature matrices or vectors

    if (!requireNamespace("pmsignature", quietly=TRUE)) {
        stop("Function getSignatureListFromEstimatedParameters requires the package pmsignature to be installed.")
    }
    
    if (class(Param) != "EstimatedParameters") {
        stop("Param must be an object of type EstimatedParameters (as produced by pmsignature's getPMSignature)")
    }

    sigList = list()

    # number of signatures
    numSigs = getElement(Param, "signatureNum")
    
    # if one of them was background, the true number of signatures is one less!
    numSigs = numSigs - as.numeric(getElement(Param, "isBackGround"))

    for (sig in 1:numSigs) {
        sigList[[length(sigList)+1]] <-
            pmsignature::getSignatureValue(Param, sig)
    }

    return(sigList)
}



# load Alexandrov signatures (COSMIC format) from a file

loadAlexandrovSignatures <- function(file="http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt") {
    # load a set of Alexandrov signatures from a tab-seperated flat file in the
    # format provided by COSMIC Mutational Signatures, available at: 
    # http://cancer.sanger.ac.uk/cancergenome/assets/
    #                   signatures_probabilities.txt (this is taken by default)
 
    if (class(file) != "character") {
        stop("Parameter 'file' must be a filename or URL!")
    }

    # load all data in one table
    sigmatrix <- as.matrix(read.table(file, header=TRUE,
                                      row.names=NULL, sep="\t"))

    sigList = list()
    sigNames = c()

    for (colId in 4:ncol(sigmatrix)) {
        if (!all(is.na(sigmatrix[,colId]))) { # signature is defined
            sigVec = as.numeric(sigmatrix[,colId])
            names(sigVec) = sigmatrix[,3]

            # make sure we sort the vector corretly (the official file has
            # another sorting): first check number of bases and presence of
            # transcription direction
            sigFeatures <- determineTypeNumBasesAndTrDir(sigVec)
            if (is.null(sigFeatures$type) || sigFeatures$type != "Alexandrov") {
                stop("Wrong number of patterns for an Alexandrov signature!")
            }
            # reorder the signature vector so that base changes stay together
            sigVec <-
                sigVec[buildSortedAlexandrovSignaturePatternList(
                    numBases=sigFeatures$numBases,
                    trDir=sigFeatures$trDir)]
            
            # the vector is fine now, save in list of multiple signatures
            sigList[[length(sigList)+1]] <- sigVec
            sigNames = c(sigNames, colnames(sigmatrix)[colId])
        }
    }

    names(sigList) = sigNames

    return(sigList)
}



# [EXPERIMENTAL] convert Alexandrov signatures to Shiraishi signatures

convertAlexandrov2Shiraishi <- function(signatures) {

    if (!class(signatures) == "list") {
        # if the user specified a single signature, make it a list
        signatures = list(signatures)
    }

    shSignatures <- lapply(signatures, convAlx2Shi)

    return(shSignatures)
}


# Plot mutational data for a signature or a genome (Shiraishi-type or
# Alexandrov-type)

plotMutationDistribution <- function(mutData) {

    if (!requireNamespace("pmsignature", quietly=TRUE)) {
        stop("Function plotMutationDistribution requires the package pmsignature to be installed.")
    }

    if (class(mutData) == "list") {
        cat("Warning: object 'mutData' passed to plotMutationDistribution is a list. Taking only the first element and trying to plot it ...\n")
        mutData = mutData[[1]]  # take first object from list
    }
    
    # build an "EstimatedParameters" object for visualization with pmsignature
    pmParam = new(Class = "EstimatedParameters")

    # try to determine the number of bases and whether the transcription
    # direction was considered
    dataFeatures <- determineTypeNumBasesAndTrDir(mutData)
    if (is.null(dataFeatures$type)) {
        stop("Object 'mutData' must be a single Shiraishi- or Alexandrov signature!")
    }

    haveType = dataFeatures$type

    pmParam@flankingBasesNum = as.integer(dataFeatures$numBases)
    pmParam@transcriptionDirection = dataFeatures$trDir


    if(haveType == "Shiraishi") {
        # Shiraishi-type

        pmParam@type = "independent" # Shiraishi-type
        pmParam@signatureFeatureDistribution =
            array(0, dim=c(1, nrow(mutData), ncol(mutData)) )
        
    } else if (haveType == "Alexandrov") {
        # Alexandrov-type
        pmParam@type = "full"
        pmParam@signatureFeatureDistribution =
            array(0, dim=c(1, length(mutData), 1) )
    }
    
    pmParam@signatureFeatureDistribution[1,,] = mutData  # data to visualize

    pmsignature::visPMSignature(pmParam, 1, isScale=TRUE)
}


plotDecomposedContribution <- function(decomposition, signatures=NULL,
                                       removeNA=TRUE) {

    # dummy declarations to avoid NOTEs by R CMD check together with
    # ggplot2 syntax
    Signatures <- NULL
    Exposures <- NULL
    # end of dummy stuff
    
    if (class(decomposition) == "list") {
        
        cat("Warning: object 'decomposition' passed to plotDecomposedContribution is a list. Taking only the first element and trying to plot ...\n")

        decomposition = decomposition[[1]]
    }

    # determine signature names

    if (!is.null(names(decomposition))) {
        
        # first guess: directly from the exposure vector
        sigNames = names(decomposition)
        
    } else if (!is.null(signatures) && !is.null(names(signatures))) {
        
        # second guess: from passed signatures
        sigNames = names(signatures)
        
    } else {
        
        # last solution: just number them
        sigNames = paste0("sign_", 1:length(decomposition))
        
    }

    
    if (length(sigNames) != length(decomposition)) {
        # might happen if the signatures aren't associated with the exposure vector
        stop("The number of exposures is different from the number of signatures!")
    }


    if(removeNA) {
        # ignore NAs; might be due to a greedy search
        sigNames <- sigNames[!is.na(decomposition)]
        decomposition <- decomposition[!is.na(decomposition)]
    }
    
    # construct data frame for exposures
    df <- data.frame(Signatures=factor(sigNames, levels=sigNames),
                     Exposures=decomposition)

    ggplot(data=df, aes(x=Signatures, y=Exposures)) + xlab(NULL) +
        ylab("Exposures (percent contribution)") +
        geom_bar(stat="identity", width=0.75, color="black", fill="steelblue") +
        geom_text(aes(label=round(Exposures,digits=2)), vjust=1.6,
                  color="white", size=2) +
        theme(panel.background = element_rect(fill="white", colour="black"),
              panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                  colour = "lightgray"),
              panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                  colour = "lightgray"),
              axis.text.x = element_text(face="bold", angle=90, vjust=0.5))

}



# compute the explained variance for the estimated contributions (exposures)
# of a set of signatures when compared to the observed genomes

computeExplainedVariance <- function(exposures, signatures, genomes) {
    if ((is.vector(genomes) & is.numeric(genomes))
        || is.matrix(genomes) || is.data.frame(genomes)) {
        # this is only one genome, use a list nonetheless for later iteration
        genomes = list(genomes)
    }
     if (class(exposures) == "numeric") {
        # same as for genomes
        exposures = list(exposures)
    }
    if (class(signatures) != "list") {
        stop("Parameter signatures must be a list object.")
    }

    # check that we have as many exposure estimates as genomes
    if(length(genomes)!=length(exposures)) {
        stop("Number of exposure vectors different from number of genomes.")
    }

    expvar <- NULL

    # for each genome, compare estimated exposures
    for (g in 1:length(genomes)) {

        if (is.null(exposures[[g]])) {
            # can be NULL if minExplainedVariance wasn't reached in a
            # subset search

            expvar <- c(expvar, NA)
            next
        }

        # set NA in exposures to 0
        exposures[[g]][is.na(exposures[[g]])] = 0
        
        # first compute estimated genome from exposures and signatures
        if (length(signatures) != length(exposures[[g]])) {
            stop("Number of exposures different from number of signatures.")
        }

        first = TRUE        
        for (s in 1:length(signatures)) {

            if (first) {
                estgenome <- signatures[[s]] * as.numeric(exposures[[g]][s])
                first = FALSE
            } else {
                estgenome <- estgenome + (signatures[[s]] *
                                              as.numeric(exposures[[g]][s]))
            }
        }

        # now compute the residual sum of squares (RSS) between estimated and
        # observed genome
        #rss <- sum((estgenome - genomes[[g]])^2)
        rss <- computeRSS(estgenome, genomes[[g]])

        # ... and the total sum of squares (TSS) between the observed genome
        # and the "average/unvaried" genome.
        # Important: for the Alexandrov format, the unvaried genome
        # corresponds to the average over all mutation frequencies (which is
        # precisely 1/96 for all mutation frquencies because the 96 frequencies
        # add up to 1!!!)
        # For Shiraishi signatures the "unvaried" genome is represented by 
        # the following exemplary matrix:
        # [ 1/6  1/6  1/6  1/6  1/6  1/6    <- 16.7% for each of 6 base changes
        #   1/4  1/4  1/4  1/4   0    0     <- 25% for each flanking base
        #   ...
        #   1/2  1/2   0    0    0    0 ]   <- 50% for each transcription dir.

        if (is.vector(genomes[[g]]) & is.numeric(genomes[[g]])) {

            tss <- sum( (genomes[[g]] - mean(genomes[[g]]))^2 )

        } else if (is.matrix(genomes[[g]]) || is.data.frame(genomes[[g]])) {

            unvargenome <- genomes[[g]]
            unvargenome[1,] = rep(1/6, 6)
            for (ii in 2:nrow(genomes[[g]])) {
                unvargenome[ii,] = c(rep(1/4, 4), rep(0, 2))
            }
            if ( (nrow(genomes[[g]])%%2) == 0) {
                # even number of rows; last one must be transcription direction
                unvargenome[nrow(genomes[[g]]),] = c(rep(1/2, 2), rep(0, 4))
            }

            tss <- sum( (genomes[[g]] - unvargenome)^2 )
            
        } else { # something unknown, complain
            stop("genomes must be either of the Alexandrov (numeric vectors) or the Shiraishi format (matrices or data frames)!")
        }

        # finally, compute the explained variance
        evar <- 1 - (rss / tss)

        # for the formulas, see for example
        # https://www.rdocumentation.org/packages/SomaticSignatures/
        #                              versions/2.8.3/topics/numberSignatures

        expvar <- c(expvar, evar)
    }

    # finally, name the explained variances like the genomes
    if (!is.null(names(genomes))) {
        names(expvar) <- names(genomes)
    }

    return(expvar)
}



plotExplainedVariance <- function(genome, signatures,
                                  minExplainedVariance=NULL,
                                  minNumSignatures=2, maxNumSignatures=NULL) {
    # input: gnome = mutation counts for a single genome as matrix or vector
    #        signatures = a list of signatures (matrices or vectors)
    #        [need to have the same format]

    if (!is.list(signatures)) {
        stop("Parameter 'signatures' must be a list of signature objects!")
    }
    if (!is.data.frame(signatures[[1]]) & !is.matrix(signatures[[1]])
        & !(is.vector(signatures[[1]]) & is.numeric(signatures[[1]]))) {
        
        stop("Signatures must be data.frames, matrices or numeric vectors!")
    }

    if (is.list(genome)) { # it's a list of genomes
        if (length(genome) == 1) {  # accept if only one!
            genome = genome[[1]]
        } else { # more than one genome
            stop("plotExplainedVariance can plot the explained variance for only one genome!")
        }
    }


    if (!is.data.frame(genome) & !is.matrix(genome)
        & !(is.vector(genome) & is.numeric(genome))) {
        stop("genome must be a data.frame, matrix or numeric vector!")
    }


    # check the genome format; same as signature format?
    if (length(genome) != length(as.vector(as.matrix(signatures[[1]])))) {
        stop("Formats of genome and signatures must match!")
    }

    # is the signatures are unnamed, name them by enumerating them
    if(is.null(names(signatures))) {
        names(signatures) <- paste0("sign_",1:length(signatures))
    }

    # if maximum number of signatures is not defined, set it to the total
    # number of signatures
    if(is.null(maxNumSignatures)) {
        maxNumSignatures = length(signatures)
    }


    bestDecompositions = list()
    
    # iterating for all possible numbers of signatures from minNumSignatures
    # to maxNumSignatures

    for (k in minNumSignatures:maxNumSignatures) {
    
        bestDecompositions[[length(bestDecompositions)+1]] <-
            getBestDecomp4Ksignatures(genome, signatures, k)
    }

        
    # now plot: 
    # k on x-axis; 
    # expl. var on y-axis; 
    # min. k with expl. var >= thres in red

    k <- unlist(lapply(bestDecompositions, function(x) { x$k } ))      # vector
    explVar <- unlist(lapply(bestDecompositions, function(x) { x$explVar } ))
                                                                       # vector
    sigList <- lapply(bestDecompositions, function(x) { x$sigList } )  # list

    minExplVarIndex = NULL
    if (!is.null(minExplainedVariance)) {
        whichExceed <- which(explVar >= minExplainedVariance)
        if (length(whichExceed) > 0) {
            minExplVarIndex = whichExceed[1]  # first to exceed the threshold
        }
    }

    yrange <- c(max(min(explVar),0), 1)
    
    # basic plot
    plot(k, explVar, xlab="number of signatures",
         ylab="highest explained variance", ylim=yrange)
    abline(h=1, lty=2, col="black")
    
    # indicate threshold for minimum explained variance?
    if (!is.null(minExplVarIndex)) {
        
        # indicate desired threshold for mexplained variance
        abline(h=minExplainedVariance, lty=2, col="red")
        y4Label = minExplainedVariance - 2*(max(explVar)-min(explVar))/100
        text(c(maxNumSignatures), c(y4Label),
             paste("min. explained variance: ",minExplainedVariance),
             cex=0.6, pos=2, col="red")

        # also highlight the signatures exceeding the threshold
        points(c(k[minExplVarIndex]), c(explVar[minExplVarIndex]), bg="red",
               col="black", pch=21)

        text(c(k[minExplVarIndex]),
             c(min(explVar)+(max(explVar)-min(explVar))/2),
             paste(sigList[[minExplVarIndex]], collapse="\n"), cex=0.6,
             pos=4, col="red")
        
        abline(v=k[minExplVarIndex], lty=3, col="red")
    }
}



# re-compose predicted genome(s) from the signatures and the exposures:

composeGenomesFromExposures <- function(exposures, signatures) {

    if (class(exposures) == "numeric") {
        # this is only one genome, use a list nonetheless for later iteration
        exposures = list(exposures)
    }

    if (class(signatures) != "list") {
        stop("Parameter signatures must be a list object.")
    }

    predGenomes <- list()
    
    for (e in 1:length(exposures)) {

        if(is.null(exposures[[e]])) {
            # didn't exceed minExplainedVariance ... can't reconstruct this
            predGenomes[[e]] = NULL
            next
        }
        
        exp <- exposures[[e]]

        if (class(exp) != "numeric") {
            stop("exposures must be numeric vectors.")
        }

        
        # make sure we replace NA by 0 (signatures that were not selected in a
        # search, i.e., they of course have 0 contribution/exposure
        exp[which(is.na(exp))] = 0
        
        # check that we have as many exposure values as signatures
        if(length(exp)!=length(signatures)) {
            stop("Number of exposure values different from number of signatures.")
        }

        # multiply exposures with the ressponding signatures and sum the results
        for (ii in 1:length(exp)) {

            if (ii==1) {
                # first: set pred
                pred = exp[ii]*signatures[[ii]]
            } else {
                # add the following
                pred = pred + exp[ii]*signatures[[ii]]
            }
        }

        predGenomes[[e]] = pred
    }

    return(predGenomes)
}



# evaluate the quality of a single decomposition
# (can also be plotted)

evaluateDecompositionQuality <- function(exposure, signatures, genome,
                                         plot=FALSE) {
    if (class(genome) == "list") {
        # must be a single genome
        stop("genome must be an the mutation frequencies of an individual genome (in Alexandrov or Shiraishi format).")
    }
    if (class(exposure) != "numeric") {
        stop("exposure must be the exposure vector of a single decomposition.")
    }
    if (class(signatures) != "list") {
        stop("signatures must be a list object.")
    }

    # check the genome format; same as signature format?
    if (length(genome) != length(as.vector(as.matrix(signatures[[1]])))) {
        stop("Formats of genome and signatures must match!")
    }
    
    if (!is.logical(plot)) {
        stop("plot must be logical (TRUE or FALSE).")
    }
    

    decQual = list()
    
    # compute explained variance
    decQual$explainedVariance <- computeExplainedVariance(exposure,
                                                          signatures, genome)

    # re-compose genome from exposure and signatures
    predGenome <- composeGenomesFromExposures(exposure, signatures)[[1]]

    # compute correlation of re-composed and original genome
    decQual$pearsonCorr <- cor(as.numeric(as.matrix(predGenome)),
                               as.numeric(as.matrix(genome)), method="pearson")
    # (note: as.numeric(as.matrix(X)) works on vectors, matrices and data frames!)

    
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



# determine the differences between a target and a set of signatures
# evaluates: Frobenius distance; residual sum of squares (squared error);
# or all distance measures that can be used with function "dist" from
# package "stats"
determineSignatureDistances <- function(target, signatures,
                                        method="euclidean") {

    if (!is.list(signatures)) {
        stop("Parameter 'signatures' must be a list of signature objects!")
    }
    if (!is.data.frame(signatures[[1]]) & !is.matrix(signatures[[1]])
        & !(is.vector(signatures[[1]]) & is.numeric(signatures[[1]]))) {
        
        stop("Signatures must be data.frames, matrices or numeric vectors!")
    }
    
    # is the signatures are unnamed, name them by enumerating them
    if(is.null(names(signatures))) {
        names(signatures) <- paste0("sign_",1:length(signatures))
    }

    if (!is.data.frame(target) & !is.matrix(target)
        & !(is.vector(target) & is.numeric(target))) {
        
        stop("Target signature must be a data.frame, matrix or numeric vector!")
    }

    # allow Frobenius distance only for data.frame or matrix
    if (method == "frobenius" &
        !(is.data.frame(target) || is.matrix(target)) ) {

        stop("Frobenius distance can be used only for Shiraishi signatures (matrix or data.frame)!")
    }
    
    distances = unlist(lapply(signatures, function(s) {
        # first, check we have the same format as the target signature
        if (as.vector(length(target)) != as.vector(length(s))) {
            stop("Formats of target and the set of signatures must match!")
        }

        # determine distance
        if (method == "frobenius") {
            # Frobenius distance? use own implementation
            d = computeFrobeniusNorm(as.matrix(target)-as.matrix(s))
        } else if (method == "rss") {
            # Residual sum of squares (squared error)
            d = computeRSS(as.vector(target),as.vector(s))
        } else {
            # other distances? use "dist" from "stats"
            d = dist(rbind(as.vector(target),as.vector(s)), method=method)
        }

        d # return the computed distance
    } ))
    
    #names(distances) = names(signatures)
    #sort(distances)

    return(distances)
}



# map a set of signatures to another set of signatures, choosing for each
# signature of the first set the best match in the second set (with or
# without replacement, i.e, allowing multiple mappings to the same signatures
# or requiring unique mappings)
mapSignatureSets <- function(fromSignatures, toSignatures,
                             method="euclidean", unique=FALSE) {
    
    if (!is.list(fromSignatures)) {
        stop("Parameter 'fromSignatures' must be a list of signature objects!")
    }
    if (!is.data.frame(fromSignatures[[1]]) & !is.matrix(fromSignatures[[1]])
        & !(is.vector(fromSignatures[[1]]) & is.numeric(fromSignatures[[1]]))) {
        
        stop("fromSignatures must be data.frames, matrices or numeric vectors!")
    }

    if (!is.list(toSignatures)) {
        stop("Parameter 'toSignatures' must be a list of signature objects!")
    }
    if (!is.data.frame(toSignatures[[1]]) & !is.matrix(toSignatures[[1]])
        & !(is.vector(toSignatures[[1]]) & is.numeric(toSignatures[[1]]))) {
        
        stop("toSignatures must be data.frames, matrices or numeric vectors!")
    }

    # from and to must be of the same format
    if (class(fromSignatures[[1]]) != class(toSignatures[[1]])
        || length(fromSignatures[[1]]) != length(toSignatures[[1]])) {
        stop("fromSignatures and toSignatures must be of the same type and format!")
    }

    # if we require a unique mapping, the number of fromSignatures must not
    # exceed the number of toSignatures
    if (unique & (length(fromSignatures) > length(toSignatures))) {
        stop("for a unique mapping the number of fromSignatures must not exceed the number of toSignatures!")
    }

    # make sure we have signature names
    if (is.null(names(fromSignatures))) {
        names(fromSignatures) <- paste0("sign_", c(1:length(fromSignatures)))
    }
    if (is.null(names(toSignatures))) {
        names(toSignatures) <- paste0("sign_", c(1:length(toSignatures)))
    }

    # keep the names of the fromSignatures, but remove them from the object
    fromSigNames <- names(fromSignatures)
    names(fromSignatures) <- NULL

    
    if (!unique) {
        
        # no unique mapping is required; for each fromSignature
        # simply choose the most similar one from toSignatures
        mapping <- sapply(fromSignatures, function(from) {
            d <- determineSignatureDistances(from, toSignatures, method=method)
            names(sort(d))[1] # sort by distance, take name of first
        } )

        names(mapping) <- fromSigNames
    } else {
        
        # we want a uniqur mapping, iteratively identify the best 
        # from->to pair (shortest distance), then remove both to find the
        # next best ...
        mapping = c()
        fromNames <- fromSigNames
        
        while (length(fromSignatures) > 0) {
            dists <- sapply(fromSignatures, function(from) {
                d <- determineSignatureDistances(from, toSignatures,
                                                 method=method)
                sort(d)[1] # sort by distance, take name and distance of first
            } )

            # get from and to signatures names
            choice <- which(dists == min(dists))
            fromS <- fromNames[choice]
            toS <- names(choice)

            # keep the mapping
            mapping <- c(mapping, toS)
            names(mapping)[length(mapping)] <- fromS

            # remove the two signatures from the respective sets
            fromNames <- fromNames[-choice]
            fromSignatures[choice] <- NULL
            toSignatures[toS] <- NULL
        }

        # reorder mapping according to fromSignature input
        mapping <- mapping[fromSigNames]
    }

    return(mapping)
}


# downgrade Shiraishi signatures (to samller patterns and/or removin
# transcription direction
downgradeShiraishiSignatures <- function(signatures, numBases=NULL,
                                         removeTrDir=FALSE) {

    if (!is.list(signatures)) {
        stop("Parameter 'signatures' must be a list of signature objects!")
    }
    if (!is.data.frame(signatures[[1]]) & !is.matrix(signatures[[1]])) {
        
        stop("Signatures must be data.frames or matrices!")
    }

    if (is.null(removeTrDir)) {
        removeTrDir = FALSE
    }
    if(!is.logical(removeTrDir)) {
        stop("Value of removeTrDir must be logical!")
    }

    if (is.null(numBases) & !removeTrDir) {
        stop("At least one of numBases and removeTrDir must be specified!")
    }

    if (!is.null(numBases)) {
        if (!is.numeric(numBases)
            || (numBases!=round(numBases)) || ((numBases%%2)==0) ) {
            stop("Value of numBases must be an odd integer!")
        }
    }

    newsigs <- lapply(signatures, function(sig) {
        # first, check current format of the signatures

        # if we have an even number of lines, the last one must be for the
        # transcription direction
        haveTrDir = !(nrow(sig)%%2)

        # number of bases in the sequence pattern?
        haveBases = nrow(sig) - as.numeric(haveTrDir)

        if(!is.null(numBases) && (numBases >= haveBases)) {
            stop("Value of numBases must be smaller than the signatures' current number of bases!")
        }

        # determine which indices to remove
        removeRows <- c()

        if(!is.null(numBases)) {
            # upstream
            removeRows <- c(removeRows,
                            c(1:((haveBases-numBases)/2))+1)
            # downstream
            removeRows <- c(removeRows,
                            c(((haveBases-(haveBases-numBases)/2)+1):haveBases))
        }

        if(removeTrDir & haveTrDir) {
            removeRows <- c(removeRows, nrow(sig))
        }

        sig[-removeRows,]
    })

    return(newsigs)
}
