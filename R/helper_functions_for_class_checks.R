
### helper functions for checking classes ###

######################
# internal functions #
######################

#' is.probability.vector (internal function)
#'
#' `is.probability.vector()` checks whether the input object is a numeric
#' vector of probabilities with a total sum of 1.
#' 
#' @param x Object to be checked.
#' @return Logical value (true or false).
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @keywords internal
is.probability.vector <- function(x) {
    if (is.vector(x) & is.numeric(x)) {
        return(isTRUE(all.equal.numeric(sum(x, na.rm=TRUE),
                                        tolerance=1e-5, 1)) &
                   all(x>=0 & x<=1, na.rm=TRUE))
    }

    return(FALSE) # not even a numeric vector
}


#' is.probability.matrix (internal function)
#'
#' `is.probability.matrix()` checks whether the input object is a numeric
#' matrix of probabilities with a total sum of 1 for every row. Each row
#' must have 6 columns (for Shiraishi format).
#' 
#' @param x Object to be checked.
#' @return Logical value (true or false).
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @keywords internal
is.probability.matrix <- function(x) {
    if(is.matrix(x) & is.numeric(x)) {
        if(ncol(x)==6) {
            return(all(apply(x, 1, is.probability.vector)))
        }
    }

    return(FALSE) # not even a numeric matrix
}


#' is.probability.data.frame (internal function)
#'
#' `is.probability.data.frame()` checks whether the input object is a numeric
#' data.frame of probabilities with a total sum of 1 for every row. Each row
#' must have 6 columns (for Shiraishi format).
#' 
#' @param x Object to be checked.
#' @return Logical value (true or false).
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @keywords internal
is.probability.data.frame <- function(x) {
    if(is.data.frame(x)) {
        if(all(apply(x, 1, is.numeric)) & ncol(x)==6) {
            return(all(apply(x, 1, function(r) {
                is.probability.vector(as.numeric(r))
            })))
        }
    }
    
    return(FALSE)  # not even a numeric data.frame
}


#' is.probability.object (internal function)
#'
#' `is.probability.object()` checks whether the input object is a numeric
#' vector, matrix of data.frame of probabilities with a total sum of 1 for
#' every row. Matrices and data.frames must have 6 columns (for Shiraishi
#' format).
#' 
#' @param x Object to be checked.
#' @return Logical value (true or false).
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @keywords internal
is.probability.object <- function(x) {
    return(is.probability.vector(x) ||
               is.probability.matrix(x) ||
                   is.probability.data.frame(x))
}



#' is.probability.vector.list (internal function)
#'
#' `is.probability.vector.list()` checks whether the input object is a list of
#' numeric vectors.
#' 
#' @param x Object to be checked.
#' @return Logical value (true or false).
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @importFrom plyr compact 
#' @keywords internal
is.probability.vector.list <- function(x) {
    if (is.list(x)) {
        x <- compact(x)  # exclude NULL elements from check (can happen
                         # for exposures when the min. explained variance
                         # was not reached!
        return(all(vapply(x, is.probability.vector, FUN.VALUE=logical(1))))
    }

    return(FALSE) # if not even a list ...
}


#' is.probability.matrix.list (internal function)
#'
#' `is.probability.matrix.list()` checks whether the input object is a list of
#' numeric matrices. Must have 6 columns (for Shiraishi format).
#' 
#' @param x Object to be checked.
#' @return Logical value (true or false).
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @importFrom plyr compact 
#' @keywords internal
is.probability.matrix.list <- function(x) {
    if (is.list(x)) {
        x <- compact(x)  # exclude NULL elements from check (can happen
                         # for exposures when the min. explained variance
                         # was not reached!
        return(all(vapply(x, is.probability.matrix, FUN.VALUE=logical(1))))
    }

    return(FALSE) # if not even a list ...
}

#' is.probability.data.frame.list (internal function)
#'
#' `is.probability.data.frame.list()` checks whether the input object is a
#' list of numeric data.frame objects. Must have 6 columns (for Shiraishi
#' format).
#' 
#' @param x Object to be checked.
#' @return Logical value (true or false).
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @importFrom plyr compact 
#' @keywords internal
is.probability.data.frame.list <- function(x) {
    if (is.list(x)) {
        x <- compact(x)  # exclude NULL elements from check (can happen
                         # for exposures when the min. explained variance
                         # was not reached!
        return(all(vapply(x, is.probability.data.frame, FUN.VALUE=logical(1))))
    }

    return(FALSE) # if not even a list ...
}


######################
# external functions #
######################

#' isShiraishiSet
#'
#' `isShiraishiSet()` checks whether the input object is a set (list) of
#' numeric objects compatible with the Shiraishi format (matrices or
#' data.frames of probabilities; 6 columns, each row sums up to 1). NOTE:
#' These can also be genomes compatible with the Shiraishi format!
#' 
#' @param x Object to be checked.
#' @return Logical value (true or false).
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @seealso \code{\link{decompTumor2Sig}}\cr
#' \code{\link{isSignatureSet}}\cr
#' \code{\link{readShiraishiSignatures}}
#' @examples
#' 
#' ### get Alexandrov signatures from COSMIC
#' signAlexandrov <- readAlexandrovSignatures()
#' 
#' ### convert them to the Shiraishi model
#' signShiraishi <- convertAlexandrov2Shiraishi(signAlexandrov)
#'
#' isShiraishiSet(signShiraishi)
#' @export isShiraishiSet
isShiraishiSet <- function(x) {
    return(is.probability.matrix.list(x) || is.probability.data.frame.list(x))
}


#' isAlexandrovSet
#'
#' `isAlexandrovSet()` checks whether the input object is a set (list) of
#' numeric objects compatible with the Alexandrov format (probability vectors;
#' sum up to 1). NOTE: These can also be genomes compatible with the
#' Alexandrov format!
#' 
#' @param x Object to be checked.
#' @return Logical value (true or false).
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @seealso \code{\link{decompTumor2Sig}}\cr
#' \code{\link{readAlexandrovSignatures}}\cr
#' \code{\link{isSignatureSet}}
#' @examples
#' 
#' ### get Alexandrov signatures from COSMIC
#' signAlexandrov <- readAlexandrovSignatures()
#' 
#' isAlexandrovSet(signAlexandrov)
#' @export isAlexandrovSet
isAlexandrovSet <- function(x) {
    return(is.probability.vector.list(x))
}


#' isSignatureSet
#'
#' `isSignatureSet()` checks whether the input object is a set (list) of
#' numeric objects compatible with either the Alexandrov format (probability
#' vectors; see \code{isAlexandrovSet}) or the Shiraishi format (matrices or
#' data.frames of probabilities; see \code{isShiraishiSet}).
#' NOTE: These can also be genomes compatible with one of the two formats!
#' 
#' @param x Object to be checked.
#' @return Logical value (true or false).
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @seealso \code{\link{decompTumor2Sig}}\cr
#' \code{\link{isAlexandrovSet}}\cr
#' \code{\link{isShiraishiSet}}
#' @examples
#' 
#' ### get Alexandrov signatures from COSMIC
#' signAlexandrov <- readAlexandrovSignatures()
#' 
#' isSignatureSet(signAlexandrov)
#' @export isSignatureSet
isSignatureSet <- function(x) {
    return(isShiraishiSet(x) || isAlexandrovSet(x))
}


#' isExposureSet
#'
#' `isExposureSet()` checks whether the input object is a set (list) of
#' numeric objects compatible with exposure output obtained from
#' \code{decomposeTumorGenomes}.
#' 
#' @param x Object to be checked.
#' @return Logical value (true or false).
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
#' ### load reference genome
#' refGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' 
#' ### read breast cancer genomes from Nik-Zainal et al (PMID: 22608084) 
#' gfile <- system.file("extdata",
#'          "Nik-Zainal_PMID_22608084-VCF-convertedfromMPF.vcf.gz", 
#'          package="decompTumor2Sig")
#' genomes <- readGenomesFromVCF(gfile, numBases=3, type="Alexandrov",
#'          trDir=FALSE, refGenome=refGenome, verbose=FALSE)
#' 
#' ### compute exposures
#' exposures <- decomposeTumorGenomes(genomes, signatures, verbose=FALSE)
#' 
#' isExposureSet(exposures)
#' @export isExposureSet
isExposureSet <- function(x) {
    return(is.probability.vector.list(x))
}


#' sameSignatureFormat
#'
#' `sameSignatureFormat()` checks whether two input object are sets (lists) of
#' numeric objects both compatible with the same signature format (probability
#' vectors for Alexandrov signatures and probability matrices or data.frames
#' for Shiraishi signatures). For Shiraishi signatures also the number of
#' flanking bases and the presence of transcription-strand information are
#' compared. For Alexandrov signatures also the number of triplet changes are
#' compared.
#' 
#' @param x First object to be checked.
#' @param y Second object to be checked.
#' @return Logical value (true or false).
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @seealso \code{\link{decompTumor2Sig}}\cr
#' \code{\link{isAlexandrovSet}}\cr
#' \code{\link{isShiraishiSet}}
#' @examples
#' 
#' ### get Alexandrov signatures from COSMIC
#' signAlexandrov <- readAlexandrovSignatures()
#' 
#' ### convert them to the Shiraishi model
#' signShiraishi <- convertAlexandrov2Shiraishi(signAlexandrov)
#'
#' sameSignatureFormat(signAlexandrov, signShiraishi)
#' @export sameSignatureFormat
sameSignatureFormat <- function(x, y) {
    
    if (isAlexandrovSet(x) & isAlexandrovSet(y)) {
        if (length(x[[1]]) == length(y[[1]])) {
            return(TRUE)
        }
    } else if (isShiraishiSet(x) & isShiraishiSet(y)) {
        if (nrow(x[[1]]) == nrow(y[[1]])) {
            return(TRUE)
        }
    }
    
    return(FALSE)
}

