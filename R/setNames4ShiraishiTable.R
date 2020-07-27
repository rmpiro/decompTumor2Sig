
#####################
# internal function #
#####################

#' setNames4ShiraishiTable (internal function)
#'
#' Set row and column names of a Shiraishi-model genome or signature table.
#'
#' @usage setNames4ShiraishiTable(table)
#' @param table The table to be named.
#' @return The same table with named columns and rows.
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @keywords internal
setNames4ShiraishiTable <- function(table) {

    # the table needs to be compatible with the Shiraishi model
    stopifnot(ncol(table) == 6)
    
    # determine whether we have transcription information and how many
    # bases the sequence patterns have [note: we don't use the function
    # 'determineTypeNumBasesAndTrDir' here, because that works only
    # for tables which respect the probability constraints (row sums are 1).
    # But we want to be able to use 'setNames4ShiraishiTable' also
    # for initialized tables which contain only zeroes or dummy values

    trDir <- TRUE
    if (nrow(table) %% 2) {
        # odd number of rows: only sequence pattern, no transcript direction
        trDir <- FALSE
    }

    numBases <- nrow(table) - as.numeric(trDir)


    # set column names
    colnames(table) <- c("[C>A]","[C>G]","[C>T]","[T>A]","[T>C]","[T>G]")

    # set row names
    rnames <- c("mut")

    if (numBases > 1) {
        rnames <- c(rnames, seq(-(numBases%/%2),-1), seq_len(numBases%/%2))
        rnames <- gsub("^([0-9]*)$", "+\\1", rnames)
    }
        
    if(trDir) {
        rnames <- c(rnames, "tr")
    }

    rownames(table) <- rnames

    # return the named table
    table
}
