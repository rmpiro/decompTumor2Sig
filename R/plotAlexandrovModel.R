
#####################
# internal function #
#####################

#' plotAlexandrovModel (internal function)
#'
#' `plotAlexandrovModel()` plots a single signature or the mutation
#' frequency data for a single genome of the Alexandrov-type model.
#'
#' @usage plotAlexandrovModel(mutData, numBases, trDir, colors = NULL,
#' strip = NULL)
#' @param mutData (Mandatory) The signature or genome mutation frequency data
#' to be plotted.
#' @param numBases The number of bases of the sequence pattern.
#' @param trDir Logical value specifying whether transcription strand
#' information is present.
#' @param colors Vector of six colors to be used for the base change data.
#' If \code{NULL} (default), the colors are set to those used by the COSMIC
#' website.
#' @param strip Background color for strip labels. If \code{NULL} (default),
#' "papayawhip" will be used.
#' @return Returns (or draws) a plot according to the Alexandrov model of
#' mutational signatures.
#' @author Rosario M. Piro\cr Politecnico di Milano\cr Maintainer: Rosario
#' M. Piro\cr E-Mail: <rmpiro@@gmail.com> or <rosariomichael.piro@@polimi.it>
#' @references \url{http://rmpiro.net/decompTumor2Sig/}\cr
#' Krueger, Piro (2019) decompTumor2Sig: Identification of mutational
#' signatures active in individual tumors. BMC Bioinformatics
#' 20(Suppl 4):152.\cr
#' @importFrom ggplot2 ggplot geom_bar scale_fill_manual theme_bw theme rel
#' element_text element_blank element_rect guides facet_grid scale_x_continuous
#' @keywords internal
plotAlexandrovModel <- function(mutData, numBases, trDir,
                                colors = NULL, strip = NULL) {

    # this function assumes all probabilities are sorted as follows:
    # first, alphabetically according to the six base changes; 
    # then, blockwise alphabetically according to the flanking bases (within
    # each base change block.
    
    # base change types and colors
    changes <- c("[C>A]","[C>G]","[C>T]","[T>A]","[T>C]","[T>G]")

    if (is.null(colors)) {
        # default colors are those used by the COSMIC website!
        colors = c("#1ebff0", "#050708", "#e62725", "#cbcacb",
            "#a1cf64", "#edc8c5")
    }

    if (is.null(strip)) {
        strip = "papayawhip"
    }
    
    # number of blocks (one per base change)
    numBlocks <- length(changes)

    # make sure we have the correct number of colors
    stopifnot(length(colors) == numBlocks)
    
    # numer of probabilies per block (depends on sequence pattern length)
    probsPerBlock <- 4^(numBases-1)

    # transcription strand information?
    strands <- c("transcription independent")
    if (trDir) {
        strands <- c("transcription strand", "opposite strand")
    }

    # quickly verify, we got the correct counts (should not happen)
    stopifnot(length(mutData) ==
                  numBlocks * probsPerBlock * (as.numeric(trDir)+1))

    
    # get list of all sequence patterns and set colors per blocks
    patterns <- buildSortedAlexandrovSignaturePatternList(numBases, FALSE)
    patterns <- gsub("\\[.{3}\\]", " ~ ", patterns)
    #patterns <- gsub("\\[", "", gsub(">.]", "", patterns))

    # build data frame with:
    # <strand> <base change> <seq. pattern> <pos in block> <probability/prob.>
    dat <- data.frame(strand = factor(rep(strands,
                          each = probsPerBlock * numBlocks),
                          levels = strands),
                      change = factor(rep(rep(changes,
                          each = probsPerBlock), length(strands)),
                          levels = changes),
                      pattern = rep(patterns, length(strands)),
                      pos = seq(probsPerBlock),
                      probability = as.vector(mutData)
                      )

    xlabels <- element_blank()
    if (numBases == 3) {
        # if only one base or too many bases (>3); don't print sequence labels
        xlabels <- element_text(face="bold", vjust=0.5,
                                size=rel(0.8), angle=90)
    }
    
    ggp <- ggplot(dat, aes_string(x = "pos", y = "probability",
                                  fill = "change")) +
            geom_bar(stat = "identity", position = "identity", width = 0.8) + 
            scale_fill_manual(values = colors) +
            scale_x_continuous(breaks = seq(probsPerBlock),
                               labels = patterns[seq(probsPerBlock)]) +
            theme_bw() + 
            theme(axis.text.x = xlabels,
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.y = element_text(size = rel(1.2)),
                  axis.title.y = element_text(size = rel(1.2)),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  strip.text = element_text(face = "bold", size = rel(1.2)),
                  strip.background = element_rect(fill = strip)) +
            guides(fill = FALSE) +
            facet_grid(strand ~ change)

    ggp$labels$y = "Mutation type frequency or probability"
    
    ggp

}
