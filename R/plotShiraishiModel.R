
#####################
# internal function #
#####################

#' plotShiraishiModel (internal function)
#'
#' `plotShiraishiModel()` plots a single signature or the mutation
#' frequency data for a single genome of the Shiraishi-type model.
#'
#' @usage plotShiraishiModel(mutData, numBases, trDir, colors = NULL)
#' @param mutData (Mandatory) The signature or genome mutation frequency data
#' to be plotted.
#' @param numBases The number of bases of the sequence pattern.
#' @param trDir Logical value specifying whether transcription strand
#' information is present.
#' @param colors Vector of four colors to be used for the base data.
#' If \code{NULL} (default), the colors are consensus colors used for
#' sequence logos.
#' @return Returns (or draws) a plot according to the Shiraishi model of
#' mutational signatures.
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
#' @importFrom ggseqlogo ggseqlogo make_col_scheme
#' @importFrom ggplot2 ggplot geom_polygon theme_bw theme scale_y_continuous
#' scale_fill_manual dup_axis labs  rel aes guides element_blank element_rect
#' element_text
#' @importFrom gridExtra grid.arrange
#' @keywords internal
plotShiraishiModel <- function(mutData, numBases, trDir, colors = NULL) {

    # To avoid R CMD complaining about "global variable"
    arrow_poly <- NULL
    
    # this function assumes the input data is a Shiraishi probability matrix
    
    # bases and colors
    bases <- c("A", "C", "G", "T")

    if (is.null(colors)) {
        # default colors are those used by ggseqlogo
        colors = c("#109648", "#255c99", "#f7b32b", "#d62839")
    }

    # other default colors
    panelBackCol <- "#eeeeee"   # background color for mutation base panel
    trColors <- c("coral3", "skyblue3")
    arrowCol <- "#a8a8a8"
    
    # set own colorschemes for bases
    colscheme <- make_col_scheme(chars=bases, cols=colors)
                                  
    # make sure we have the correct number of colors
    stopifnot(length(colors) == length(bases))
    
    # numer of probabilies per block (=base)
    probsPerBlock <- 4

    # quickly verify, we got the correct counts (should not happen)
    stopifnot(ncol(mutData) == 6 &&
                  nrow(mutData) == numBases + as.numeric(trDir))

    # construct dummy sequence data for ggseqlogo::geom_logo
    # maximum 0.1% accurracy is fair enough (i.e., 1000 sequences)
    numSeqs <- 1000

    # probabilities of mutated bases C and T:
    probMut <- c(sum(mutData[1, seq_len(3)]), sum(mutData[1, seq(from=4,to=6)]))

    
    # sequence logos for mutated bases

    baseBoth <- ggseqlogo("CT", method="probability", col_scheme=colscheme) +
        theme_bw() + 
        labs(caption="mutated base", size=(4)) +
        theme(axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.background = element_rect(fill=panelBackCol),
              plot.caption = element_text(hjust=0.5, size=rel(1.2))) +
        geom_text(data = data.frame(prob=paste(sprintf("%.1f",
                                        probMut*100),"%")),
                  aes_string(x=seq_len(2), y=-0.1, label="prob"), size=rel(4),
                  color="black") +
        geom_vline(xintercept = 1.5, linetype="dotted")
        
    # sequence logos for variant bases

    # first base
    seqs <- rep(bases[-2], roundIntegerSum(numSeqs * mutData[1, seq_len(3)],
                                           targetSum=numSeqs))
    # add second base
    seqs <- paste0(seqs,
                   rep(bases[-4],
                       roundIntegerSum(numSeqs * mutData[1, seq(from=4,
                                                                to=6)],
                                       targetSum=numSeqs)))

    mutBoth <- ggseqlogo(seqs, method="probability", col_scheme=colscheme) +
        theme_bw() + 
        labs(caption="variant base", size=(4)) +
        theme(axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.background = element_rect(fill=panelBackCol),
              plot.caption = element_text(hjust=0.5, size=rel(1.2))) +
         geom_vline(xintercept = 1.5, linetype="dotted")

    
    # sequence logos for upstream and downstream bases

    flankBases <- list()

    if (numBases > 1) {
        rows <- seq(2, numBases)
        basePos <- c(0, seq(from=-(numBases%/%2), to=-1),   # used for labeling
                     seq(from=1, to=numBases%/%2))
        basePos <- gsub("+-", "-", gsub("^", "+", as.character(basePos)),
                        fixed=TRUE)

        flankBases <- list()

        for (r in rows) {
            seqs <- rep(bases, round(numSeqs * mutData[r, seq_len(4)]))

            ggp <- ggseqlogo(seqs, method="probability", col_scheme=colscheme) +
                labs(caption=paste("pos", basePos[r])) + 
                theme_bw() +
                theme(axis.title = element_blank(),
                      axis.text.x = element_blank(),
                      axis.ticks = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.grid.major.x = element_blank(),
                      plot.caption = element_text(hjust=0.5, size=rel(1.2)))
            if (r == 2) {
                # first: ticks left
                #ggp <- ggp + theme(axis.text.y = element_blank())
                ggp <- ggp + scale_y_continuous(sec.axis = dup_axis()) +
                    theme(axis.text.x = element_blank(),
                          axis.text.y.right = element_blank())
            } else if (r == numBases) {
                # last: ticks right
                ggp <- ggp + scale_y_continuous(sec.axis = dup_axis()) +
                    theme(axis.text.x = element_blank(),
                          axis.text.y.left = element_blank())
            } else {
                # others: no ticks
                ggp <- ggp + scale_y_continuous(sec.axis = dup_axis()) +
                    theme(axis.text = element_blank())
            }
        
            flankBases[[length(flankBases)+1]] <- ggp
        }
    }

    # dummy logo, only with probability labels on right side
    labelR <- ggseqlogo(rep("C", 10), method="probability",
                        col_scheme=make_col_scheme("C", cols="#ffffff")) +
        scale_y_continuous(sec.axis = dup_axis()) +
        labs(caption=" ", size=(4)) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y.left = element_blank(),
              axis.title=element_blank(),
              plot.caption = element_text(hjust=0.5, size=rel(1.2)))

    if (trDir) {
        trPlot <- ggplot(data=data.frame(x=factor(c("+", "-"),
                                                  levels=c("+", "-")),
                             probability=mutData[nrow(mutData), seq_len(2)]),
                  aes_string(x="x", y="probability")) +
                  geom_bar(stat="identity", fill=trColors) +
                  geom_text(aes_string(label="x"), vjust=-0.2, size=rel(6)) +
                  scale_y_continuous(sec.axis = dup_axis(), limits=c(0,1)) +
                  labs(caption="transcription\nstrand", size=(4)) +
                  theme_bw() +
                  theme(panel.grid = element_blank(),
                        panel.border = element_blank(),
                        axis.ticks = element_blank(),
                        axis.title = element_blank(),
                        #axis.text.x = element_text(size=rel(1.2)),
                        axis.text.x = element_blank(),
                        axis.text.y.left = element_blank(),
                        plot.caption = element_text(hjust=0.5, size=rel(1.2)))
        
    } else {
        # empty plot
        trPlot <- ggplot() + theme(panel.background=element_blank())
    }
    
    # arrow
    xs <- c(1/3, 2/3, 2/3, 5/6, 1/2, 1/6, 1/3, 1/3)
    xs2 <- c(1/3, 2/3, 2/3, 5/6, 1/2, 1/6, 1/3, 1/3) + 0.75
    ys <- c(1/4, 1/4, 1/2, 1/2, 3/4, 1/2, 1/2, 1/4)
    vs <- rep("arrow", length(xs))

    arrow <- ggplot() +
        geom_polygon(data = arrow_poly,
                     aes(x = xs, y = ys, fill = vs)) +
        geom_polygon(data = arrow_poly,
                     aes(x = xs2, y = ys, fill = vs)) +
        theme(axis.title=element_blank(),
              axis.text=element_blank(),
              axis.line=element_blank(),
              axis.ticks=element_blank(),
              panel.background=element_blank(),
              plot.background=element_blank()) +
        guides(fill=FALSE) +
        scale_fill_manual(values=c(arrow=arrowCol))

    
    # get list of plots for grid.arrange
    if (numBases > 1) {
        plots <- c(list(labelR, mutBoth, trPlot, arrow),
                   flankBases[seq(numBases%/%2)],
                   list(baseBoth),
                   flankBases[seq(numBases%/%2+1, numBases-1)])
        nGridCols <- numBases
    } else {
        # numBases = 1
        plots <- list(labelR, mutBoth, trPlot, arrow, baseBoth)
        nGridCols <- 3   # we need at least 3 columns
    }
    
    # define layout for grid.arrange
    heights <- c(2,1,2)

    widths <- rep(1, nGridCols)
    widths[1] <- 1.25
    widths[nGridCols] <- 1.25
    widths[nGridCols%/%2+1] <- 1.8

    layout <- matrix(NA, nrow=3, ncol=nGridCols)
    layout[1, nGridCols%/%2] <- 1
    layout[1, nGridCols%/%2+1] <- 2
    layout[1, nGridCols] <- 3
    layout[2, nGridCols%/%2+1] <- 4
    if (numBases > 1) {
        layout[3, ] <- seq(from=5, to=numBases+4)
    } else {
        layout[3, nGridCols%/%2+1] <- 5
    }
    
    # call grid.arrange
    args <- list(nrow=3, ncol=nGridCols, widths=widths,
                 heights=heights, layout_matrix=layout)

    do.call("grid.arrange", c(plots, args))

}
