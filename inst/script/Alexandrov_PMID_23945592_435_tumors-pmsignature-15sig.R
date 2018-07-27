# produces the file:
# Alexandrov_PMID_23945592_435_tumors-pmsignature-15sig.Rdata

###
### This is run by Alexandrov_PMID_23945592_435_tumors-pmsignature-15sig.sh !!!
###

library(pmsignature)
library(decompTumor2Sig)

NUMSIG <- 15

numBases <- 5
trDir <- TRUE
sigType <- "independent"   # Shiraishi

# get reference genome and transcript annotation
refGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
transcriptAnno <-
    TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene


# read mutation data by pmsignature

gfile <- "Somatic_SNVs_Alexandrov_genomes_only-noMT-min100snvs.tsv.gz"

G <- readMPFile(infile=gfile, trDir=trDir, numBases=numBases, type=sigType, bs_genome=refGenome, txdb_transcript=transcriptAnno)


# get signatures (and exposures) by de novo signature inference (pmsignature)
### IMPORTANT: this is a stochastic process, so signatures can slightly change
### every time this is executed; above all: there is no guarantee that they
### will have the same order as before!
Param <- getPMSignature(G, K=NUMSIG)


# convert the signatures to the decompTumor2Sig format
signatures <- getSignatureListFromEstimatedParameters(Param)

save(signatures, file="Alexandrov_PMID_23945592_435_tumors-pmsignature-15sig.Rdata")

