# produces the files:
# Nik-Zainal_PMID_22608084-pmsignature-G.Rdata
# Nik-Zainal_PMID_22608084-pmsignature-Param.Rdata
# Nik-Zainal_PMID_22608084-pmsignature-sig1.tsv
# Nik-Zainal_PMID_22608084-pmsignature-sig2.tsv
# Nik-Zainal_PMID_22608084-pmsignature-sig3.tsv
# Nik-Zainal_PMID_22608084-pmsignature-sig4.tsv
# Nik-Zainal_PMID_22608084-pmsignature-sigall.Rdata

### ATTENTION: every run will be slightly different, due to stochastic 
### processes in pmsignature; see below!

library(pmsignature)
library(decompTumor2Sig)


NUMSIG <- 4

numBases <- 5
trDir <- TRUE
sigType <- "independent"  # Shiraishi


# get reference genome and transcript annotation
refGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
transcriptAnno <-
    TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene

# read mutation data by pmsignature
G <- readMPFile(infile="inst/extdata/Nik-Zainal_PMID_22608084-MPF.txt.gz", trDir=trDir, numBases=numBases, type=sigType, bs_genome=refGenome, txdb_transcript=transcriptAnno)

# save the G object as Rdata file
save(G, file="Nik-Zainal_PMID_22608084-pmsignature-G.Rdata")

# get signatures (and exposures) by de novo signature inference (pmsignature)
### IMPORTANT: this is a stochastic process, so signatures can slightly change
### every time this is executed; above all: there is no guarantee that they
### will have the same order as before!
Param <- getPMSignature(G, K=NUMSIG)

# save the G and Param file as Rdata file
save(G, Param, file="Nik-Zainal_PMID_22608084-pmsignature-Param.Rdata")


# convert the signatures to the decompTumor2Sig format
signatures <- getSignatureListFromEstimatedParameters(Param)


# save the signatures in individual TSV flat files
for (ii in seq_along(signatures)) {
    write.table(signatures[[ii]], file=paste0("Nik-Zainal_PMID_22608084-pmsignature-sig",ii,".tsv"), sep="\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
}


# now read them by decompTumor2Sig 
# (signatures will then be named by the file name)
sigfiles <- dir(pattern="Nik-Zainal_PMID_22608084-pmsignature-sig")

signatures <- loadShiraishiSignatures(files=sigfiles)

# save them again as Rdata object
save(signatures, file="Nik-Zainal_PMID_22608084-pmsignature-sigall.Rdata")
