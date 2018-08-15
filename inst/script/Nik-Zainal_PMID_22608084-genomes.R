# produces the files:
# Nik-Zainal_PMID_22608084-genomes-Alexandrov_3bases.Rdata
# Nik-Zainal_PMID_22608084-genomes-Shiraishi_5bases_trDir.Rdata

library(decompTumor2Sig)

refGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

gfile <- "inst/extdata/Nik-Zainal_PMID_22608084-MPF.txt.gz"


# read according to Alexandrov model (3 bases, no transcription direction)
genomes <- readGenomesFromMPF(gfile, numBases=3, type="Alexandrov",
                              trDir=FALSE, refGenome=refGenome,
                              verbose=FALSE)

save(genomes, file="Nik-Zainal_PMID_22608084-genomes-Alexandrov_3bases.Rdata")

# read accoring to Shiraishi model (5 bases, transcript direction)
transcriptAnno <-
    TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene

genomes <- readGenomesFromMPF(gfile, numBases=5, type="Shiraishi",
                              trDir=TRUE, refGenome=refGenome,
                              transcriptAnno=transcriptAnno, verbose=FALSE)

save(genomes, file="Nik-Zainal_PMID_22608084-genomes-Shiraishi_5bases_trDir.Rdata")

