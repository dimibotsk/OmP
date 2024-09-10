library(GenomicFeatures)
library(openxlsx)

homePath <- "/data/fleming/fousteri-te"
gtfFile <- file.path(homePath,"teanalysis","hg19_rmsk_TE.gtf")

# Create TxDb object from GTF
txdb <- suppressWarnings(makeTxDbFromGFF(gtfFile))
save(txdb,file=file.path(homePath,"teanalysis","hg19_rmsk_TE.rda"))

# Use GenomicFeatures package facilities to extract elements of interest
g <- genes(txdb,single.strand.genes.only=FALSE)
G <- unlist(g)
G$name <- names(G)
G <- as.data.frame(unname(G))
G <- G[,c(1:3,6,5)]
names(G) <- c("chromosome","start","end","name","strand")
write.table(G,file=file.path(homePath,"teanalysis","hg19_TE_coords.txt"),
    sep="\t",quote=FALSE,row.names=FALSE)
