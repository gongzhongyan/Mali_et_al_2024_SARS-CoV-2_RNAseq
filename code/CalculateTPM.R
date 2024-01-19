# script to convert to TPM
library(ggplot2)
library(rtracklayer)

GTFfile = "/Mus_musculus.GRCm39.107.gtf" # gtf file acquired from Ensembl
GTF <- import.gff(GTFfile, format="gtf", feature.type="exon")
grl <- GenomicRanges::reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))
elementMetadata(reducedGTF)$widths <- width(reducedGTF)
calc_length <- function(x) {
  sum(elementMetadata(x)$widths)
}
f1output <- data.frame(sapply(split(reducedGTF, 
                                    elementMetadata(reducedGTF)$gene_id), 
                              calc_length))
colnames(f1output) <- c("Length")
write.csv(f1output,"MusMusculusGeneLength.csv",row.names = TRUE, 
          col.names = TRUE)

# read the MusMusculus Gene Length file since it is already created.
f1output <- read.csv("file1output/MusMusculusGeneLength.csv", row.names = 1)

# read the RNA counts matrix
filepath<-"/RNA/Matrix/metaData/directory"
expMatrix_infectedr <- read.csv(paste0(filepath,
                                       "/expMatrix_infectedMouseTG_WT.csv"),
                                row.names = 1, header = TRUE, 
                                stringsAsFactors=FALSE)
f1_dds_len <- merge(expMatrix_infectedr,f1output,by=0)
rownames(f1_dds_len) <- f1_dds_len$Row.names
f1_dds_len <- f1_dds_len[-1]

#calculate tpm
f1_dds_len2 <- f1_dds_len / f1_dds_len$Length
#remove last column
f1_dds_len2 <- f1_dds_len2[-length(f1_dds_len2)]
lib_size_scale=data.frame(row.names = rownames(f1_dds_len2))
rownames(lib_size_scale) <- rownames(f1_dds_len2)
for (i in 1:length(colnames(f1_dds_len2))) {
  lib_size_scale[i] <- sum(f1_dds_len2[i])/1e6 # transcripts per million
}
colnames(lib_size_scale) <- colnames(f1_dds_len2)

f1_TPM=data.frame(row.names = rownames(f1_dds_len2))
rownames(f1_TPM) <- rownames(f1_TPM)
rownames(lib_size_scale) <- rownames(f1_dds_len2)
for (i in 1:length(colnames(f1_dds_len2))){
  f1_TPM[i] <- f1_dds_len2[i] / lib_size_scale[i]
}
colnames(f1_TPM) <- colnames(f1_TPM)
write.csv(f1_TPM,"MusMusculus_INFECTED_tpm.csv",row.names = TRUE, 
          col.names = TRUE)

