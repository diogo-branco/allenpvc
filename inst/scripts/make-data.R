
library(GEOquery)
library(SingleCellExperiment)

asc_num <- "GSE71585"
suppl <- getGEOSuppFiles(asc_num)
for (f in dir(path=asc_num, pattern="*\\.gz")) {
  gunzip(file.path(asc_num, f), overwrite = TRUE)
}

counts <- read.csv(file.path(asc_num, "GSE71585_RefSeq_counts.csv"), row.names = 1)
rpkm <- read.csv(file.path(asc_num, "GSE71585_RefSeq_RPKM.csv"), row.names = 1)
tpm <- read.csv(file.path(asc_num, "GSE71585_RefSeq_TPM.csv"), row.names = 1)
clusters <- read.csv(file.path(asc_num, "GSE71585_Clustering_Results.csv"), row.names = 1)
# Fix a minor bad formatting on some cell namings
rownames(clusters) <- sub("-", ".", rownames(clusters))

# Add spike-in data
ercc_tomato_counts <- read.csv(file.path(asc_num, "GSE71585_ERCC_and_tdTomato_counts.csv"), row.names = 1)
ercc_tomato_rpkm <- read.csv(file.path(asc_num, "GSE71585_ERCC_and_tdTomato_RPKM.csv"), row.names = 1)
# Spike-in not available in tpm, create dummy
ercc_tomato_tpm <- apply(ercc_tomato_counts, c(1, 2), function(x) NA)

counts <- rbind(counts, ercc_tomato_counts)
rpkm <- rbind(rpkm, ercc_tomato_rpkm)
tpm <- rbind(tpm, ercc_tomato_tpm)
rm(ercc_tomato_counts, ercc_tomato_rpkm, ercc_tomato_tpm)

counts <- as.matrix(counts)
rpkm <- as.matrix(rpkm)
tpm <- as.matrix(tpm)
clusters <- as.matrix(clusters)

# Order colData properly, according to the assay matrix
cell_names <- factor(rownames(clusters), levels = colnames(tpm))
clusters_ordered <- clusters[order(cell_names),]

allenpvc <- SingleCellExperiment(
  assays = list(counts = counts, tpm = tpm, rpkm = rpkm),
  colData = clusters_ordered)

# Mark spike-ins
isSpike(allenpvc, "ERCC") <- grepl("^ERCC-", rownames(allenpvc))
isSpike(allenpvc, "tdTomato") <- grepl("^tdTomato", rownames(allenpvc))

save(allenpvc, file = "allenpvc.rda")
