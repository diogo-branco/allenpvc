# GSE71585 metadata

meta <- data.frame(
    Title="GEO accession data GSE71585 as a SingleCellExperiment",
    Description="Celular taxonomy of the primary visual cortex in adult mice
        based on single cell RNA-sequencing from a study performed by the Allen
        Institute for Brain Science. In said study 49 transcriptomic cell types
        are identified.",
    BiocVersion="3.6",
    Genome="mm10",
    SourceType=".csv.gz",
    SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71585",
    SourceVersion="Jan 04 2016",
    Species="Mus musculus",
    TaxonomyId="10090",
    Coordinate_1_based=TRUE,
    DataProvider="GEO",
    Maintainer="Diogo P. P. Branco <diogo.pp.branco@gmail.com>",
    RDataClass="SingleCellExperiment",
    DispatchClass="Rda",
    RDataPath="allenpvc/allenpvc.rda"
)

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
