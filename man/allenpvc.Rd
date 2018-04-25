\name{allenpvc}

\alias{allenpvc}

\title{Adult mouse cortical cell taxonomy by single cell transcriptomics}

\description{
Celular taxonomy of the primary visual cortex in adult mice based on single
cell RNA-sequencing from a study performed by the Allen Institute for Brain
Science. In said study 49 transcriptomic cell types are identified.
This data set is the supplementary data from GEO accession
\href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71585}{GSE71585}
encapsulated in a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
}

\details{
See the vignette for examples of using these data in differential gene
expression analysis.

browseVignettes("allenpvc")

Details of how these data were creates are in the scripts/ directory of the
source package.
}

\references{
Tasic, Bosiljka, et al. Adult mouse cortical cell taxonomy revealed by single
cell transcriptomics. Nature neuroscience 19.2 (2016): 335.
}

\examples{
library(ExperimentHub)
hub <- ExperimentHub()
x <- query(hub, "allenpvc")
x
}
