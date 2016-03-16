#Get a list of BAM files containing libraries for cells from the same organism, aligned to the same genome
#In this case these are the example BAM files provided with the package (hence the call to system.file);
data("exampleDividedChr")
data("exampleWCMatrix")
data("exampleReadCounts")

library(BiocParallel)

example.dir <- file.path(system.file(package='contiBAIT'), 'extdata')

chrGrange <- exampleDividedChr[which(exampleDividedChr$name %in% rownames(exampleWCMatrix))]

relatedLibList <- lapply(seq_len(length(unique(seqnames(chrGrange)))), function(x) findSimilarLibraries(exampleWCMatrix, exampleReadCounts, chrGrange, x, cluster=1, clusterParam=MulticoreParam()))

show(relatedLibList)
