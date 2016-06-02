#Get a list of BAM files containing libraries for cells from the same organism, aligned to the same genome
#In this case these are the example BAM files provided with the package (hence the call to system.file);
data("exampleLibList")


library(BiocParallel)

example.dir <- file.path(system.file(package='contiBAIT'), 'extdata')

exampleRange <- thoroughBed(example.dir, exampleLibList)

show(exampleRange)