#Get a list of BAM files containing libraries for cells from the same organism, aligned to the same genome
#In this case these are the example BAM files provided with the package (hence the call to system.file);

example.dir <- file.path(system.file(package='contiBAIT'), 'extdata')
bam.files <- dir(example.dir, full.names=TRUE)

strand.freq <- strandSeqFreqTable(bam.files, pairedEnd = FALSE)

show(strand.freq[[1]])
show(strand.freq[[2]])

