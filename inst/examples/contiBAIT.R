#Get a list of BAM files containing libraries for cells from the same organism, aligned to the same genome
#In this case these are the example BAM files provided with the package (hence the call to system.file);
data(exampleChrTable)

example.dir <- file.path(system.file(package='contiBAIT'), 'extdata')
dividedChr <- divideMyChr(exampleChrTable, splitBy=1000000)

orderedContigs <- contiBAIT(path=example.dir, filter=dividedChr, pairedEnd=FALSE) 

show(orderedContigs)
