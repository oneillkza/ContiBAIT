#Get a list of BAM files containing libraries for cells from the same organism, aligned to the same genome
#In this case these are the example BAM files provided with the package (hence the call to system.file);

example.dir <- file.path(system.file(package='contiBAIT'), 'extdata')

orderedContigs <- contiBAIT(path=example.dir, pairedEnd=F) 

show(orderedContigs)
