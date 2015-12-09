bam.files <- sapply(dir('extdata'), function(single.file){system.file('extdata', single.file, package='contiBAIT')})
strand.freq <- strandSeqFreqTable(bam.files, pairedEnd = FALSE)

show(strand.freq[[1]])
show(strand.freq[[2]])
