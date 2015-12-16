#Get an example BAM file and generate a chromosome table featuring fragment names and lengths

example.bam <- list.files(file.path(system.file(package='contiBAIT'), 'extdata'), full.names=TRUE)[1]

chrTable <- makeChrTable(example.bam) 

show(chrTable)