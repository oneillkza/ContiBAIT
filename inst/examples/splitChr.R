#Get an example BAM file and generate a chromosome table featuring fragment names and lengths, then split into 50kb fragments

data(exampleChrTable)
chrTable <- makeChrTable(example.bam, verbose=FALSE) 
dividedChr <- divideMyChr(chrTable, splitBy=50000)
show(dividedChr)