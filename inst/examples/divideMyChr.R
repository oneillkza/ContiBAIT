#Get an example BAM file and generate a chromosome table featuring fragment names and lengths, then split into 50kb fragments

data("exampleChrTable")
dividedChr <- divideMyChr(exampleChrTable, splitBy=1000000)