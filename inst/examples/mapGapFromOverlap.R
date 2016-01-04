#Find intersection between gaps in the assembly and strand state changes
#' @import GenomicRanges
#' @import rtracklayer
#import the gap file for chr1-4 from UCSC
data("exampleChrTable")
gapFile <- import.bed("http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=465319523_SLOtFPExny48YZFaXBh4sSTzuMcA&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_gap&hgta_ctDesc=table+browser+query+on+gap&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbDownBases=200&hgta_doGetBed=get+BED", genome="hg38")
sceFile <- GRanges(rep('chr4',4), IRanges(c(1410000, 1415000, 1420000, 1425000), c(1430000, 1435000, 1430000, 1435000)))

overlappingFragments <- mapGapFromOverlap(sceFile,  gapFile, exampleChrTable, overlapNum=4)
show(overlappingFragments)