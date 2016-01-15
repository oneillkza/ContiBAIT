data("exampleWatsonFreq")
data("exampleCrickFreq")
data('exampleDividedChr')

singleWatsonLibrary <- new('StrandReadMatrix', exampleWatsonFreq[,2, drop=FALSE])
singleCrickLibrary <- new('StrandReadMatrix', exampleCrickFreq[,2, drop=FALSE]) 

ideogramPlot(singleWatsonLibrary, singleCrickLibrary, exampleDividedChr)