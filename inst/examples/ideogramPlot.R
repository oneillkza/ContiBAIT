data("exampleWatsonFreq")
data("exampleCrickFreq")
data('exampleDividedChr')

singleWatsonLibrary <- StrandReadMatrix(exampleWatsonFreq[,2, drop=FALSE])
singleCrickLibrary <- StrandReadMatrix(exampleCrickFreq[,2, drop=FALSE])

ideogramPlot(singleWatsonLibrary, singleCrickLibrary, exampleDividedChr)
