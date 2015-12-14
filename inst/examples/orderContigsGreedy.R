#Get a data.frame of ordered contigs from cells from the same organism, aligned to the same genome
 
data(exampleLGList)
data(exampleWCMatrix)
data(exampleStrandFreq)
data(exampleReadCounts)

orderedGroup <- orderContigsGreedy(exampleLGList, exampleWCMatrix, exampleStrandFreq, exampleReadCounts, 1)
show(orderedGroup[[3]])

