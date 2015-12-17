#Get a data.frame of ordered contigs from cells from the same organism, aligned to the same genome
 
data(exampleLGList)
data(exampleWCMatrix)
data(exampleStrandFreq)
data(exampleReadCounts)
contigOrder <- orderAllLinkageGroups(exampleLGList, exampleWCMatrix, exampleStrandFreq, exampleReadCounts)

show(contigOrder)

