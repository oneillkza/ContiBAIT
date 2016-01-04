data("exampleWCMatrix")
 
clusteredContigs <- clusterContigs(exampleWCMatrix, verbose=FALSE)
show(clusteredContigs)
show(clusteredContigs[[1]])

LGOrientations <- reorientLinkageGroups(clusteredContigs, exampleWCMatrix)
reorientedMatrix <- reorientStrandTable(exampleWCMatrix, linkageGroups = clusteredContigs, orientation = LGOrientations)
mergedLinkageGroups <- mergeLinkageGroups(clusteredContigs,reorientedMatrix)
