data("exampleWCMatrix")
 
clusteredContigs <- clusterContigs(exampleWCMatrix, verbose=FALSE)
show(clusteredContigs)
show(clusteredContigs[[1]])

reorientedMatrix <- reorientLinkageGroups(clusteredContigs, exampleWCMatrix)
mergedLinkageGroups <- mergeLinkageGroups(clusteredContigs,reorientedMatrix[[1]])
