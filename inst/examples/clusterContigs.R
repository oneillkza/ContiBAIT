# data(exampleWWCCMatrix)
# data(exampleWCMatrix)
# 
# clusteredContigs <- clusterContigs(exampleWWCCMatrix, verbose=FALSE)
# show(clusteredContigs)
# show(clusteredContigs[[1]])
# 
# LGOrientations <- reorientLinkageGroups(clusteredContigs, exampleWWCCMatrix)
# reorientedMatrix <- reorientStrandTable(exampleWCMatrix, linkageGroups = clusteredContigs, orientation = LGOrientations)
# mergedLinkageGroups <- mergeLinkageGroups(clusteredContigs,reorientedMatrix)
