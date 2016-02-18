data(exampleWCMatrix)
clusteredContigs <- clusterContigs(exampleWCMatrix, randomise=FALSE)

reorientedMatrix <- reorientLinkageGroups(clusteredContigs,
										  exampleWCMatrix)

exampleLGList <- mergeLinkageGroups(clusteredContigs,
									reorientedMatrix[[1]])
