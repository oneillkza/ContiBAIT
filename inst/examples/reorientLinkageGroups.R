data(exampleWCMatrix)
clusteredContigs <- clusterContigs(exampleWCMatrix, randomise=FALSE)

reorientedMatrix <- reorientLinkageGroups(clusteredContigs,
										  exampleWCMatrix)

# Note that in this example data, everything is correctly oreiented to
# to begin with, so all contigs come out as + orientation
