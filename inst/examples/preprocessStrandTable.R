data(exampleStrandFreq)

strandStates <- preprocessStrandTable(exampleStrandFreq, lowQualThreshold=0.8)

show(strandStates[[1]]) # WW-WC-CC matrix
show(strandStates[[2]]) # WW-CC only matrix
