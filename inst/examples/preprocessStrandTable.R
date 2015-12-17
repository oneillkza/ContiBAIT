data(exampleStrandFreq)

strandStates <- preprocessStrandTable(exampleStrandFreq)

show(strandStates[[1]]) # WW-WC-CC matrix
show(strandStates[[2]]) # WW-CC only matrix
