expandMergedContigs <- function(contigsInOrder, contigMap,linkage.group)
{
	newOrder <- vector()
	for (i in contigsInOrder)
		newOrder <- append(newOrder, contigMap[[i]])
	
	return(intersect(newOrder,linkage.group))
}