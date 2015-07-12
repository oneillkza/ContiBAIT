computeBarPlotMatrix <- function(linkageGroups, assemblyBED)
{
	#if assemblyBED already has a length column, do nothing. Otherwise treat as 3 column bed file (chr, start, end)
	if(!"length" %in% colnames(assemblyBED))
	{
		assemblyBED$length <- apply(assemblyBED, 1, function(x){as.numeric(x[3]) - as.numeric(x[2])})
	}

	linkage.lengths <- lapply(linkageGroups, function(x){assemblyBED[x, 'length']}) 
	linkage.chr <- lapply(linkageGroups, function(x){as.character(assemblyBED[x, colnames(assemblyBED[1])])})
	
	complete.list <- unique(unlist(linkage.chr))
	
	#Calculate length of each chromosome represented for one linkage group:
	calcOneGroupChr <- function(lg.num, linkage.chr, linkage.lengths, complete.list)
	{
		lg.chr <- linkage.chr[[lg.num]]
		lg.lengths <- linkage.lengths[[lg.num]]
		chr.vector <- rep(0, length(complete.list))
		names(chr.vector) <- complete.list
		chr.represented <- unique(lg.chr)
		chr.lengths <- sapply(chr.represented, function(chr.name){sum(lg.lengths[which(lg.chr==chr.name)])})
		chr.vector[chr.represented] <- chr.lengths
		chr.vector
	}
	
	chr.table <- sapply(1:length(linkage.chr), calcOneGroupChr, linkage.chr, linkage.lengths, complete.list)
	chr.table2 <- matrix(unlist(chr.table), nrow=nrow(chr.table))
	rownames(chr.table2) <- rownames(chr.table)
	colnames(chr.table2) <- c(paste('LG', 1:ncol(chr.table), sep=""))
	chr.table2 <- chr.table2[order(rownames(chr.table2)),]
	chr.table2 <- chr.table2 / 10^6
	chr.table2
}