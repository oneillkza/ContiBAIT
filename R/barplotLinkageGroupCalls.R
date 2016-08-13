barplotLinkageGroupCalls.func <- function(object, chrTable, by='lg', bySize=TRUE, returnTable=NULL, percentage=NULL)
{

	#Calculate length of each chromosome represented for one linkage group:
	calcOneGroupChr <- function(lg.num, linkage.chr, linkage.lengths, complete.list)
	{
		lg.chr <- linkage.chr[[lg.num]]
		lg.lengths <- linkage.lengths[[lg.num]]
		chr.vector <- rep(0, length(complete.list))
		names(chr.vector) <- complete.list
		chr.represented <- unique(lg.chr)
		chr.lengths <- sapply(chr.represented, 
							  function(chr.name){sum(lg.lengths[which(lg.chr==chr.name)])})
		chr.vector[chr.represented] <- chr.lengths
		chr.vector
	}

	reCompileObject <- function(linkage.lengths)
	{
		totLength <- lapply(linkage.lengths, sum)
		newOrderNames <- names(totLength[order(unlist(totLength), decreasing=TRUE)])
		object <- LinkageGroupList(object[newOrderNames],
								  names= newOrderNames
								  )
	}

	if(bySize){
		linkage.lengths <- lapply(object, 
								  function(x){ width(chrTable[chrTable$name %in% x])  }) 
		object <- reCompileObject(linkage.lengths)
	}
	linkage.lengths <- lapply(object, 
							  function(x){ width(chrTable[chrTable$name %in% x])  }) 
	linkage.chr <- lapply(object, 
						  function(x){as.character(seqnames(chrTable)[chrTable$name %in% x])  })

	complete.list <- unique(unlist(linkage.chr))

	chr.table <- sapply(1:length(linkage.chr), 
						calcOneGroupChr, linkage.chr, linkage.lengths, complete.list)
	chr.table2 <- matrix(unlist(chr.table), nrow=nrow(chr.table))
	rownames(chr.table2) <- rownames(chr.table)
	colnames(chr.table2) <- paste('LG', 1:ncol(chr.table), sep="")
	chr.table2 <- chr.table2[order(rownames(chr.table2)),]
	chr.table <- chr.table2 / 10^6

	chr.table <- chr.table[mixedsort(rownames(chr.table)),]

	if(!(is.null(percentage)))
	{
		if(by=='lg')
		{
			perLength <- unlist(lapply(linkage.lengths, sum))/10^6
			names(perLength) <- c(paste('LG', 1:length(perLength), sep=""))
			perFrame <- data.frame(sapply(seq_len(length(perLength)), function(x) round(chr.table[,which(colnames(chr.table) == names(perLength)[x])]/perLength[x]*100, digits=1)))
			perFrame <- apply(perFrame, 2, sort, decreasing=TRUE)
			perFrame <- perFrame[which(apply(perFrame, 1, 
											 function(x){length(which(x == 0))} != nrow(perFrame) )),]
			rownames(perFrame) <- paste('ch', 1:(nrow(perFrame)), sep="")
			colnames(perFrame) <- paste('LG', 1:ncol(chr.table), sep="")
		}else{
			perLength <- sapply(complete.list, function(x) sum(width(chrTable[which(seqnames(chrTable) == x)]))/10^6)
			names(perLength) <- complete.list
			perFrame <- data.frame(t(sapply(seq_len(length(perLength)), function(x) round(chr.table[which(rownames(chr.table) == names(perLength)[x]),]/perLength[x]*100, digits=1)))) 
			perFrame <- t(apply(perFrame, 1, sort, decreasing=TRUE))
			perFrame <- perFrame[,which(apply(perFrame, 2, 
											 function(x){length(which(x == 0))} != nrow(perFrame) ))]
			perFrame <- cbind(perFrame, apply(perFrame, 1, function(x) 100-sum(x)))
			colnames(perFrame) <- c(paste(1:(ncol(perFrame)-1), sep=""), "empty")
			rownames(perFrame) <- complete.list
			perFrame <- as.matrix(perFrame[mixedsort(rownames(perFrame)),])
		}
		chromoFrame <- melt(perFrame)
	}else{
		chromoFrame <- melt(chr.table)
	}
	colnames(chromoFrame) <- c('chr', 'LG', 'count')
	chromoFrame$LG <- factor(chromoFrame$LG, levels=mixedsort(levels(chromoFrame$LG)) )
	chromoFrame$chr <- factor(chromoFrame$chr, levels=mixedsort(levels(chromoFrame$chr)))


	#Plot by linkage group
	if(by=='lg')
	{	
		if(!(is.null(percentage)))
		{
			jColors <- c("#000000", rep("#989898", nrow(perFrame)-1))
			names(jColors) <- levels(chromoFrame$chr)
			print(ggplot(chromoFrame, aes_string("LG", "count"))+
			geom_bar(stat="identity", aes_string(fill="chr"), colour='grey80')+
			scale_fill_manual(values=jColors)+		
			theme(axis.text.x = element_text(angle = 90, hjust = 1))+
			labs(x="Chromosome", y="DNA Represented in Chromosome (%)")+
			theme(legend.position='none')+
			ggtitle(paste("Barplot of ", length(unique(chromoFrame$LG)), 
						  " contigs clustered into ", 
						  length(unique(chromoFrame$chr)), 
						  " chromosomes",  
						  sep="")))
		}else{
			if(length(unique(chromoFrame$chr)) > 50){leg='none'}else{leg='right'}

			print(ggplot(chromoFrame, aes_string("LG", "count"))+
			geom_bar(stat="identity", aes_string(fill="chr"), colour='black')+
			theme(axis.text.x = element_text(angle = 90, hjust = 1))+
			labs(x="Linkage Group", y="DNA Represented in Linkage Groups (Mb)")+
			theme(legend.position=leg)+
			ggtitle(paste("Barplot of ", length(unique(chromoFrame$chr)), 
						  " contigs clustered into ", 
						  length(unique(chromoFrame$LG)), 
						  " linkage groups",  
						  sep="")))
		}
	}
	
	#Alternately, plot by chromosome:
	if(by=='chr')
	{
		if(!(is.null(percentage)))
		{
			jColors <- c("#000000", rep("#989898", ncol(perFrame)-2), "#ffffff")
			names(jColors) <- levels(chromoFrame$LG)
			print(ggplot(chromoFrame, aes_string("chr", "count"))+
			geom_bar(stat="identity", aes_string(fill="LG"), colour='grey80')+
			scale_fill_manual(values=jColors)+		
			theme(axis.text.x = element_text(angle = 90, hjust = 1))+
			labs(x="Chromosome", y="DNA Represented in Chromosome (%)")+
			theme(legend.position='none')+
			ggtitle(paste("Barplot of ", length(unique(chromoFrame$LG)), 
						  " linkage groups clustering into ", 
						  length(unique(chromoFrame$chr)), " chromosomes",  sep="")))
		}else{
			if(length(unique(chromoFrame$LG)) > 50){leg='none'}else{leg='right'}

			print(ggplot(chromoFrame, aes_string("chr", "count"))+
			geom_bar(stat="identity", aes_string(fill="LG"), colour='black')+
			theme(axis.text.x = element_text(angle = 90, hjust = 1))+
			labs(x="Chromosome", y="DNA Represented in Chromosome (Mb)")+
			theme(legend.position=leg)+
			ggtitle(paste("Barplot of ", length(unique(chromoFrame$LG)), 
						  " linkage groups clustering into ", 
						  length(unique(chromoFrame$chr)), " chromosomes",  sep="")))
		}
	}

	if(!(is.null(returnTable))){return(chr.table)}
}


####################################################################################################
#' Bar plot all linkage groups, with the true chromosomes of contigs coloured.
#' @param object LinkageGroupList, as generated by clusterContigs
#' @param chrTable GRanges object containing assembly information about the contigs, including a meta column called 'name' that has names matching the object. 
#' Note that the rownames of chrTable should be the contig names, as they are used in object, and the first column (chromosome name) will used to order by chromosome if 'chr' option used in by parameter. 
#' To use a bam file header, the product of makeChrTable(bamFile) is suitable for input
#' @param by whether to plot by linkage group (if 'lg') or chromosomes ('chr')
#' @param bySize logical value to return barplot either with LGs sorted by number of contigs or size (in Mb). Default is TRUE.
#' @param returnTable Logical to return chromosome length matrix. Default is NULL
#' @param percentage Logical that returns the percentage of different chromosomes or LG within the barplot. Default is NULL
#' Note to include legend, use legend=rownames(chr.table) for by='lg', and legend=colnames(chr.table) for by='chr'
#' 
#' @return a matrix of lengths of each chromosome (rows) in each linkage group (columns)
#' 
#' @aliases barplotLinkageGroupCalls barplotLinkageGroupCalls,LinkageGroupList,LinkageGroupList-method,ChrTable,ChrTable-method
#' @example inst/examples/barplotLinkageGroupCalls.R
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt 
#' @importFrom gtools mixedsort
####################################################################################################

setMethod('barplotLinkageGroupCalls',
          signature = signature(object = 'LinkageGroupList', chrTable='ChrTable'),
          definition = barplotLinkageGroupCalls.func
)
