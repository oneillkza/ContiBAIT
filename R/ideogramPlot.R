ideogramPlot.func <- function(WatsonFreqList,
 							  CrickFreqList, 
 							  chrTable, 
 							  plotBy='lib', 
 							  showPage=FALSE, 
 							  orderFrame=FALSE, 
 							  verbose=TRUE)
{

	capOffPlots <- function(object, capper)
	{
		object$WatsonPlot[which(object$WatsonPlot > capper)] <- capper
		object$CrickPlot[which(object$CrickPlot < -capper)] <- -capper

		pointFrameW <- data.frame(bin=1, chr=object$chr[1], max=0)
		pointFrameC <- pointFrameW
		if(length(which(object$WatsonPlot == capper)) > 0) {pointFrameW <- data.frame(bin=object$bin[which(object$WatsonPlot == capper)], 
																					  chr=object$chr[which(object$WatsonPlot == capper)], 
																					  max=capper, 
																					  lib=object$lib[which(object$WatsonPlot == capper)])}
		if(length(which(object$CrickPlot == capper)) > 0) {pointFrameC <- data.frame(bin=object$bin[which(object$CrickPlot == -capper)], 
																					 chr=object$chr[which(object$CrickPlot == -capper)], 
																					 max=-capper, 
																					 lib=object$lib[which(object$CrickPlot == capper)])}
		return(list(object, pointFrameW, pointFrameC))
	}

	WatsonFreqList <- WatsonFreqList[which(rownames(WatsonFreqList) %in% rownames(chrTable)),,drop=FALSE]
	CrickFreqList <- CrickFreqList[which(rownames(CrickFreqList) %in% rownames(chrTable)),, drop=FALSE]

	if(length(orderFrame) != 1)
	{
		WatsonFreqList <- WatsonFreqList[orderFrame[,2],,drop=FALSE]
		CrickFreqList <- CrickFreqList[orderFrame[,2],,drop=FALSE]
		chrTable <- chrTable[orderFrame$contig,]
		chrTable$chr <- sub("\\..*", "", orderFrame$LG)

	}

	if(plotBy == 'chr')
	{
		allLibraryDataFrame <- data.frame(WatsonPlot=vector(), 
							   CrickPlot=vector(), 
							   bin=vector(), 
							   chr=vector(), 
							   lib=vector())
	}

	for(lib in seq(1,ncol(WatsonFreqList)))
	{
		if(verbose){message('-> Generating plotting data [', lib, '/', ncol(WatsonFreqList), ']' )}

		WFreqs <- as.data.frame(WatsonFreqList[,lib, drop=FALSE])
		CFreqs <- as.data.frame(CrickFreqList[,lib, drop=FALSE]*-1)
		binNums <- table(chrTable$chr)
		binNums <- binNums[which(binNums >0)]
		findMax <- max(binNums)
		maxCap <- vector()
		totalReads <- 0
		allChrDataFrame <- data.frame(WatsonPlot=vector(), 
							 CrickPlot=vector(), 
							 bin=vector(), 
							 chr=vector(), 
							 lib=vector())

		for(i in unique(chrTable$chr))
		{
			WatsonPlot <- WFreqs[rownames(chrTable[which(chrTable$chr == i),]),]
			plotOffset <- findMax-length(WatsonPlot)
			WatsonPlot <- c(rep(0, plotOffset), WatsonPlot)

			CrickPlot <- CFreqs[rownames(chrTable[which(chrTable$chr == i),]),]
			plotOffset <- findMax-length(CrickPlot)
			CrickPlot <- c(rep(0, plotOffset), CrickPlot)

			chrDataFrame <- data.frame(WatsonPlot, 
						     CrickPlot, 
						     bin=seq(1, length(WatsonPlot)), 
						     chr=i, 
						     lib=colnames(WatsonFreqList)[lib])
			allChrDataFrame <- rbind(allChrDataFrame, chrDataFrame)

			binCounts <-abs(CrickPlot)+ WatsonPlot
			readsPerChr <- sum(binCounts)
			capOff <- mean(binCounts[which(binCounts > 0)])+sd(binCounts[which(binCounts > 0)])
			maxCap <- c(maxCap, capOff)
			totalReads <- totalReads+readsPerChr
		}
		ideos <- data.frame(a=-1, b=1, c=findMax-binNums, d=findMax, chr=unique(chrTable$chr))

		if(plotBy == 'chr')
		{
			allLibraryDataFrame <- rbind(allLibraryDataFrame, allChrDataFrame)

		}else{

			maxCap <- max(maxCap, na.rm=TRUE)
			plotList <- capOffPlots(allChrDataFrame, maxCap)
			levels(plotList[[1]]$chr) <- mixedsort(levels(plotList[[1]]$chr))
			levels(ideos) <- levels(plotList[[1]]$chr)
			print(ggplot()+
			geom_ribbon(data=plotList[[1]], aes(x=bin, ymin=0, ymax=WatsonPlot), fill='paleturquoise4')+
			geom_ribbon(data=plotList[[1]], aes(x=bin, ymin=0, ymax=CrickPlot), fill='sandybrown')+
			geom_point(data=plotList[[2]], aes(bin, max), colour='paleturquoise4', size=0.5)+
			geom_point(data=plotList[[3]], aes(bin, max), colour='sandybrown', size=0.5)+
			geom_rect(data=ideos, mapping=aes(ymin=a, ymax=b, xmin=c, xmax=d), fill='grey70')+

			coord_flip()+
			scale_x_reverse()+
			ggtitle(paste("Library ", colnames(WatsonFreqList)[lib],  sep=""))+
			facet_wrap( ~ chr, nrow=2, switch='x')+
		
			theme(legend.position="none", 
				  axis.ticks.y = element_blank(), 
				  axis.text.y=element_blank(), 
				  axis.title.y = element_blank(), 
				  axis.text.x=element_blank(), 
				  axis.title.x=element_blank(), 
				  axis.ticks.x = element_blank()))
		}
	}

	if(plotBy == 'chr')
	{
		for(chr in unique(allLibraryDataFrame$chr))
		{
			subsetChr <- allLibraryDataFrame[allLibraryDataFrame$chr == chr,]

			for(page in seq(1, ceiling(length(unique(subsetChr$lib))/25)))
			{
				if(showPage != FALSE){page <- showPage}
				elementStart <-	seq(0,(page*30),by=30)[page]+1
				elementEnd <- seq(0,(page*30), by=30)[page+1]
				subsetLib <- subsetChr[which(subsetChr$lib %in% unique(subsetChr$lib)[elementStart:elementEnd]),]

				counter=1
				while(length(unique(subsetLib$lib)) < 30)
				{
					levels(subsetLib$lib) <- c(levels(subsetLib$lib), counter)
					subsetLib <- rbind(subsetLib, data.frame(WatsonPlot=0, CrickPlot=0, bin=0, chr=1, lib=counter))
					counter=counter+1
				}

				binCounts <-abs(subsetLib$CrickPlot)+subsetLib$WatsonPlot
				maxCap <- mean(binCounts[which(binCounts > 0)])+sd(binCounts[which(binCounts > 0)])

				plotList <- capOffPlots(subsetLib, maxCap)
				print(ggplot()+
				geom_ribbon(data=plotList[[1]], aes(x=bin, ymin=0, ymax=WatsonPlot), fill='paleturquoise4')+
				geom_ribbon(data=plotList[[1]], aes(x=bin, ymin=0, ymax=CrickPlot), fill='sandybrown')+
				geom_point(data=plotList[[2]], aes(bin, max), colour='paleturquoise4', size=0.5)+
				geom_point(data=plotList[[3]], aes(bin, max), colour='sandybrown', size=0.5)+
				#geom_rect(data=bob, mapping=aes(ymin=a, ymax=b, xmin=c, xmax=d), fill='grey70')+

				coord_flip()+
				scale_x_reverse()+
				ggtitle(paste("Chromosome ", chr, "  (Page", page, ")", sep=""))+
				facet_wrap( ~ lib, nrow=2, switch='x')+
			
				theme(legend.position="none", 
					  axis.ticks.y = element_blank(), 
					  axis.text.y=element_blank(), 
					  axis.title.y = element_blank(), 
					  axis.text.x=element_blank(), 
					  axis.title.x=element_blank(), 
					  axis.ticks.x = element_blank()))
				if(page == showPage){break}
			}	
		}
	}
}

# Copyright (c) 2016, Mark Hills
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

#    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


####################################################################################################
#' ideogramPlot -- plots BAIT-like ideograms 
#' 
#' @param WatsonFreqList data.frame of Watson calls. Product of strandSeqFreqTable[[3]] when BAITtables=TRUE
#' @param CrickFreqList data.frame of Crick calls. Product of strandSeqFreqTable[[4]] when BAITtables=TRUE
#' @param chrTable  A data.frame consisting of chromosomes and lengths. Generated by makeChrTable(). Note rownames equal to chromosome names are required
#' @param plotBy Whether to generate a plot for each library ('lib') or a plot for each chromosome ('chr')
#' @param showPage Integer specifying which page to plot if plotBy='chr' selected. Useful when not plotting to a file. Default is FALSE
#' @param orderFrame ordered data.frame of contigs (produced by orderAllLinkageGroups). Default is FALSE, where plots will be made from elements in chrTable.
#' @param verbose prints messages to the terminal (default is TRUE)
#' 
#' @return ordered contigs in bed format. Depending on options, intermediate files and plots will also be generated
#' @import ggplot2
#' @importFrom gtools mixedsort chr
#' @importFrom S4Vectors DataFrame
#' @aliases ideogramPlot ideogramPlot,StrandReadMatrix,StrandReadMatrix-method, ChrTable, ChrTable-method 
#' @export
#' @example inst/examples/ideogramPlot.R
#' @include AllClasses.R
####################################################################################################

setMethod('ideogramPlot',
          signature = signature(WatsonFreqList='StrandReadMatrix', CrickFreqList='StrandReadMatrix', chrTable='ChrTable'),
          definition = ideogramPlot.func
)
