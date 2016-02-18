

plotWCdistribution.func <- function(object, filterThreshold=0.8)
{
	object <- object[, colSums(!(is.na(object))) > 0]
	breaks <-  round( seq(-1,1, by=2/40),digits=2)
	matrixHistogram <- sapply(seq(1, ncol(object)), function(x) 
		hist(object[,x], breaks=breaks, plot=FALSE)$count)
	matrixHistogram <- t(matrixHistogram)
	avContNum <- mean(apply(matrixHistogram, 1, sum))
	matrixTable <- table(c(rbinom(100000, 40, 0.5), 1:40)) /(100000/(avContNum/2)) 
	matrixTable[1] <- avContNum/4
	matrixTable[length(matrixTable)] <- avContNum/4
	colnames(matrixHistogram) <- breaks[-20]
	colors <- rep('white', 40)
	colors[which(abs(as.numeric(colnames(matrixHistogram))) >= filterThreshold)] <- rgb(1,0,0,0.5)
	nameLabs <- c(-1, rep("", 18), 0, rep("", 18), 1 )
	nameLabs[round((filterThreshold+((1-filterThreshold)/2))*40, digits=0)] <- filterThreshold
	nameLabs[round((1-filterThreshold-((1-filterThreshold)/2)) *40, digits=0) ] <- -filterThreshold

	boxplot(matrixHistogram, 
			col=colors, 
			xlab="WC frequency per contig (W-C)/(W+C)", 
			ylab="Number of contigs", 
			main=paste("WC distributions from", nrow(matrixHistogram), 
					   "libraries\nwith an average of", round(avContNum, digits=0), 
					   "fragments", sep=" ") )
	lines(names(matrixTable), matrixTable, col='mediumblue', lwd=2)

	#find the average height of WW and CC calls to compare to modelled data
	AvHomoHeight <- 
		apply(matrixHistogram[,which(abs(as.numeric(colnames(matrixHistogram))) >= filterThreshold)], 2, median)
	WWLoc <- AvHomoHeight[which(as.numeric(names(AvHomoHeight)) >= filterThreshold)]
	CCLoc <- AvHomoHeight[which(as.numeric(names(AvHomoHeight)) <= filterThreshold)]
	lines(c(1,length(WWLoc)), c(sum(CCLoc),sum(CCLoc)), col='green4', lwd=2)
	lines(c(41-length(WWLoc), 40), c(sum(WWLoc),sum(WWLoc)), col='green4', lwd=2)
	text(5, sum(CCLoc), pos=3, 'Av. CC contigs', cex=0.7, col='green4')
	text(36, sum(WWLoc), pos=3, 'Av. WW contigs', cex=0.7, col='green4')

}

####################################################################################################
#' Creates median distribution boxplots across all libraries and contigs
#' @param object object of class StrandFreqMatrix (product of strandSeqFreqTable)
#' @param filterThreshold numeric value used in assessing the threshold for homozygous strand calls. Default is 0.8.
#' @aliases plotWCdistribution plotWCdistribution,StrandFreqMatrix,StrandFreqMatrix-method 
#' @example inst/examples/plotWCdistribution.R
#' @return nothing, just plots.
#' @export
####################################################################################################

setMethod('plotWCdistribution',
          signature = signature(object = 'StrandFreqMatrix'),
          definition = plotWCdistribution.func
)
