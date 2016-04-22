####################################################################################################
#' Plot ordering of contigs within a single linkage group.
#' @param contigOrder matrix from orderAllContigs with the subdivided linkage groups and the names of the contigs to plot
#' @param lg Integer specifying the linkage group by which to plot. Default is all
#' @param verbose prints messages to the terminal (default is TRUE)
#' @import ggplot2
#' @example inst/examples/plotContigOrder.R
#' @return A ggplot object (which will be plotted automatically if not assigned).
#' @export
####################################################################################################

plotContigOrder <- function(contigOrder, lg='all', verbose=TRUE)
{

	masterGroups <- sapply(1:nrow(contigOrder), function(x) strsplit(as.character(contigOrder[,1]), "\\.")[[x]][1])

	if(lg == 'all'){lg <- seq(1:length(unique(masterGroups)))}
	for(link in lg)
	{
		if(verbose){message(' -> Processing ', link)}
		contigOrderGrp <- contigOrder[grep(paste(unique(masterGroups)[link],"\\.", sep=""), contigOrder[,1]),]

		contigChr <- sub(':.*', '', contigOrderGrp[,2])
		primaryContigChr <- names(sort(table(contigChr), decreasing=TRUE))[1]

		contigLengths <- sub('.*:', '', contigOrderGrp[,2])
		contigStarts <- sub('-.*', '', contigLengths)

		orderedLocation <- unlist(sapply(1:length(unique(contigOrderGrp[,1])), function(x) rep(x, length(contigOrderGrp[,1][which(contigOrderGrp[,1] == unique(contigOrderGrp[,1])[x])]))))

		contigOrderFrame <- data.frame(chr=contigChr, start=as.numeric(contigStarts)/10^6, bin=orderedLocation)

		rsquareOrient <- summary(lm(bin ~ start, data=contigOrderFrame))$r.squared
		rsquare <- summary(lm(bin ~ -1+start, data=contigOrderFrame))$r.squared
		if(rsquare < rsquareOrient)
		{
			contigOrderFrame[,1:2] <- contigOrderFrame[nrow(contigOrderFrame):1, 1:2]
			rsquare <- summary(lm(bin ~ 1+start, data=contigOrderFrame))$r.squared
		}

		rsquare <- round(rsquare, digits=2)
		print(ggplot(contigOrderFrame, aes_string("bin", "start") )+
		geom_point(aes_string(x="bin", y="start" , colour="chr"), size=2)+
		labs(x="contiBAIT predicted location of contigs", y="Assembly ordered location of contigs (Mb)")+
		geom_smooth(method="lm", se=FALSE, formula=y~x-1)+
		ggtitle(paste("Plot of ", length(contigChr), " fragments (", length(unique(contigOrderFrame$bin)), " sub-linkage groups)\nR-squared = ", rsquare,  sep="")))
	}
}