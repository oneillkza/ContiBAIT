#Function to plot WC proportion of each contig vs contig representation across libraries
# strandTable should be the full WC-WW-CC strand table

plotWCvsQuality <- function(strandTable, ...)
{
	contigs.wc <- apply(strandTable, 1, function(x){length(which(x==2)) / length(which(!is.na(x)))})
	contigs.libs <- apply(strandTable, 1, function(x){length(which(!is.na(x)))})
	plot(contigs.libs, contigs.wc, pch=16, cex=0.7, col=densCols(contigs.libs, contigs.wc, colramp=colorRampPalette(rainbow(50))), ...)
}