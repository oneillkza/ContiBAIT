computeQCMatrix <- function(strandMatrix)
{
	contig.NAs <- lapply(1:nrow(strandMatrix), function(x){which(!is.na(strandMatrix[x,]))})
	
	computeOneContig <- function(contig.num)
	{
		contig.Pairs <- contig.num:nrow(strandMatrix)
		calcOnePair <- function(c1, c2)
		{
			length( intersect( contig.NAs[[c1]], contig.NAs[[c2]] ) )
			#length(which(!is.na(strandMatrix[c1,]) & !is.na(strandMatrix[c2,])))
		}
		contig.QC <- sapply(contig.Pairs, calcOnePair, contig.num)
		contig.QC
	}
	
	qc.list <- lapply(1:nrow(strandMatrix), computeOneContig)
	qc.mat <- matrix(rep(0, nrow(strandMatrix)^2),  nrow=nrow(strandMatrix), ncol=nrow(strandMatrix))
	
	for(i in 1:length(qc.list))
	{
		
		qc.mat[i,i:nrow(qc.mat)] <- qc.list[[i]]
		qc.mat[i:nrow(qc.mat),i] <- qc.list[[i]]
	}
	
	qc.mat
}