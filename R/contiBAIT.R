# Copyright (c) 2015, Mark Hills
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

#    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


####################################################################################################
#' contiBAIT -- master function to process strand-seq libraries into ordered linkage groups 
#' 
#' @param path  String denoting location of Strand-seq bam files (default is ".")
#' @param cluster  Integer denoting the number of reclusterings to be performed for creating linkage groups (default is 1)
#' @param clusNum  Number of parallel processors to use when clustering contigs. Default is 1. 
#' @param saveName  String denoting the file name for saved data. If FALSE, no intermediate files are saved (default is FALSE)
#' @param filter  additional file to split chromosomes based on locations. If this parameter is blank,
#' a filter table will be automatically generated from the header of the first file in bamFileList
#' @param readQual Integer dictating the minimal mapping quality required for a read to be accepted. Default is 10. 
#' @param readLimit  Minimum number of reads on a contig to make a strand call. Default is 10
#' @param pairedEnd  Whether the bam files being read are in paired end format. Default is TRUE. Note,
#' since paired reads will be the same direction, only first mate read of pair is used in output
#' @param makePlots Logical determining whether plots should be created. Default is TRUE
#' @param verbose prints messages to the terminal (default is TRUE)
#' 
#' @return ordered contigs in bed format. Depending on options, intermediate files and plots will also be generated
#' @import snow
#' @importFrom S4Vectors DataFrame
#' @example inst/examples/contiBAIT.R
#' @export
#' @include AllClasses.R
####################################################################################################


contiBAIT <- function(path=".", 
                        cluster=1, 
                        clusNum=1, 
                        saveName=FALSE, 
                        filter=FALSE, 
                        readQual=10,
                        readLimit=10, 
                        pairedEnd=TRUE,
                        makePlots=FALSE,
                        verbose=TRUE)
{

  #Create directory to store all the files
  bamFileList <- list.files(path=path, pattern=".bam$", full.names=TRUE)

  if(verbose){message('RUNNING CONTIBAIT ON ', length(bamFileList), ' BAM FILES!\n\n PARAMETERS FOR BAM ANALYSIS: \n----------------------------------\n     -> paired end data=', pairedEnd, '\n     -> mapping quality=', readQual, '\n     -> bed filter=', if(length(filter)==1){'FALSE'}else{'TRUE'}, '\n     -> minimal reads required=', readLimit, '\n\n PARAMETERS FOR CLUSTERING: \n----------------------------------\n     -> number of reclusters=', cluster, '\n     -> number of cores to use=', clusNum, '\n\n ADDITIONAL PARAMETERS:', '\n----------------------------------\n     -> saving intermediate files=', if(saveName==FALSE){'FALSE'}else{'TRUE'}, '\n     -> creating analysis plots=', makePlots, '\n----------------------------------')}

  if(verbose){message('-> Creating read table from bam files [1/6]')}
  strandFrequencyList <- strandSeqFreqTable(bamFileList, filter=filter, qual=readQual, pairedEnd=pairedEnd, BAITtables=TRUE)

  # subset data with: strandFrequencyList[[1]][which(strandFrequencyList[[2]] < 100)] <- NA
  if(saveName != FALSE){save(strandFrequencyList, file=paste(saveName, '_table.Rd', sep="")) }
  
  #subset data to only include data above readlimit
  strandFrequencyList[[1]][which(strandFrequencyList[[2]] < readLimit)] <- NA 
  

  if(verbose){message('-> Processing table and filtering data [2/6]')}
  strandStateMatrixList <- preprocessStrandTable(strandFrequencyList[[1]], lowQualThreshold=0.8)
  if(saveName != FALSE){save(strandStateMatrixList, file=paste(saveName, '_strands.Rd', sep=""))}

  #create weighting criteria; median of read depth 
  filtWeight <- strandFrequencyList[[2]][which(rownames(strandFrequencyList[[2]]) %in% rownames(strandStateMatrixList[[1]])  ),]
  libWeight <- apply(filtWeight, 1, median)

  if(verbose){message('-> Clustering data ', cluster, 'x using ', clusNum, ' cores [3/6]')}      
  slaveNum <- makeCluster(clusNum)
  #linkage.groups <- clusterContigs(strandStateMatrixList[[1]], randomWeight=libWeight, snowCluster=slaveNum, recluster=cluster, randomise=TRUE, minimumLibraryOverlap=10, similarityCutoff=0.8)
  linkage.groups <- clusterContigs(strandStateMatrixList[[1]], randomWeight=libWeight, snowCluster=slaveNum, recluster=cluster, randomise=TRUE)

   
  stopCluster(slaveNum)
  
  if(saveName != FALSE){ save(linkage.groups, file=paste(saveName, '_LG_', cluster, 'x_reclust.Rd', sep="") ) }	
if(saveName != FALSE){ save(linkage.groups2, file=paste(saveName, '_LG_homo_', cluster, 'x_reclust.Rd', sep="") ) } 

   # make orientation calls for each group; WW and CC only
  if(verbose){message('-> Reorienting discordant fragments [4/6]')}
  reorientedTable <- reorientLinkageGroups(linkage.groups, strandStateMatrixList[[1]])

  if(saveName != FALSE){save(reorientedTable, file=paste(saveName, '_', cluster, 'x_reclust_reoriented.Rd', sep=""))}

  if(verbose){message('-> Merging related linkage groups [5/6]')}
  linkage.merged <- mergeLinkageGroups(linkage.groups, reorientedTable[[1]])
 
  reorientedTable <- reorientLinkageGroups(linkage.merged, reorientedTable[[1]], verbose=FALSE)

  if(saveName != FALSE){save(linkage.merged, file=paste(saveName, '_', cluster, 'x_reclust_merged.Rd', sep="")  )}

  if(verbose){message('-> Checking for high quality sex chromosome fragments [6/6]')}

  if(nrow(strandStateMatrixList[[2]]) > 2)
  {
    if(verbose){message(paste('  -> ', nrow(strandStateMatrixList[[2]]), ' found. Processing.', sep=""))}
    filtWeightSex <- strandFrequencyList[[2]][which(rownames(strandFrequencyList[[2]]) %in% rownames(strandStateMatrixList[[2]])  ),]
    libWeightSex <- apply(filtWeightSex, 1, median)
    # cluster the sex groups if any (should only be present in males, assuming either C or W)
    slaveNum <- makeCluster(clusNum)
    linkage.groups.sex <- clusterContigs(strandStateMatrixList[[2]], randomWeight=libWeightSex, snowCluster=slaveNum, recluster=cluster, randomise=TRUE)
    stopCluster(slaveNum)
    reorientedGroups.sex <- reorientLinkageGroups(linkage.groups.sex, strandStateMatrixList[[2]])
    linkage.merged.sex <- mergeLinkageGroups(linkage.groups.sex, reorientedGroups.sex[[1]])
    if(saveName != FALSE){save(linkage.merged.sex, file=paste(saveName, '_', cluster, 'reclust_merged_sex.Rd', sep="")  )}
  } else {
    if(verbose){message('  -> None found')}
  }

  if(makePlots != TRUE)
  {
    contigOrder <- orderAllLinkageGroups(linkage.merged, reorientedTable[[1]], strandFrequencyList[[1]], strandFrequencyList[[2]])
    if(saveName != FALSE){save(contigOrder, file=paste(saveName, '_', cluster, 'reclust_ordered_LGs.Rd', sep="")  )}
  }else{
    if(saveName == FALSE){saveName = 'contiBAIT'}

    pdf(paste(saveName, '_all_plots.pdf', sep='_'))

    plotWCdistribution(strandFrequencyList[[1]], filterThreshold=0.8)

    if(length(filter) < 3){
      chrTable <- makeChrTable(bamFileList[1])
    }else{
      chrTable <- filter
    }
    makeBoxPlot(chrTable, linkage.merged)

    plotLGDistances(linkage.merged, strandStateMatrixList[[1]])

    barplotLinkageGroupCalls(linkage.merged, chrTable)
    barplotLinkageGroupCalls(linkage.merged, chrTable, by='chr')
    contigOrder <- orderAllLinkageGroups(linkage.merged, reorientedTable[[1]], strandFrequencyList[[1]], strandFrequencyList[[2]], saveOrdered=TRUE)
    ideogramPlot(strandFrequencyList[[3]], strandFrequencyList[[4]], chrTable, plotBy='chr', orderFrame=contigOrder, verbose=TRUE)

    dev.off()
  }

  return(contigOrder)

}

 #find groups that are smaller than they should be by finding the SD of linkage group size
  #sdLinkageSize <- ceiling(sd(do.call(rbind, lapply(linkage.groups, function(x) length(x)))))
  #failedLinks <- unlist(linkage.groups[sdLinkageSize:length(linkage.groups)])

  #bob <- lapply(1:length(linkage.groups), function(x) gsub("_.", "", linkage.groups[[x]]))                


#unlist(sapply(seq(1, length(result[[2]])), function(x) rep(names(result[[2]][x]), length(result[[2]][[x]]) ) ))
#orderedGroups <- data.frame(LG=vector(), name=vector())
#runContiBAIT()

### plotting functions

#if(length(filter) < 3)
#{
#chrTable <- makeChrTable(list.files(path=path, pattern=".bam$", full.names=TRUE)[1])
#}else{
#chrTable <- data.frame(chr=as.character(filter[,4]), length=filter[,3]-filter[,2], stringsAsFactors=F)
#}

#assemblyBED <- chrTable
#assemblyBED$chr <- matrix(unlist(strsplit(chrTable$chr, ':')), ncol=2, byrow=T)[,1]
#plotLGDistances(linkage.merged, strandStateMatrixList[[1]], saveFile=paste(dataNames, '_LGdistances_', cluster, 'x_clusters', sep=''))

#makeBoxPlot(chrTable, linkage.merged, saveFile=paste(dataNames, '_included', sep=''))

# OR
#  makeBoxPlot(chrTable, linkage.merged, saveFile=paste(dataNames, '_included', sep=''), sex.contigs=linkage.merged.sex)
#assemblyBED <- chrTable
#library(stringr)
# nam <- str_split_fixed(assemblyBED$chr, "_", 3)
# nam <- nam[,1]
# assemblyBED$chr <- nam

#barplotLinkageGroupCalls(linkage.merged, assemblyBED, saveFile=paste(dataNames, '_barplot_LG', sep=''))
#barplotLinkageGroupCalls(linkage.merged, assemblyBED, by='chr', saveFile=paste(dataNames, '_barplot_chr', sep=''))

#save(strandFrequencyList, strandStateMatrixList, linkage.groups, linkage.merged, reorientedTable, filtWeight, libWeight, filtWeightSex, libWeightSex, file=paste(dataNames, '_all_data_.Rd', sep=""))

#chrTable$chr <- as.character(matrix(unlist(strsplit(chrTable$chr, ':')), ncol=2, byrow=T)[,1])



#lengthPerLG <- sapply(seq(1,length(linkage.merged)), function(x) {sum(as.numeric(matrix(unlist(strsplit(matrix(unlist(strsplit(linkage.merged[[x]], ':')), ncol=2, byrow=T)[,2], '-')), ncol=2, byrow=T)[,2])-as.numeric(matrix(unlist(strsplit( matrix(unlist(strsplit(linkage.merged[[x]], ':')), ncol=2, byrow=T)[,2], '-')), ncol=2, byrow=T)[,1]))})
#lengthPerLG <- lengthPerLG/1000000
#LGname <- c(paste('LG', 1:length(linkage.merged), sep=""))
#maxX <- max(lengthPerLG)
#barplot(lengthPerLG, names.arg=LGname, las=2, xlab='Linkage group', ylab='Length in Mb', ylim=c(0,maxX))
