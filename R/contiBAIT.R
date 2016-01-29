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
#' @import diagram
#' @import methods
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

  #subset data to only include data above readlimit
  strandFrequencyList[[1]][which(strandFrequencyList[[2]] < readLimit)] <- NA 
  

  if(verbose){message('-> Processing table and filtering data [2/6]')}
  strandStateMatrixList <- preprocessStrandTable(strandFrequencyList[[1]], lowQualThreshold=0.8)

  #create weighting criteria; median of read depth 
  filtWeight <- strandFrequencyList[[2]][which(rownames(strandFrequencyList[[2]]) %in% rownames(strandStateMatrixList[[1]])  ),]
  libWeight <- apply(filtWeight, 1, median)

  if(verbose){message('-> Clustering data ', cluster, 'x using ', clusNum, ' cores [3/6]')}      
  slaveNum <- makeCluster(clusNum)
  linkage.groups <- clusterContigs(strandStateMatrixList[[1]], randomWeight=libWeight, snowCluster=slaveNum, recluster=cluster, randomise=TRUE)
  stopCluster(slaveNum)
  
   # make orientation calls for each group; WW and CC only
  if(verbose){message('-> Reorienting discordant fragments [4/6]')}
  reorientedTable <- reorientLinkageGroups(linkage.groups, strandStateMatrixList[[1]])

  if(verbose){message('-> Merging related linkage groups [5/6]')}
  linkage.merged <- mergeLinkageGroups(linkage.groups, reorientedTable[[1]])
 
  reorientedTable <- reorientLinkageGroups(linkage.groups, reorientedTable[[1]], previousOrient=reorientedTable[[2]], verbose=FALSE)

  if(makePlots != TRUE)
  {
    contigOrder <- orderAllLinkageGroups(linkage.merged, reorientedTable[[1]], strandFrequencyList[[1]], strandFrequencyList[[2]])
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
    contigOrder <- orderAllLinkageGroups(linkage.merged, reorientedTable[[1]], strandFrequencyList[[1]], strandFrequencyList[[2]], saveOrdered=TRUE, randomAttempts=500)
    ideogramPlot(strandFrequencyList[[3]], strandFrequencyList[[4]], chrTable, plotBy='chr', orderFrame=contigOrder, orientationData=reorientedTable[[2]], verbose=TRUE)

    dev.off()
  }

 if(saveName != FALSE){save(strandFrequencyList,strandStateMatrixList,linkage.groups,linkage.merged, reorientedTable, contigOrder, file=paste(saveName, '_', cluster, 'clusters_linkgeData.Rd', sep=""))}

  return(contigOrder)

}
