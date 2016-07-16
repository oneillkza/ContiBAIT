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
#' @param clusterParam  Number of parallel processors to use when clustering contigs. Default is NULL. 
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
#' @import diagram
#' @importFrom S4Vectors DataFrame
#' @example inst/examples/contiBAIT.R
#' @export
#' @include AllClasses.R
####################################################################################################


contiBAIT <- function(path=".", 
                        cluster=1, 
                        clusterParam=NULL, 
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

 	this.message <- paste('RUNNING CONTIBAIT ON ', 
                      length(bamFileList), 
                      ' BAM FILES!\n\n PARAMETERS FOR BAM ANALYSIS: \n----------------------------------\n     -> paired end data=', 
                      pairedEnd, 
                      '\n     -> mapping quality=', 
                      readQual, 
                      '\n     -> bed filter=', 
                      if(length(filter)==1){'FALSE'}else{'TRUE'}, 
                      '\n     -> minimal reads required=', 
                      readLimit, 
                      '\n\n PARAMETERS FOR CLUSTERING: \n----------------------------------\n     -> number of reclusters=',
                       cluster)
  	
  if(!is.null(clusterParam))
  {

  	this.message <- paste(this.message,
    					  '\n     -> number of cores to use=', 
  						  clusterParam$workers)
  }
    this.message <- paste(this.message,
                   '\n\n ADDITIONAL PARAMETERS:', 
                   '\n----------------------------------\n     -> saving intermediate files=', 
                   if(saveName==FALSE){'FALSE'}else{'TRUE'}, 
                   '\n     -> creating analysis plots=', 
                    makePlots, 
                    '\n----------------------------------')

  if(verbose){message(this.message)}

  if(verbose){message('-> Creating read table from bam files [1/7]')}

  if(length(bamFileList) ==0)
  {
    warning('\n####################\nWARNING! NO BAM FILES PRESENT IN ', path, '. CHECK PATH LOCATION! \n####################')
    break
  }
  strandFrequencyList <- strandSeqFreqTable(bamFileList, 
                                            filter=filter, 
                                            qual=readQual, 
                                            pairedEnd=pairedEnd, 
                                            BAITtables=TRUE)

  #subset data to only include data above readlimit
  strandFrequencyList[[1]][which(strandFrequencyList[[2]] < readLimit)] <- NA 
  

  if(verbose){message('-> Processing table and filtering data [2/7]')}
  strandStateMatrixList <- preprocessStrandTable(strandFrequencyList[[1]], 
                                                 lowQualThreshold=0.8)

  #create weighting criteria; median of read depth 
  filtWeight <- strandFrequencyList[[2]][which(rownames(strandFrequencyList[[2]]) %in% rownames(strandStateMatrixList[[1]])  ),]
  libWeight <- apply(filtWeight, 1, median)

  if(verbose){message('-> Clustering data ', cluster, 'x using ',
  					if(is.null(clusterParam)){1} else {clusterParam$workers},
  					   ' cores [3/7]')}      

    linkage.groups <- clusterContigs(strandStateMatrixList[[1]], 
                                     randomWeight=libWeight, 
                                     clusterParam=clusterParam, 
                                     recluster=cluster, 
                                     randomise=TRUE)

   # make orientation calls for each group; WW and CC only
  if(verbose){message('-> Reorienting discordant fragments [4/7]')}
  reorientedTable <- reorientAndMergeLGs(linkage.groups, 
                                          cluster=cluster,
                                          clusterParam=clusterParam,
                                          strandStateMatrixList[[1]])

  linkage.merged <- reorientedTable[[3]]

   if(makePlots != TRUE)
  {
    if(verbose){message('-> Ordering contigs [7/7]')}

    contigOrder <- orderAllLinkageGroups(linkage.merged, 
                                         reorientedTable[[1]], 
                                         strandFrequencyList[[1]], 
                                         strandFrequencyList[[2]])
  }else{
    if(saveName == FALSE){saveName = 'contiBAIT'}
    if(verbose){message('-> Creating plots and ordering contigs [7/7]')}

    pdf(paste(saveName, '_all_plots.pdf', sep='_'))

    plotWCdistribution(strandFrequencyList[[1]], 
                       filterThreshold=0.8)

    if(length(filter) < 3){
      chrTable <- makeChrTable(bamFileList[1])
    }else{
      chrTable <- filter
    }
    makeBoxPlot(chrTable, 
                linkage.merged)

    plotLGDistances(linkage.merged, 
                    strandStateMatrixList[[1]])

    barplotLinkageGroupCalls(linkage.merged, 
                             chrTable)

    barplotLinkageGroupCalls(linkage.merged, 
                             chrTable, 
                             by='chr')

    contigOrder <- orderAllLinkageGroups(linkage.merged, 
                                         reorientedTable[[1]], 
                                         strandFrequencyList[[1]], 
                                         strandFrequencyList[[2]], 
                                         saveOrdered=TRUE, 
                                         randomAttempts=100)
    ideogramPlot(strandFrequencyList[[3]], 
                 strandFrequencyList[[4]], 
                 chrTable, 
                 plotBy='chr', 
                 orderFrame=contigOrder[[1]], 
                 orientationData=reorientedTable[[2]], 
                 verbose=TRUE)

    dev.off()
  }

 if(saveName != FALSE){save(this.message,
                            filter, 
                            strandFrequencyList,
                            strandStateMatrixList,
                            linkage.groups, 
                            reorientedTable, 
                            contigOrder, 
                            file=paste(saveName, '_', cluster, 'clusters_linkgeData.Rd', sep=""))}
#sapply(seq_len(length(lgs)), function(x) sum(width(filter[which(filter$name %in% lgs[[x]])])))
  return(contigOrder)

}