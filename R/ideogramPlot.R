ideogramPlot.func <- function(WatsonFreqList,
                              CrickFreqList,
                              chrTable,
                              plotBy = 'lib',
                              sizeProportional = NULL,
                              showPage = NULL,
                              orderFrame = NULL,
                              orientationData = NULL,
                              verbose = TRUE)
{
  capOffPlots <- function(object, capper)
  {
    object$WatsonPlot[which(object$WatsonPlot > capper)] <- capper
    object$CrickPlot[which(object$CrickPlot < -capper)] <- -capper
    
    pointFrameW <- data.frame(bin = 1,
                              chr = object$chr[1],
                              max = 0)
    pointFrameC <- pointFrameW
    if (length(which(object$WatsonPlot == capper)) > 0) {
      pointFrameW <-
        data.frame(
          bin = object$bin[which(object$WatsonPlot == capper)],
          chr = object$chr[which(object$WatsonPlot == capper)],
          max = capper,
          lib = object$lib[which(object$WatsonPlot == capper)]
        )
    }
    if (length(which(object$CrickPlot == capper)) > 0) {
      pointFrameC <-
        data.frame(
          bin = object$bin[which(object$CrickPlot == -capper)],
          chr = object$chr[which(object$CrickPlot == -capper)],
          max = -capper,
          lib = object$lib[which(object$CrickPlot == capper)]
        )
    }
    return(list(object, pointFrameW, pointFrameC))
  }
  
  roorientBAITtables <-
    function(WatsonFreqList,
             CrickFreqList,
             orientationFrame)
    {
      toFlip <- orientationFrame[which(orientationFrame[, 2] == '-'), 1]
      tempWatson <- WatsonFreqList
      WatsonFreqList[which(rownames(WatsonFreqList) %in% toFlip), ] <-
        CrickFreqList[which(rownames(CrickFreqList) %in% toFlip), ]
      CrickFreqList[which(rownames(CrickFreqList) %in% toFlip), ] <-
        tempWatson[which(rownames(tempWatson) %in% toFlip), ]
      return(list(WatsonFreqList, CrickFreqList))
    }
  
  
  
  if (!(is.null(orientationData)))
  {
    flippedBAITtables <-
      roorientBAITtables(WatsonFreqList, CrickFreqList, orientationData)
    WatsonFreqList <- flippedBAITtables[[1]]
    CrickFreqList <- flippedBAITtables[[2]]
  }
  
  if (!is.null(orderFrame))
  {
    WatsonFreqList <- WatsonFreqList[orderFrame[, 2], , drop = FALSE]
    CrickFreqList <- CrickFreqList[orderFrame[, 2], , drop = FALSE]
    
    #Create new GRange in order
    chrTable <- as.data.frame(chrTable)
    rownames(chrTable) <- chrTable$name
    chrTable <- chrTable[orderFrame[, 2], ]
    chrTable <- GRanges(
      sub("\\..*", "", orderFrame[, 1]),
      IRanges(start = chrTable$start, end = chrTable$end),
      name = chrTable$name
    )
  }
  
  #Subset list so only analyzing pages required
  if (!(is.null(showPage)) && plotBy == 'lib')
  {
    WFreqSub <- WatsonFreqList[, showPage, drop = FALSE]
    CFreqSub <- CrickFreqList[, showPage, drop = FALSE]
  } else if (!(is.null(showPage)) && plotBy == 'chr') {
    chrTable <-
      chrTable[which(seqnames(chrTable) == unique(seqnames(chrTable))[showPage])]
    WFreqSub <- WatsonFreqList[chrTable$name, , drop = FALSE]
    CFreqSub <- CrickFreqList[chrTable$name, , drop = FALSE]
  } else{
    WFreqSub <- WatsonFreqList
    CFreqSub <- CrickFreqList
  }
  
  if (!(is.null(sizeProportional)))
  {
    namesToSplit <- rownames(WFreqSub)
    widths <- sub(".*:", '', namesToSplit)
    widths <- data.table(w = widths)
    widths <-
      widths[, c("w", "w2") := tstrsplit(w, "-")][, diff := as.numeric(w2) -
                                                    as.numeric(w)]$diff
    widthBins <-  ceiling(widths / sizeProportional)
    WFreqSub <- WFreqSub / widthBins
    CFreqSub <- CFreqSub / widthBins
    
    WFreqSub <-
      WFreqSub[rep(seq_len(nrow(WFreqSub)), widthBins), , drop = FALSE]
    CFreqSub <-
      CFreqSub[rep(seq_len(nrow(CFreqSub)), widthBins), , drop = FALSE]
    OldNames <- rownames(WFreqSub)
    rownames(WFreqSub) <- make.unique(rownames(WFreqSub))
    rownames(CFreqSub) <- make.unique(rownames(CFreqSub))
  }
  allLibraryDataFrameList <- list()
  
  
  if (!(is.null(sizeProportional)))
  {
    binNums <-
      sapply(unique(as.character(droplevels(
        seqnames(chrTable)
      ))), function(x)
        sum(ceiling(width(chrTable[seqnames(chrTable) == x]) / sizeProportional)))
    names(binNums) <- unique(droplevels(seqnames(chrTable)))
  } else{
    binNums <- table(droplevels(seqnames(chrTable)))
  }
  #find longest chromosome. use only first element in instance where two chromosomes are the same size
  findMax <- max(binNums)[1]
  
  for (lib in seq_len(ncol(WFreqSub)))
  {
    if (verbose) {
      message('-> Generating plotting data [',
              lib,
              '/',
              ncol(WFreqSub),
              ']')
    }
    
    WFreqs <- as.data.frame(WFreqSub[, lib, drop = FALSE])
    CFreqs <- as.data.frame(CFreqSub[, lib, drop = FALSE] * -1)
    
    
    maxCap <- vector()
    
    allChrDataFrameList <- list()
    for (i in unique(seqnames(chrTable)))
    {
      if (!(is.null(sizeProportional)))
      {
        WFreqs$oldName <- OldNames
        CFreqs$oldName <- OldNames
        WatsonPlot <-
          WFreqs[which(WFreqs$oldName %in% chrTable$name[which(seqnames(chrTable) == i)]), 1]
        CrickPlot <-
          CFreqs[which(CFreqs$oldName %in% chrTable$name[which(seqnames(chrTable) == i)]), 1]
      } else{
        WatsonPlot <- WFreqs[chrTable$name[which(seqnames(chrTable) == i)], ]
        CrickPlot <-
          CFreqs[chrTable$name[which(seqnames(chrTable) == i)], ]
      }
      plotOffset <- findMax - length(WatsonPlot)
      WatsonPlot <- c(rep(0, plotOffset), WatsonPlot)
      
      plotOffset <- findMax - length(CrickPlot)
      CrickPlot <- c(rep(0, plotOffset), CrickPlot)
      
      chrDataFrame <- data.frame(
        WatsonPlot,
        CrickPlot,
        bin = seq(1, length(WatsonPlot)),
        chr = i,
        lib = colnames(WFreqs)[1]
      )
      
      allChrDataFrameList <-
        c(allChrDataFrameList, list(chrDataFrame))
      
      binCounts <- abs(CrickPlot) + WatsonPlot
      readsPerChr <- sum(binCounts)
      capOff <-
        mean(binCounts[which(binCounts > 0)]) + sd(binCounts[which(binCounts > 0)])
      maxCap <- c(maxCap, capOff)
      
    }
    
    
    
    if (plotBy == 'chr')
    {
      allLibraryDataFrameList <-
        c(allLibraryDataFrameList, allChrDataFrameList)
      
    }	else{
      ideos <- data.frame(
        a = -1,
        b = 1,
        c = as.vector(findMax - binNums),
        d = as.vector(findMax),
        chr = unique(seqnames(chrTable))
      )
      allChrDataFrame <-
        as.data.frame(rbindlist(allChrDataFrameList, use.names = T))
      maxCap <- max(maxCap, na.rm = TRUE)
      plotList <- capOffPlots(allChrDataFrame, maxCap)
      levels(plotList[[1]]$chr) <-
        mixedsort(levels(plotList[[1]]$chr))
      levels(ideos) <- levels(plotList[[1]]$chr)
      print(
        ggplot() +
          geom_ribbon(
            data = plotList[[1]],
            aes_string(
              x = "bin",
              ymin = 0,
              ymax = "WatsonPlot"
            ),
            fill = 'paleturquoise4'
          ) +
          geom_ribbon(
            data = plotList[[1]],
            aes_string(
              x = "bin",
              ymin = 0,
              ymax = "CrickPlot"
            ),
            fill = 'sandybrown'
          ) +
          geom_point(
            data = plotList[[2]],
            aes_string("bin", "max"),
            colour = 'paleturquoise4',
            size = 0.5
          ) +
          geom_point(
            data = plotList[[3]],
            aes_string("bin", "max"),
            colour = 'sandybrown',
            size = 0.5
          ) +
          geom_rect(
            data = ideos,
            mapping = aes_string(
              ymin = "a",
              ymax = "b",
              xmin = "c",
              xmax = "d"
            ),
            fill = 'grey70'
          ) +
          
          coord_flip() +
          scale_x_reverse() +
          ggtitle(paste("Library ", colnames(WFreqs)[lib],  sep = "")) +
          facet_wrap(~ chr, nrow = 2, switch = 'x') +
          
          theme(
            legend.position = "none",
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank()
          )
      )
    }
  }
  
  
  
  if (plotBy == 'chr')
  {
    allLibraryDataFrame <-
      as.data.frame(rbindlist(allLibraryDataFrameList, use.names = T))
    for (chr in unique(allLibraryDataFrame$chr))
    {
      print(chr)
      subsetChr <-
        allLibraryDataFrame[allLibraryDataFrame$chr == chr, ]
      
      for (page in seq_len(ceiling(length(unique(
        subsetChr$lib
      )) / 26)))
      {
        elementStart <-	seq(0, (page * 26), by = 26)[page] + 1
        elementEnd <- seq(0, (page * 26), by = 26)[page + 1]
        subsetLib <-
          subsetChr[which(subsetChr$lib %in% unique(subsetChr$lib)[elementStart:elementEnd]), ]
        
        #If less than 26 chr, create empty spaces
        padIdeos <- 26 - length(unique(subsetLib$lib))
        if (padIdeos > 0)
        {
          levels(subsetLib$lib) <- c(levels(subsetLib$lib), seq(1, padIdeos))
          padData <-
            data.frame(
              WatsonPlot = rep(0, padIdeos),
              CrickPlot = rep(0, padIdeos),
              bin = rep(1, padIdeos),
              chr = rep('chr', padIdeos),
              lib = seq(1, padIdeos)
            )
          subsetLib <- rbind(subsetLib, padData)
        }
        
        binCounts <- abs(subsetLib$CrickPlot) + subsetLib$WatsonPlot
        maxCap <-
          mean(binCounts[which(binCounts > 0)]) + sd(binCounts[which(binCounts > 0)])
        
        plotList <- capOffPlots(subsetLib, maxCap)
        print(
          ggplot() +
            geom_ribbon(
              data = plotList[[1]],
              aes_string(
                x = "bin",
                ymin = 0,
                ymax = "WatsonPlot"
              ),
              fill = 'paleturquoise4'
            ) +
            geom_ribbon(
              data = plotList[[1]],
              aes_string(
                x = "bin",
                ymin = 0,
                ymax = "CrickPlot"
              ),
              fill = 'sandybrown'
            ) +
            geom_point(
              data = plotList[[2]],
              aes_string("bin", "max"),
              colour = 'paleturquoise4',
              size = 0.5
            ) +
            geom_point(
              data = plotList[[3]],
              aes_string("bin", "max"),
              colour = 'sandybrown',
              size = 0.5
            ) +
            
            coord_flip() +
            scale_x_reverse() +
            ggtitle(paste(
              "Chromosome ", chr, "  (Page", page, ")", sep = ""
            )) +
            facet_wrap(~ lib, nrow = 2, switch = 'x') +
            
            theme(
              legend.position = "none",
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank()
            )
        )
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
#' @param sizeProportional The plotter will divide each fragment into a discrete number of bins of size sizeProportional and distribute reads across them. For example, when set to 200000, the
#' all contigs in each linkage group will be split into 200kb fragments and each fragment will have an equal proportion of the reads derived from that fragment. This ultimately rescales the ideograms
#' such that LGs with the greatest DNA content are the biggest, as opposed to teh LGs with the most contigs. Default value is NULL
#' @param showPage Integer specifying which LG (if plotBy='chr') or libraries (if plotBy='lib') to plot. Useful when not plotting to a file, or when wishing to subset data. Default is NULL
#' @param orderFrame ordered data.frame of contigs (produced by orderAllLinkageGroups). Default is FALSE, where plots will be made from elements in chrTable.
#' @param orientationData data.frame of contig orientations of type OrientationFrame telling which reads to flip Watson and Crick counts
#' @param verbose prints messages to the terminal (default is TRUE)
#'
#' @return ordered contigs in bed format. Depending on options, intermediate files and plots will also be generated
#' @import ggplot2
#' @importFrom data.table rbindlist
#' @importFrom data.table data.table
#' @importFrom data.table tstrsplit
#' @importFrom data.table :=
#' @importFrom gtools mixedsort chr
#' @importFrom S4Vectors DataFrame
#' @aliases ideogramPlot ideogramPlot,StrandReadMatrix,StrandReadMatrix-method,ChrTable,ChrTable-method
#' @export
#' @example inst/examples/ideogramPlot.R
#' @include AllClasses.R
####################################################################################################

setMethod(
  'ideogramPlot',
  signature = signature(
    WatsonFreqList = 'StrandReadMatrix',
    CrickFreqList = 'StrandReadMatrix',
    chrTable = 'ChrTable'
  ),
  definition = ideogramPlot.func
)
