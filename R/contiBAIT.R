#library(contiBAIT)
#path="."
#cluster=12
#dataNames='contiBAIT'
#clusNum=6
#verbose=TRUE
#saveFiles=TRUE
#filter=read.table('ferret_merged_sce_events.bed')
runContiBAIT <- function(path=".", cluster=1, dataNames='contiBAIT', clusNum=1, verbose=TRUE, saveFiles=TRUE)
{
  #Create directory to store all the files
  bamFileList <- list.files(path=path, pattern=".bam$", full.names=TRUE)

  if(verbose){message('RUNNING CONTIBAIT ON ', length(bamFileList), ' BAM FILES!')}


 if(verbose){message('-> Creating read table from bam files [1/6]')}
 # animal.tab <- readStrandCountsFromBAM(path, field=1, dups=TRUE, readLimit=20, freq=TRUE, filter=filter)
 animal.tab <- strandSeqFreqTable(bamFileList, filter=filter, qual=10)

  # subset data with: animal.tab[[1]][which(animal.tab[[2]] < 100)] <- NA
  if(saveFiles){save(animal.tab, file=paste(dataNames, '_table.Rd', sep="")) }
  
  #subset data to only include data above readlimit
  animal.tab[[1]][which(animal.tab[[2]] < readlimit)] <- NA 
  

  if(verbose){message('-> Processing table and filtering data [2/6]')}
  animal.strands <- preprocessStrandTable(animal.tab[[1]], lowQualThreshold=0.8, minLib=6)
  if(saveFiles){save(animal.strands, file=paste(dataNames, '_strands.Rd', sep=""))}

###ADD FUNCTION TO ANALYZE CONTIGS THAT ARE ALWAYS WC!

###RERUN readStrandCounts but pulling out the reads (added filter option)
#lapply(animal.strands[[6]], )

###This is code to do the fused contigs stuff which you shouldn't need for this data and takes a long time, so I've commented it out:
#totalSce <- splitFusedContigs('.', animal.strands[[6]], animal.tab[[1]])
#animal.tab2 <- list(rbind(animal.tab[[1]], totalSce[[1]]), rbind(animal.tab[[2]], totalSce[[2]]))
#animal.strands2 <- preprocessStrandTable(animal.tab2[[1]], lowQualThreshold=0.8)

#2. Analyze reads using CNV caller across all libraries.
#3. for events occuring within ~200kb across > 3? libraries, split ALL libraries into 2 fragments, then rerun readStrandCounts and preprocess to get new calls   

  #create weighting criteria; median of read depth 
  filtWeight <- animal.tab[[2]][which(rownames(animal.tab[[2]]) %in% rownames(animal.strands[[1]])  ),]
  libWeight <- apply(filtWeight, 1, median)
 
  if(verbose){message('-> Clustering data ', cluster, 'x using ', clusNum, ' cores [3/6]')}      
  library(snow)
  slaveNum <- makeCluster(clusNum)
  linkage.groups <- clusterContigs(animal.strands[[2]], randomWeight=libWeight, snowCluster=slaveNum, recluster=cluster, randomise=TRUE, minimumLibraryOverlap=10, similarityCutoff=0.9)
  stopCluster(slaveNum)

  #order linkage groups by biggest first
  linkage.groups <- linkage.groups[order(sapply(linkage.groups, length), decreasing=T)]

  if(saveFiles){ save(linkage.groups, file=paste(dataNames, '_LG_', cluster, 'x_reclust.Rd', sep="") ) }	

  #find groups that are smaller than they should be by finding the SD of linkage group size
  #sdLinkageSize <- ceiling(sd(do.call(rbind, lapply(linkage.groups, function(x) length(x)))))
  #failedLinks <- unlist(linkage.groups[sdLinkageSize:length(linkage.groups)])


  #bob <- lapply(1:length(linkage.groups), function(x) gsub("_.", "", linkage.groups[[x]]))                
  # make orientation calls for each group; WW and CC only
  if(verbose){message('-> Reorienting discordant fragments [4/6]')}
  reorientedGroups <- reorientLinkageGroups(linkage.groups, animal.strands[[2]])
                 
  # perform reorientation of linkage groups that belong together but are misoriented with each other
  # convert the strand table to account for the reorientations
  reorientedTable <- reorientStrandTable(animal.strands[[1]], linkage.groups, reorientedGroups)
  if(saveFiles){save(reorientedTable, file=paste(dataNames, '_reclust_reoriented.Rd', sep=""))}

  if(verbose){message('-> Merging related linkage groups [5/6]')}
  linkage.merged <- mergeLinkageGroups(linkage.groups, reorientedTable)
  linkage.merged <- linkage.merged[order(sapply(linkage.merged, length), decreasing=TRUE)]
  save(linkage.merged, file=paste(dataNames, '_', cluster, 'reclust_merged.Rd', sep="")  )

  if(verbose){message('-> Checking for high quality sex chromosome fragments [6/6]')}
  if(nrow(animal.strands[[3]]) > 1)
  {
    if(verbose){message(paste('  -> ', nrow(animal.strands[[3]]), ' found. Processing.', sep=""))}
    filtWeightSex <- animal.tab[[2]][which(rownames(animal.tab[[2]]) %in% rownames(animal.strands[[3]])  ),]
    libWeightSex <- apply(filtWeightSex, 1, median)
    # cluster the sex groups if any (should only be present in males, assuming either C or W)
    slaveNum <- makeCluster(clusNum)

    linkage.groups.sex <- clusterContigs(animal.strands[[3]], randomWeight=libWeightSex, snowCluster=slaveNum, recluster=cluster, randomise=TRUE, minimumLibraryOverlap=10, similarityCutoff=0.8)
    stopCluster(slaveNum)
    reorientedGroups.sex <- reorientLinkageGroups(linkage.groups.sex, animal.strands[[3]])
    reorientedTable.sex <- reorientStrandTable(animal.strands[[3]], linkage.groups.sex, reorientedGroups.sex)
    linkage.merged.sex <- mergeLinkageGroups(linkage.groups.sex, reorientedTable.sex)
    save(linkage.merged.sex, file=paste(dataNames, '_', cluster, 'reclust_merged_sex.Rd', sep="")  )
  } else {
    if(verbose){message('  -> None found')}
  }

  for( lg in seq(1, length(linkage.merged)))
  {
    if(verbose){message(paste('  -> Ordering fragments in LG', lg, sep=""))}
    outOfOrder <-  orderWithinGroup(linkage.merged[[lg]], animal.strands[[1]], animal.tab[[1]])
    orderFrame <- cbind(LG=rep(lg, length(outOfOrder[[1]])), name=outOfOrder[[1]])
    orderedGroups <- rbind(orderedGroups, orderFrame)
  }
}

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
#plotLGDistances(linkage.merged, animal.strands[[1]], saveFile=paste(dataNames, '_LGdistances_', cluster, 'x_clusters', sep=''))

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

#save(animal.tab, animal.strands, linkage.groups, linkage.merged, reorientedTable, filtWeight, libWeight, filtWeightSex, libWeightSex, file=paste(dataNames, '_all_data_.Rd', sep=""))

#chrTable$chr <- as.character(matrix(unlist(strsplit(chrTable$chr, ':')), ncol=2, byrow=T)[,1])



#lengthPerLG <- sapply(seq(1,length(linkage.merged)), function(x) {sum(as.numeric(matrix(unlist(strsplit(matrix(unlist(strsplit(linkage.merged[[x]], ':')), ncol=2, byrow=T)[,2], '-')), ncol=2, byrow=T)[,2])-as.numeric(matrix(unlist(strsplit( matrix(unlist(strsplit(linkage.merged[[x]], ':')), ncol=2, byrow=T)[,2], '-')), ncol=2, byrow=T)[,1]))})
#lengthPerLG <- lengthPerLG/1000000
#LGname <- c(paste('LG', 1:length(linkage.merged), sep=""))
#maxX <- max(lengthPerLG)
#barplot(lengthPerLG, names.arg=LGname, las=2, xlab='Linkage group', ylab='Length in Mb', ylim=c(0,maxX))
