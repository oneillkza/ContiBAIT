####################################################################################################
#' Bar plot all linkage groups, with the true chromosomes of contigs coloured.
#' @param path Path to a directory containing BAM files; every BAM file is read
#' 
#' @return a matrix of lengths of each chromosome (rows) in each linkage group (columns)
#' @export
#' @importFrom colorspace rainbow_hcl
####################################################################################################


contiBAIT <- function(path=".", cluster=1, clusNum=1, verbose=TRUE, saveFiles=TRUE)
{
  #Create directory to store all the files
 
  bamFileLength <- length(list.files(path=path, pattern=".bam$", full.names=TRUE))
  if(verbose){message('RUNNING CONTIBAIT ON ', bamFileLength, ' BAM FILES!')}


 if(verbose){message('-> Creating read table from bam files [1/6]')}
  animal.tab <- readStrandCountsFromBAM(path, field=1, dups=TRUE, readLimit=4, freq=TRUE)
  
  if(saveFiles){save(animal.tab, file='contiBAIT_table.Rd') }
   
  if(verbose){message('-> Processing table and filtering data [2/6]')}
  animal.strands <- preprocessStrandTable(animal.tab[[1]], lowQualThreshold=0.8)
  if(saveFiles){save(animal.strands, file='contiBAIT_strands.Rd')}

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
  slaveNum <- makeCluster(clusNum)
  linkage.groups <- clusterContigs(animal.strands[[2]], randomWeight=libWeight, snowCluster=slaveNum, recluster=cluster, randomise=TRUE, minimumLibraryOverlap=10, similarityCutoff=0.8)
  stopCluster(slaveNum)

  if(saveFiles){ save(linkage.groups, file=paste('contiBAIT_LG_', cluster, 'x_reclust208.Rd', sep="") ) }	

#bob <- lapply(1:length(linkage.groups), function(x) gsub("_.", "", linkage.groups[[x]]))                
  # make orientation calls for each group; WW and CC only
  if(verbose){message('-> Reorienting discordant fragments [4/6]')}
  reorientedGroups <- reorientLinkageGroups(linkage.groups, animal.strands[[2]])
                 
  # perform reorientation of linkage groups that belong together but are misoriented with each other
  # convert the strand table to account for the reorientations
  reorientedTable <- reorientStrandTable(animal.strands[[1]], linkage.groups, reorientedGroups)
  if(saveFiles){save(reorientedTable, file='contiBAIT_reclust_reoriented.Rd')}

  if(verbose){message('-> Merging related linkage groups [5/6]')}
  linkage.merged <- mergeLinkageGroups(linkage.groups, reorientedTable)
  save(linkage.merged, file=paste('contiBAIT_', cluster, 'reclust_merged.Rd', sep="")  )

  if(verbose){message('-> Checking for high quality sex chromosome fragments [6/6]')}
  if(nrow(animal.strands[[3]]) > 1)
  {
    if(verbose){message(paste('  -> ', nrow(animal.strands[[3]]), ' found. Processing.', sep=""))}
    filtWeight <- animal.tab[[2]][which(rownames(animal.tab[[2]]) %in% rownames(animal.strands[[3]])  ),]
    libWeight <- apply(filtWeight, 1, median)`
    # cluster the sex groups if any (should only be present in males, assuming either C or W)
    slaveNum <- makeCluster(clusNum)
    if(verbose){message(paste('  -> ','Check 1.', sep=""))}

    linkage.groups.sex <- clusterContigs(animal.strands[[3]], randomWeight=libWeight, snowCluster=slaveNum, recluster=cluster, randomise=TRUE, minimumLibraryOverlap=10, similarityCutoff=0.8)
    stopCluster(slaveNum)
    reorientedGroups.sex <- reorientLinkageGroups(linkage.groups.sex, animal.strands[[3]])
    reorientedTable.sex <- reorientStrandTable(animal.strands[[3]], linkage.groups.sex, reorientedGroups.sex)
    if(verbose){message(paste('  -> ','Check 2.', sep=""))}
    linkage.merged.sex <- mergeLinkageGroups(linkage.groups.sex, reorientedTable.sex)
    save(linkage.merged.sex, file=paste('contiBAIT_', cluster, 'reclust_merged_sex.Rd', sep="")  )
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

