
####################################################################################################
#' extractReads -- function to make bed file from all reads occupying uncalled regions from createOrientationBed
#' 
#' @param locations  data.frame of strand frequencies created by createReadLocations
#' @param uncalled  data.frame of uncalled regions created createOrientationBed
#' 
#' @return data.frame in bed format of all uncalled reads 
#' 
#' @export
####################################################################################################


extractReads <- function(locations, uncalled)
{
  uncallBed <- data.frame(chromosome=vector(), start=vector(), end=vector(), strand=vector(), stringsAsFactors=FALSE)
  for( row in seq(1:nrow(uncalled)))
  {
    locations <- locations[order(locations$location),]
    uncalledReads <- locations[which(locations$location > uncalled$start[row] & locations$location < uncalled$end[row]),]
    uncallBed <- rbind(uncallBed, uncalledReads)
  }
  return(uncallBed)
}
