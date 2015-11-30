####################################################################################################
#' otherWriteBed -- function to write contig order to bed file
#' @param fileName bed file name to write
#' @param contigOrder a contig ordering
#' @param ligWeight average quality across all libraries for a contig
#' 
#' @return void -bed file written to working directory names by fileName, chromosme and state.
#' 
#' @export
####################################################################################################

otherWriteBed <- function(fileName,contigOrder,libWeight ){

  rtracklayer::export.bed(con = fileName,GRanges(contigOrder, score  =libWeight[contigOrder]))
}