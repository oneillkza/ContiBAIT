####################################################################################################
#' function to write contig order to BED file
#' @param fileName bed file name to write
#' @param contigOrder a contig ordering
#' @param libWeight average quality across all libraries for a contig
#' @importFrom rtracklayer export.bed
#' @import GRanges
#' @return void; BED file written to working directory names by fileName, chromosome and state.
#' 
#' @export
####################################################################################################

writeBed <- function(fileName,contigOrder,libWeight ){

  export.bed(con = fileName,GRanges(contigOrder, score  =libWeight[contigOrder]))
}