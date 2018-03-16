#ifndef _contiBAIT_computeSim_
#define _contiBAIT_computeSim_

RcppExport  SEXP computeSim (   //Inputs:
								SEXP ScontigStrand, // Strand state for contig
								SEXP SlinkageStrand, // Strand state for linkage group to compare to
								// Both should be numeric vectors
								SEXP SminimumLibraryOverlap //if overlap is less than this, NA will be returned
								//SEXP Ssimilarity // Space for result

);


#endif
