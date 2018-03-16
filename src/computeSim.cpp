#include <R_ext/Print.h>
#include <Rcpp.h>
#include <boost/dynamic_bitset.hpp>
#include "computeSim.h" 

using namespace Rcpp ;
using namespace boost;

/*****************************************************************************
 * Function to compute similarity between two contigs
 * 
 ****************************************************************************/
 
 
RcppExport  SEXP computeSim (   //Inputs:
								SEXP ScontigStrand, // Strand state for contig
								SEXP SlinkageStrand, // Strand state for linkage group to compare to
								// Both should be numeric vectors, with 1 denoting WW and 2 denoting WC
								SEXP SminimumLibraryOverlap //if overlap is less than this, NA will be returned
								//SEXP Ssimilarity // Space for result

)
{
	//Convert from SEXPs to R data types:
	IntegerVector contigStrand (ScontigStrand);
	IntegerVector linkageStrand (SlinkageStrand);
	double minimumLibraryOverlap = as<double>(SminimumLibraryOverlap);
	
	int vectorSize = contigStrand.size();
	
	dynamic_bitset<> contigBits(vectorSize);
	dynamic_bitset<> linkageBits(vectorSize);
	dynamic_bitset<> contigNAs(vectorSize);
	dynamic_bitset<> linkageNAs(vectorSize);
	
	//Fill contig bitsets:
	for (int i=0; i<vectorSize; i++)
	{
		if(contigStrand[i]==2)
			contigBits[i] = 1;
			
		if(contigStrand[i]==NA_INTEGER)
			contigNAs[i] = 1;	
	}
	
	//Fill linkage bitsets:
	for (int i=0; i<vectorSize; i++)
	{
		if(linkageStrand[i]==2)
			linkageBits[i] = 1;
			
		if(linkageStrand[i]==NA_INTEGER)
			linkageNAs[i] = 1;	
		
	}
	
	//Compute common bits:
	dynamic_bitset<> commonBits = ~(contigNAs|linkageNAs);
	//Rprintf("Linkage NAs=%i, Contig NAs=%i, Common NAs=%i, Common bits=%i\n", linkageNAs.count(), contigNAs.count(), numNAs.count(), commonBits.count());
	dynamic_bitset<> equalBits = ~(contigBits^linkageBits);
	
	//Rprintf("Linkage bits=%i, Contig bits=%i, Equal bits=%i\n", linkageBits.count(), contigBits.count(), equalBits.count());
	
	int numCommon = commonBits.count();
	
	if(numCommon < minimumLibraryOverlap)
		return wrap(NA_REAL);
	
	int numEqual = (equalBits&commonBits).count();
	
	
	double similarity = (double) numEqual / (double) numCommon;
	//Rprintf("numEqual=%i, numCommon=%i, similarity=%f\n", numEqual, numCommon, similarity);
	return wrap(similarity);
	
}

 
 