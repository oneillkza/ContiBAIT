#include <Rcpp.h>
using namespace Rcpp;

struct LinkageGroups {
	struct Node { int row, prev; } *nodes;
	int *end, *len, *strand;
	int nRow, nCol, nNodes, nGroups;

	LinkageGroups(int row, int col): nRow(row), nCol(col), nNodes(0), nGroups(0) {
		nodes = new Node[nRow];
		end = new int[nRow];
		len = new int[nRow];
		strand = new int[nRow * nCol];
	}

	~LinkageGroups() {
		delete[] nodes;
		delete[] end;
		delete[] len;
	}

	void addGroup(int row) {
		nodes[nNodes].row = row;
		nodes[nNodes].prev = -1;
		end[nGroups] = nNodes;
		len[nGroups] = 1;
		++nNodes;
		++nGroups;
	}

	void addRow(int group, int row) {
		nodes[nNodes].row = row;
		nodes[nNodes].prev = end[group];
		end[group] = nNodes;
		len[group] += 1;
		++nNodes;
	}

	int* groupStrand(int group) {
		return strand + group * nCol;
	}
};

RcppExport SEXP buildLinkageGroups(SEXP s_contigs, SEXP s_similarityCutoff, SEXP s_minimumLibraryOverlap, SEXP s_verbose, SEXP s_row_names) {
	IntegerMatrix contigs(s_contigs);	
	double similarityCutoff = as<double>(s_similarityCutoff);
	int minimumLibraryOverlap = as<int>(s_minimumLibraryOverlap);
	bool verbose = as<bool>(s_verbose);
	StringVector row_names(s_row_names);

	LinkageGroups LG(contigs.nrow(), contigs.ncol());

	// initial state
	LG.addGroup(0);
	for(int i = 0; i < LG.nCol; ++i) {
		LG.strand[i] = contigs(0, i);
	}

	if(verbose) {
		Rcout << "Initializing contig " << row_names[0] << " [1/" << LG.nRow << "] as LG1\n";
	}  
	
	// clustering loop
	for(int row = 1; row < LG.nRow; ++row) {
		if(verbose) {
			Rcout << "Clustering contig " << row_names[row] << " [" << row+1 << "/" << LG.nRow << "]\n";
		}

		int bestGroup = -1;
		double bestSim = -1.0;
		for(int group = 0; group < LG.nGroups; ++group) {
			// inlined computeSim() (which also does strand[which(strand==3)] <- 2):
			int numCommon = 0, numEqual = 0, *strand = LG.groupStrand(group);
			for(int col = 0; col < LG.nCol; ++col) {
				const int C = contigs(row, col);
				const int L = strand[col];
				const bool common = (C != NA_INTEGER) && (L != NA_INTEGER);
				const bool equal = (C == 2 || C == 3) == (L == 2 || L == 3); 
				numCommon += (int)common;
				numEqual += (int)(common && equal);
			}

			double sim = (numCommon > 0) ? (double)numEqual / (double)numCommon : 0.0;
			if(numCommon >= minimumLibraryOverlap && sim > bestSim) {
				bestSim = sim;
				bestGroup = group;
			}
		}

		if(bestGroup < 0 || bestSim < similarityCutoff) {
			// no good match, make this contig founder of new linkage group
			LG.addGroup(row);
			int *strand = LG.groupStrand(LG.nGroups-1);
			for(int i = 0; i < LG.nCol; ++i) {
				strand[i] = contigs(row, i);
			}
		} else {
			// otherwise add to best matched group and recompute strand state
			if(verbose) {
				Rcout << "  -> Adding " << row_names[row] << " to LG" << bestGroup << " for a cluster of " << LG.len[bestGroup]+1 << "\n";
			}

			LG.addRow(bestGroup, row);
			int *strand = LG.groupStrand(bestGroup);
			// inlined computeConsensus():
			const double minSupport = 0.05;
			for(int col = 0; col < LG.nCol; ++col) {
				int count[4] = { 0, 0, 0, 0 };
				int node = LG.end[bestGroup];
				while(node >= 0) {
					int C = contigs(LG.nodes[node].row, col);
					count[(1 <= C && C <= 3) ? C : 0] += 1;
					node = LG.nodes[node].prev;
				}

				double qcScore = 1.0 - (double)(count[0]) / double(LG.len[bestGroup]);
				if(qcScore < minSupport) {
					strand[col] = NA_INTEGER;
				} else {
					strand[col] = (count[1] >= count[2]) ? 
						((count[1] >= count[3]) ? 1 : 3) :
						((count[2] >= count[3]) ? 2 : 3);
				}
			}
		}		
	}

	// build return list of linkage group vectors
	List ret;
	for(int g = 0; g < LG.nGroups; ++g) {
		IntegerVector group(LG.len[g]);
	
		// link list goes backward so reverse when copying
		int node = LG.end[g], i = LG.len[g] - 1;
		while(node >= 0) {
			group[i--] = LG.nodes[node].row + 1;
			node = LG.nodes[node].prev;
		}

		ret.push_back(group);
	}

	return ret;
}
