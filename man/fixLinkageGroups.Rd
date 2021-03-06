% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fixLinkageGroups.R
\docType{methods}
\name{fixLinkageGroups,ContigOrdering,StrandStateList,LinkageGroupList-method}
\alias{fixLinkageGroups,ContigOrdering,StrandStateList,LinkageGroupList-method}
\alias{fixLinkageGroups}
\alias{fixLinkageGroups-ContigOrdering-StrandStateList-method}
\title{fixLinkageGroups -- searches for discrepancies within ordered contigs to highlight erroneously merged fragments.}
\usage{

  \S4method{fixLinkageGroups}{ContigOrdering,StrandStateList,LinkageGroupList}(contigOrdering,
  orderFrame, linkageGroupList, whichLG = NULL, relatedCutOff = 0.6,
  verbose = TRUE)
}
\arguments{
\item{contigOrdering}{a data.frame of ordered contigs with linkage group names of class ContigOrdering}

\item{orderFrame}{a list of StrandStateMatrix elements of class StrandStateList}

\item{linkageGroupList}{List of vectors containing names of contigs belonging to each LG of type LinkageGroupList.}

\item{whichLG}{vector of integers specifying the element(s) of linkageGroupList to be ordered (i.e. which specific linkage groups to try to order). Default is all LGs.}

\item{relatedCutOff}{The minimal dissimilarity between adjacent contigs to subset a linkage group into multiple smaller groups. Default is 0.6}

\item{verbose}{Outputs information to the terminal. Default is TRUE.}
}
\value{
a LinkageGroupList with erroneously clustered contigs seperated into their own groups
}
\description{
fixLinkageGroups -- searches for discrepancies within ordered contigs to highlight erroneously merged fragments.
}
