% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/local_generator.R
\name{local_generator}
\alias{local_generator}
\title{Construct a locally-constrained infinitesimal generator matrix}
\usage{
local_generator(
  locs,
  row_edges,
  col_edges,
  tolocs_by_edge,
  fromlocs_by_edge,
  outedges_by_loc,
  loc_start,
  Xloc,
  betaLoc,
  Xdir,
  betaDir,
  W,
  betaAR
)
}
\arguments{
\item{locs}{vector of locations associated with \code{edges}}

\item{row_edges}{vector of edges to build generator matrix for}

\item{col_edges}{vector of edges to build generator matrix for}

\item{tolocs_by_edge}{location at which each edge in domain ends}

\item{fromlocs_by_edge}{location from which each edge in domain begins}

\item{outedges_by_loc}{edges in domain sorted by their starting location}

\item{loc_start}{index of first out-edge for a location in 
\code{outedges_by_loc}}

\item{Xloc}{covariates for location-based drivers of movement across domain}

\item{betaLoc}{coefficients for location-based drivers of movement}

\item{Xdir}{(edge-based) covariates for directional drivers of movement 
across domain}

\item{betaDir}{coefficients for directional drivers of movement}

\item{W}{unit-vector matrix for whole domain; rows define each edge's 
direction of movement}

\item{betaAR}{scalar coefficient for directional persistence in movement}
}
\description{
Construct a locally-constrained infinitesimal generator matrix
}
