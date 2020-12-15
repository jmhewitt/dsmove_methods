//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>
#include "RookNeighborhood.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix TestRookNeighborhood(std::vector<int> dims, std::vector<int> x) {
    // Parameters:
    //   dims - number of locations along each dimension

    RookNeighborhood<int, std::vector<int>> rn(dims);
    rn.setCenter(x);

    NumericMatrix out = NumericMatrix(rn.neighborhoodSize(), dims.size());

    for(int i=0; i < out.nrow(); ++i) {
        std::vector<int> nbr = rn.nextNeighbor();
        for(int j=0; j < dims.size(); ++j) {
            out(i,j) = nbr[j];
        }
    }

    // return neighborhood
    return out;
}