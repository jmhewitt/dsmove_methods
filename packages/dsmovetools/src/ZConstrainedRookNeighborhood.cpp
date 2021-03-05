//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>
#include "ZConstrainedRookNeighborhood.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix TestZConstrainedRookNeighborhood(std::vector<int> dims,
                                               std::vector<int> x,
                                               std::vector<double> zfield,
                                               std::vector<double> zvals) {
    // Parameters:
    //   dims - number of locations along each dimension

    ZConstrainedRookNeighborhood<int, std::vector<int>> rn(dims,
                                                           zfield.data(),
                                                           zvals.data());
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