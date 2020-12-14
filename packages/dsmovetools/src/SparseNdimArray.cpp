//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>
#include "SparseNdimArray.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix TestSparseNdimArrayReadWrite(
        NumericMatrix coords, NumericVector values) {

    // initialize container w/lexicographically sorted iterator via STL defaults
    SparseNdimArray<std::vector<double>,
                    double,
                    std::map<std::vector<double>, double>> array;

    // fill container
    for(int i=0; i < coords.nrow(); i++) {
        // munge R coordinate into array's storage format
        std::vector<double> c(coords.ncol());
        for(int j=0; j< coords.ncol(); ++j)
            c[j] = coords(i,j);
        // insert coord/value pair into sparse array
        array.set(c, values(i));
    }

    // package output
    NumericMatrix out = NumericMatrix(array.data.size(), coords.ncol() + 1);
    int i=0;
    for(auto iter = array.data.begin(); iter != array.data.end(); ++iter) {
        // extract coordinate info from array's storage format
        std::vector<double> c = iter->first;
        for(int j=0; j < coords.ncol(); ++j) {
            out(i,j) = c[j];
        }
        // extract value information
        out(i++, out.ncol() - 1) = iter->second;
    }

    // return lexicographically sorted coord/val pairs
    return out;
}