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

// [[Rcpp::export]]
NumericMatrix TestBivariateSparseNdimArrayReadWrite(
        NumericMatrix coords1, NumericMatrix coords2, NumericVector values) {

    typedef std::vector<double> CoordVec;
    typedef std::pair<CoordVec, CoordVec> CoordVecPair;

    // initialize container w/lexicographically sorted iterator via STL defaults
    SparseNdimArray<CoordVecPair, double, std::map<CoordVecPair, double>> array;

    // fill container
    for(int i=0; i < coords1.nrow(); i++) {
        // munge R coordinate into array's storage format
        std::vector<double> c1(coords1.ncol());
        std::vector<double> c2(coords2.ncol());
        for(int j=0; j< coords1.ncol(); ++j)
            c1[j] = coords1(i,j);
        for(int j=0; j< coords2.ncol(); ++j)
            c2[j] = coords2(i,j);
        // insert coord/value pair into sparse array
        CoordVecPair cp(c1, c2);
        array.set(cp, values(i));
    }

    // package output
    NumericMatrix out = NumericMatrix(array.data.size(),
                                      coords1.ncol() + coords2.ncol() + 1);
    int i=0;
    for(auto iter = array.data.begin(); iter != array.data.end(); ++iter) {
        // extract coordinate info from array's storage format
        CoordVecPair c = iter->first;
        for(int j=0; j < coords1.ncol(); ++j) {
            out(i,j) = c.first[j];
        }
        for(int j=0; j < coords2.ncol(); ++j) {
            out(i, coords1.ncol() + j) = c.second[j];
        }
        // extract value information
        out(i++, out.ncol() - 1) = iter->second;
    }

    // return lexicographically sorted coord/val pairs
    return out;
}

// [[Rcpp::export]]
NumericMatrix TestSparseNdimArrayHash(
        NumericMatrix coords1, NumericMatrix coords2,
        std::vector<double> values) {

    typedef std::vector<int> CoordVec;
    typedef std::pair<CoordVec, CoordVec> CoordVecPair;

    struct VectorHasher {
        int operator()(const CoordVecPair &V) const {
            int hash = V.first.size() + V.second.size();
            for(auto &i : V.first) {
                hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            }
            for(auto &i : V.second) {
                hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            }
            return hash;
        }
    };

    // initialize container w/lexicographically sorted iterator via STL defaults
    SparseNdimArray<CoordVecPair,
            double,
            std::unordered_map<CoordVecPair, double,
            VectorHasher>> array;

    // fill container
    for(int i=0; i < coords1.nrow(); i++) {
        // munge R coordinate into array's storage format
        std::vector<int> c1(coords1.ncol());
        std::vector<int> c2(coords2.ncol());
        for(int j=0; j< coords1.ncol(); ++j)
            c1[j] = coords1(i,j);
        for(int j=0; j< coords2.ncol(); ++j)
            c2[j] = coords2(i,j);
        // insert coord/value pair into sparse array
        CoordVecPair cp(c1, c2);
        array.set(cp, values[i]);
    }

    // package output
    NumericMatrix out = NumericMatrix(array.data.size(),
                                      coords1.ncol() + coords2.ncol() + 1);
    int i=0;
    for(auto iter = array.data.begin(); iter != array.data.end(); ++iter) {
        // extract coordinate info from array's storage format
        CoordVecPair c = iter->first;
        for(int j=0; j < coords1.ncol(); ++j) {
            out(i,j) = c.first[j];
        }
        for(int j=0; j < coords2.ncol(); ++j) {
            out(i, coords1.ncol() + j) = c.second[j];
        }
        // extract value information
        out(i++, out.ncol() - 1) = iter->second;
    }

    // return lexicographically sorted coord/val pairs
    return out;
}