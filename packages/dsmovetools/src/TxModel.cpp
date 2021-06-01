//
// Created by Joshua Hewitt on 3/15/21.
//

#include "TxModel.h"
#include "RookDirections.h"
#include "RookNeighborhood.h"

using namespace Rcpp;

typedef RookNeighborhood<unsigned int, std::vector<unsigned int>> RN;
typedef RookDirections<unsigned int, std::vector<unsigned int>> RD;

// [[Rcpp::export]]
NumericMatrix TxModelParams(
        std::vector<unsigned int> cur_loc, std::vector<unsigned int> prev_loc,
        std::vector<unsigned int> dims, double betaAR
) {

    typedef unsigned int CoordVal;
    typedef std::vector<CoordVal> CoordVec;

    RN rn(dims);
    RD rd(cur_loc.size());

    TxModel<RN, RD, CoordVec, CoordVal> txm(rn, rd);

    txm.setBetaAR(betaAR);

    txm.constructProbs(cur_loc, prev_loc);

    NumericMatrix out = NumericMatrix(rn.neighborhoodSize(), dims.size() + 1);

    std::vector<double> lp = txm.logProbs();

    std::vector<CoordVec> nbrs = txm.neighbors();

    for(int i=0; i < out.nrow(); ++i) {
        CoordVec nbr = nbrs[i];
        for(int j=0; j < dims.size(); ++j) {
            out(i,j) = nbr[j];
        }
        out(i,dims.size()) = lp[i];
    }

    return out;
}

// [[Rcpp::export]]
std::vector<unsigned int> TxModelSample(
        std::vector<unsigned int> cur_loc, std::vector<unsigned int> prev_loc,
        std::vector<unsigned int> dims, double betaAR
) {

    RN rn(dims);
    RD rd(cur_loc.size());

    TxModel<RN, RD, std::vector<unsigned int>, unsigned int> txm(rn, rd);

    txm.setBetaAR(betaAR);

    txm.constructProbs(cur_loc, prev_loc);

    return txm.sampleNeighbor();
}

// [[Rcpp::export]]
double TxModelLd(
        std::vector<unsigned int> cur_loc, std::vector<unsigned int> prev_loc,
        std::vector<unsigned int> dims, double betaAR,
        std::vector<unsigned int> dst_loc
) {

    RN rn(dims);
    RD rd(cur_loc.size());

    TxModel<RN, RD, std::vector<unsigned int>, unsigned int> txm(rn, rd);

    txm.setBetaAR(betaAR);

    txm.constructProbs(cur_loc, prev_loc);

    return txm.ld(dst_loc);
}