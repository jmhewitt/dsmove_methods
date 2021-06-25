//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>
#include "RookNeighborhood.h"
#include "RookDirections.h"
#include "ZConstrainedRookNeighborhood.h"
#include "TxModel.h"
#include "log_complement.h"
#include "ffar_cached.h"

using namespace Rcpp;

typedef unsigned int IndexType;
// rook neighborhoods on integer grids
typedef RookNeighborhood<IndexType, DualSparseCoordVec::Coord> RN;
typedef ZConstrainedRookNeighborhood<IndexType, DualSparseCoordVec::Coord> ZRN;
typedef RookDirections<IndexType , DualSparseCoordVec::Coord> RD;
typedef TxModel<ZRN, RD, DualSparseCoordVec::Coord, double> TXM;


// [[Rcpp::export]]
NumericMatrix FFRWLightLogConstrainedSelfTxARCached(
        NumericMatrix a0, NumericMatrix a0_prev_coords,
        NumericVector log_a0val,
        std::vector<unsigned int> dims, int steps,
        std::vector<double> surface_heights,
        std::vector<double> domain_heights, double log_self_tx, double betaAR) {

    // initial probability container
    DualSparseCoordVec log_a0;

    // fill initial probability container
    for (int i = 0; i < a0.nrow(); i++) {
        // munge R coordinate into array's storage format
        DualSparseCoordVec::Coord c(a0.ncol());
        DualSparseCoordVec::Coord c0(a0_prev_coords.ncol());
        for (int j = 0; j < a0.ncol(); ++j) {
            c[j] = a0(i, j);
            c0[j] = a0_prev_coords(i, j);
        }
        // insert coord/value pair into sparse array
        DualSparseCoordVec::CoordPair cp(c, c0);
        log_a0.addToActive(cp, log_a0val(i));
    }
    log_a0.swapActive();

    // initialize neighborhood and transition structures
    ZRN zrn(dims, surface_heights.data(), domain_heights.data());
    RD rd(a0.ncol());
    TXM txm(zrn, rd);
    txm.setBetaAR(betaAR);


    class UniformLik {
    public:
        double ll(const DualSparseCoordVec::CoordPair &coord) {
            return 0;
        }
    };

    UniformLik ulik;

    // diffuse mass
    for(IndexType step_cur = 0; step_cur < steps; ++step_cur) {
        // forward-filter all mass from the most recently diffused vector
        diffuseMassSelfTxARDual<IndexType, ZRN, TXM, UniformLik>(
            &log_a0, &zrn, &txm, log_self_tx, &ulik
        );
        // swap state
        log_a0.swapActive();
    }

    // extract final, diffused probability
    NumericMatrix out = NumericMatrix(
            log_a0.data.size(), a0.ncol() + a0_prev_coords.ncol() + 1
    );
    int i = 0;
    for (auto iter = log_a0.data.begin(); iter != log_a0.data.end(); ++iter) {
        // extract coordinate info from array's storage format
        DualSparseCoordVec::Coord c = iter->first.first;
        DualSparseCoordVec::Coord c0 = iter->first.second;
        for (int j = 0; j < a0.ncol(); ++j) {
            out(i, j) = c[j];
            out(i, j + a0.ncol()) = c0[j];
        }
        // extract value information
        out(i++, out.ncol() - 1) = log_a0.inactiveValue(iter->second);
    }

    // return lexicographically sorted coord/val pairs
    return out;
}
