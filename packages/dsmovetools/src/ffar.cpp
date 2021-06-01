//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>
#include "SparseNdimArray.h"
#include "RookNeighborhood.h"
#include "RookDirections.h"
#include "ZConstrainedRookNeighborhood.h"
#include "TxModel.h"

#include "log_complement.h"

using namespace Rcpp;

// let's use ints to index coordinates
typedef unsigned int IndexType;
// for specifying integer coordinates with grid dimension assigned at runtime
typedef std::vector<IndexType> VectorI;
// for specifying pairs of integer coordinates
typedef std::pair<VectorI, VectorI> VectorIPair;
// sparse spatial storage implemented via trees
typedef SparseNdimArrayLog<VectorIPair,
                           double,
                           std::map<VectorIPair, double>> LogARMap;
// rook neighborhoods on integer grids
typedef RookNeighborhood<IndexType, VectorI> RN;
typedef ZConstrainedRookNeighborhood<IndexType, VectorI> ZRN;
typedef RookDirections<IndexType , VectorI> RD;
typedef TxModel<ZRN, RD, VectorI, double> TXM;


template <typename size_type, typename Neighborhood, typename TxMod>
void diffuseMassSelfTxAR(LogARMap *src, LogARMap *dst, Neighborhood *nbhd,
                         TxMod *txmod, double log_self_tx) {
    // Parameters
    //  log_self_tx - log of self-transition probability

    // compute probability of making a transition
    double log_tx = log_complement(log_self_tx);

    // forward-filter all mass from the src vector to the dst vector
    auto src_mass_entry = src->data.begin();
    auto src_mass_end = src->data.end();
    for(src_mass_entry; src_mass_entry != src_mass_end; ++src_mass_entry) {
        // find reachable nodes for current location
        nbhd->setCenter(src_mass_entry->first.first);
        size_type nnbrs = nbhd->neighborhoodSize();
        // add mass for self-transition to dst vector
        dst->addScaled(src_mass_entry->first, src_mass_entry->second,
                       log_self_tx);
        // compute and extract stepwise transition probabilities
        txmod->constructProbs(src_mass_entry->first.first,
                              src_mass_entry->first.second);
        std::vector<VectorI> nbrs = txmod->neighbors();
        std::vector<double> lp = txmod->logProbs();
        // weight stepwise tx. probs. by non self-tx prob and mass of atom
        for(size_type i = 0; i < nnbrs; ++i) {
            dst->addScaled(
                std::pair<VectorI, VectorI>(nbrs[i],
                                            src_mass_entry->first.first),
                src_mass_entry->second,
                lp[i] + log_tx
            );
        }
    }
}

/**
 * Forward filtering a random walk along a grid, without storing all
 * intermediate distributions.  Forward filtering allows self-transitions.
 *
 * @param dims specify the number of locations along each dimension in grid
 * @param a0 initial probability mass vector
 * @param steps number of forward-diffusion steps to take
 * @param nbhd Class that defines the neighborhood for arbitrary locations
 * @return (sparse) diffused mass vectors
 */
template<typename size_type, typename Neighborhood, typename TxMod>
LogARMap ffrw_light_selftx_ar(const VectorI &dims,
                              const LogARMap &a0,
                              const IndexType steps, Neighborhood &nbhd,
                              TxMod &txmod, double log_self_tx) {

    // initialize forward filtering vectors and initial mass
    LogARMap cur, prev;
    prev = a0;

    // diffuse mass
    for(IndexType step_cur = 0; step_cur < steps; ++step_cur) {
        // forward-filter all mass from the most recently diffused vector
        diffuseMassSelfTxAR<size_type, Neighborhood, TxMod>(
                &prev, &cur, &nbhd, &txmod, log_self_tx
        );
        // swap state
        prev.data.swap(cur.data);
        cur.data.clear();
    }

    return prev;
};

// [[Rcpp::export]]
NumericMatrix FFRWLightLogConstrainedSelfTxAR(
        NumericMatrix a0, NumericMatrix a0_prev_coords,
        NumericVector log_a0val,
        std::vector<unsigned int> dims, int steps,
        std::vector<double> surface_heights,
        std::vector<double> domain_heights, double log_self_tx, double betaAR) {


    // initial probability container
    LogARMap log_a0;

    // fill initial probability container
    for (int i = 0; i < a0.nrow(); i++) {
        // munge R coordinate into array's storage format
        VectorI c(a0.ncol());
        VectorI c0(a0_prev_coords.ncol());
        for (int j = 0; j < a0.ncol(); ++j) {
            c[j] = a0(i, j);
            c0[j] = a0_prev_coords(i, j);
        }
        // insert coord/value pair into sparse array
        log_a0.set(VectorIPair(c, c0), log_a0val(i));
    }

    // initialize neighborhood and transition structures
    ZRN zrn(dims, surface_heights.data(), domain_heights.data());
    RD rd(a0.ncol());
    TXM txm(zrn, rd);
    txm.setBetaAR(betaAR);

    // diffuse initial probability
    LogARMap log_af = ffrw_light_selftx_ar<IndexType, ZRN, TXM>(
            dims, log_a0, steps, zrn, txm, log_self_tx
    );

    // extract final, diffused probability
    NumericMatrix out = NumericMatrix(
            log_af.data.size(), a0.ncol() + a0_prev_coords.ncol() + 1
    );
    int i = 0;
    for (auto iter = log_af.data.begin(); iter != log_af.data.end(); ++iter) {
        // extract coordinate info from array's storage format
        VectorI c = iter->first.first;
        VectorI c0 = iter->first.second;
        for (int j = 0; j < a0.ncol(); ++j) {
            out(i, j) = c[j];
            out(i, j + a0.ncol()) = c0[j];
        }
        // extract value information
        out(i++, out.ncol() - 1) = iter->second;
    }

    // return lexicographically sorted coord/val pairs
    return out;
}