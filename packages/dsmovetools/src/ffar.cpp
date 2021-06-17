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


template <typename size_type, typename Neighborhood, typename TxMod,
        typename SrcLik>
void diffuseMassSelfTxAR(LogARMap *src, LogARMap *dst, Neighborhood *nbhd,
                         TxMod *txmod, double log_self_tx, SrcLik *srclik) {
    // Parameters
    //  log_self_tx - log of self-transition probability

    // compute probability of making a transition
    double log_tx = log_complement(log_self_tx);

    // log-mass of dst vector
    double dst_mass = -std::numeric_limits<double>::infinity();

    // forward-filter all mass from the src vector to the dst vector
    auto src_mass_entry = src->data.begin();
    auto src_mass_end = src->data.end();
    for(src_mass_entry; src_mass_entry != src_mass_end; ++src_mass_entry) {
        // get log-likelihood for current location
        double ll = srclik->ll(src_mass_entry->first);
        // forward-diffuse from location if observation likelihood is finite
        if(std::isfinite(ll)) {
            // find reachable nodes for current location
            nbhd->setCenter(src_mass_entry->first.first);
            size_type nnbrs = nbhd->neighborhoodSize();
            // add mass for self-transition to dst vector
            double m = src_mass_entry->second + log_self_tx;
            dst->add(src_mass_entry->first, m);
            // aggregate total output mass
            dst_mass = log_add(dst_mass, m);
            // compute and extract stepwise transition probabilities
            txmod->constructProbs(src_mass_entry->first.first,
                                  src_mass_entry->first.second);
            std::vector<VectorI> nbrs = txmod->neighbors();
            std::vector<double> lp = txmod->logProbs();
            // weight stepwise tx. probs. by non self-tx prob and mass of atom
            for(size_type i = 0; i < nnbrs; ++i) {
                // mass of transitioning to dst vector
                double m = src_mass_entry->second + lp[i] + log_tx;
                dst->add(
                    std::pair<VectorI, VectorI>(
                        nbrs[i], src_mass_entry->first.first
                    ),
                    m
                );
                // aggregate total output mass
                dst_mass = log_add(dst_mass, m);
            }
        }
    }
    // standardize distribution
    if(dst_mass != 0) {
        auto dst_mass_entry = dst->data.begin();
        auto dst_mass_end = dst->data.end();
        for(dst_mass_entry; dst_mass_entry != dst_mass_end; ++dst_mass_entry) {
            dst_mass_entry->second -= dst_mass;
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
template<typename size_type, typename Neighborhood, typename TxMod,
        typename SrcLik>
LogARMap ffrw_light_selftx_ar(const VectorI &dims,
                              const LogARMap &a0,
                              const IndexType steps, Neighborhood &nbhd,
                              TxMod &txmod, double log_self_tx,
                              SrcLik &srclik) {

    // initialize forward filtering vectors and initial mass
    LogARMap cur, prev;
    prev = a0;

    // diffuse mass
    for(IndexType step_cur = 0; step_cur < steps; ++step_cur) {
        // forward-filter all mass from the most recently diffused vector
        diffuseMassSelfTxAR<size_type, Neighborhood, TxMod>(
                &prev, &cur, &nbhd, &txmod, log_self_tx, &srclik
        );
        // swap state
        prev.data.swap(cur.data);
        cur.data.clear();
    }

    return prev;
};

template<typename Index>
class UniformLik {
    public:
        double ll(const Index &coord) {
            return 0;
        }
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

    UniformLik<VectorIPair> ulik;

    // diffuse initial probability
    LogARMap log_af = ffrw_light_selftx_ar<IndexType, ZRN, TXM,
        UniformLik<VectorIPair>>(
            dims, log_a0, steps, zrn, txm, log_self_tx, ulik
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

class ObsIndicatorLik {
private:
    VectorI obsloc;
    VectorI oob;
public:
    ObsIndicatorLik(const VectorI &dims) {
        oob = dims;
        for(auto it = oob.begin(); it != oob.end(); ++it) {
            *(it) = *(it) + 1;
        }
    }
    void setObsLoc(const VectorI &coord) {
        obsloc = coord;
    }
    double ll(const VectorIPair &coordPair) {
        // do not match observations with nan; these are unconstrained data
        if(obsloc == oob) {
            return 0;
        }
        // return -Inf if coordinates do not match
        if(coordPair.first == obsloc) {
            return 0;
        } else {
            return -std::numeric_limits<double>::infinity();
        }
    }
};

// [[Rcpp::export]]
double ARFilteredLL(
        NumericMatrix a0, NumericMatrix a0_prev_coords,
        NumericMatrix obs_coords, NumericVector log_a0val,
        std::vector<unsigned int> dims, std::vector<double> surface_heights,
        std::vector<double> domain_heights, double log_self_tx, double betaAR) {

    // number of steps to integrate over in likelihood
    unsigned int nsteps = obs_coords.nrow();

    // initialize log-likelihood
    double ll = 0;

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

    // initialize likelihood
    ObsIndicatorLik lik(dims);


    /*
     * diffuse initial probability while aggregating likelihood
     */

    LogARMap cur, prev;
    prev = log_a0;

    // diffuse mass across likelihood steps
    for(IndexType step_cur = 0; step_cur < nsteps; ++step_cur) {

        // set observation likelihood, for diffusion
        VectorI coord(a0.ncol());
        for(int j = 0; j < a0.ncol(); ++j) {
            if(std::isnan(obs_coords(step_cur,j))) {
                coord[j] = dims[j] + 1;
            } else {
                coord[j] = obs_coords(step_cur,j);
            }
        }
        lik.setObsLoc(coord);

        // aggregate marginal likelihood over last diffused probability state
        double ll_step = -std::numeric_limits<double>::infinity();
        for (auto iter = prev.data.begin(); iter != prev.data.end(); ++iter) {
            // get log-likelihood for diffused location
            double ll_state = lik.ll(iter->first);
            if(std::isfinite(ll_state)) {
                ll_step = log_add(ll_step, iter->second + ll_state);
            }
        }
        ll += ll_step;

        // diffuse mass (i.e., update prediction distribution)
        diffuseMassSelfTxAR<IndexType, ZRN, TXM>(
                &prev, &cur, &zrn, &txm, log_self_tx, &lik
        );

        // swap state, to prepare for next diffusion
        prev.data.swap(cur.data);
        cur.data.clear();
    }

    return ll;
}

// [[Rcpp::export]]
Rcpp::List ARPredDist(
        NumericMatrix a0, NumericMatrix a0_prev_coords,
        NumericMatrix obs_coords, NumericVector log_a0val,
        std::vector<unsigned int> dims, std::vector<double> surface_heights,
        std::vector<double> domain_heights, double log_self_tx, double betaAR,
        std::vector<unsigned int> pred_steps) {

    // number of steps to integrate over in likelihood
    unsigned int nsteps = obs_coords.nrow();

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

    // initialize likelihood
    ObsIndicatorLik lik(dims);

    // initialize container for predictive distributions
    std::vector<LogARMap> pred_distns;
    pred_distns.reserve(pred_steps.size());

    /*
     * diffuse initial probability while aggregating likelihood
     */

    LogARMap cur, prev;
    prev = log_a0;

    auto pred_step_iter = pred_steps.begin();
    auto pred_step_end = pred_steps.end();

    // diffuse mass across likelihood steps
    for(IndexType step_cur = 0; step_cur < nsteps; ++step_cur) {

        // set observation likelihood, for diffusion
        VectorI coord(a0.ncol());
        for(int j = 0; j < a0.ncol(); ++j) {
            if(std::isnan(obs_coords(step_cur,j))) {
                coord[j] = dims[j] + 1;
            } else {
                coord[j] = obs_coords(step_cur,j);
            }
        }
        lik.setObsLoc(coord);

        // save prediction distribution, and check end condition
        if(pred_step_iter != pred_step_end) {
            if(*pred_step_iter == step_cur) {
                pred_distns.push_back(prev);
                pred_step_iter++;
            }
        } else {
            break;
        }

        // diffuse mass (i.e., update prediction distribution)
        diffuseMassSelfTxAR<IndexType, ZRN, TXM>(
                &prev, &cur, &zrn, &txm, log_self_tx, &lik
        );

        // swap state, to prepare for next diffusion
        prev.data.swap(cur.data);
        cur.data.clear();
    }

    //
    // package results
    //

    Rcpp::List res(pred_steps.size());

    Rcpp::Rcout << pred_distns.size() << std::endl;

    for(IndexType it = 0; it < pred_distns.size(); ++it) {

        LogARMap pred_dist = pred_distns[it];

        Rcpp::Rcout << " " << pred_dist.data.size() << std::endl;

        NumericMatrix out = NumericMatrix(pred_dist.data.size(),
                                          a0.ncol() * 2 + 1);

        int i =0;
        auto iter = pred_dist.data.begin();
        auto end = pred_dist.data.end();
        for(iter; iter != end; ++iter) {
            // extract coordinate info from array's storage format
            VectorIPair c = iter->first;
            for(int j=0; j < a0.ncol(); ++j) {
                out(i,j) = c.first[j];
            }
            for(int j=0; j < a0.ncol(); ++j) {
                out(i, a0.ncol() + j) = c.second[j];
            }
            // extract value information
            out(i++, out.ncol() - 1) = iter->second;
        }

        res[it] = out;
    }

    return res;
}