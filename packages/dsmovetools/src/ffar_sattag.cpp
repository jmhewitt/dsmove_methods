//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>
#include "RookNeighborhood.h"
#include "RookDirections.h"
#include "ZConstrainedRookNeighborhood.h"
#include "TxModel.h"
#include "log_complement.h"
#include "ffar.h"
#include "DepthLik.h"
#include "GpsLik.h"

using namespace Rcpp;

// rook neighborhoods on integer grids
typedef RookNeighborhood<IndexType, VectorI> RN;
typedef ZConstrainedRookNeighborhood<IndexType, VectorI> ZRN;
typedef RookDirections<IndexType , VectorI> RD;
typedef TxModel<ZRN, RD, VectorI, double> TXM;


// composite likelihood class
class LocationDepthLik {

    private:

        GpsLikGridded *gps_lik;
        DepthLik *depth_lik;

    public:

        LocationDepthLik(GpsLikGridded &g, DepthLik &d) :
            gps_lik(&g), depth_lik(&d) { };

        double ll(const VectorIPair &coordPair) {
            double g = gps_lik->ll(coordPair.first[1], coordPair.first[0]);
            double d = depth_lik->ll(coordPair.first);
            return g + d;
        };

        void setLikToObs(unsigned int ind) {
            gps_lik->setLikToObs(ind);
            depth_lik->setLikToObs(ind);
        };

};


// [[Rcpp::export]]
double SattagFilteredLL(
    NumericMatrix a0, NumericMatrix a0_prev_coords, NumericVector log_a0val,
    double gps_trunc_alpha, std::vector<double> obs_lons,
    std::vector<double> obs_lats, std::vector<double> obs_semi_majors,
    std::vector<double> obs_semi_minors, std::vector<double> obs_orientations,
    std::vector<unsigned int> obs_depth_bins,
    std::vector<double> lon_gridvals, std::vector<double> lat_gridvals,
    std::vector<unsigned int> dims, std::vector<double> surface_heights,
    std::vector<double> domain_heights, double log_self_tx, double betaAR) {

    // number of steps to integrate over in likelihood
    unsigned int nsteps = obs_lats.size();

    std::vector<unsigned int> dims2 = dims;
    dims2[2] = domain_heights.size();

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
    ZRN zrn2(dims2, surface_heights.data(), domain_heights.data());
    RD rd(a0.ncol());
    TXM txm(zrn, rd);
    txm.setBetaAR(betaAR);

    // initialize GPS, depth, and composite likelihoods
    GpsLikGridded gps_lik(gps_trunc_alpha, obs_lons, obs_lats, obs_semi_majors,
                          obs_semi_minors, obs_orientations, lon_gridvals,
                          lat_gridvals);
    DepthLik depth_lik(obs_depth_bins, zrn2);
    LocationDepthLik sattag_lik(gps_lik, depth_lik);


    /*
     * diffuse initial probability while aggregating likelihood
     */

    LogARMap cur, prev;
    prev = log_a0;

    bool firstError = true;

    IndexType checkmark_interval = .1 * nsteps;

    // diffuse mass across likelihood steps
    for(IndexType step_cur = 0; step_cur < nsteps; ++step_cur) {

        sattag_lik.setLikToObs(step_cur);

        // aggregate marginal likelihood over last diffused probability state
        double ll_step = -std::numeric_limits<double>::infinity();
        for (auto iter = prev.data.begin(); iter != prev.data.end(); ++iter) {
            // get log-likelihood for diffused location
            double ll_state = sattag_lik.ll(iter->first);
            if(std::isfinite(ll_state)) {
                ll_step = log_add(ll_step, iter->second + ll_state);
            }
        }

        if(!std::isfinite(ll_step)) {
            if(firstError) {
                Rcpp::Rcout <<
                    "first infinite marginal likelihood in step_cur: " <<
                    step_cur << std::endl;
                firstError = false;
            }
        }

        ll += ll_step;

//        if(step_cur % checkmark_interval == 0) {
//            std::time_t result = std::time(nullptr);
//            Rcpp::Rcout << "step_cur " << step_cur << " of " << nsteps <<
//                " at " << std::ctime(&result) << std::endl;
//        }
//
//        Rcpp::checkUserInterrupt();

        // diffuse mass (i.e., update prediction distribution)
        diffuseMassSelfTxAR<IndexType, ZRN, TXM, LocationDepthLik>(
                &prev, &cur, &zrn, &txm, log_self_tx, &sattag_lik
        );

        // swap state, to prepare for next diffusion
        prev.data.swap(cur.data);
        cur.data.clear();
    }

    return ll;
}

// [[Rcpp::export]]
NumericMatrix SattagExpandNeighborhood(
        NumericMatrix a0, NumericMatrix a0_prev_coords, NumericVector log_a0val,
        double gps_trunc_alpha, std::vector<double> obs_lons,
        std::vector<double> obs_lats, std::vector<double> obs_semi_majors,
        std::vector<double> obs_semi_minors, std::vector<double> obs_orientations,
        std::vector<unsigned int> obs_depth_bins,
        std::vector<double> lon_gridvals, std::vector<double> lat_gridvals,
        std::vector<unsigned int> dims, std::vector<double> surface_heights,
        std::vector<double> domain_heights, unsigned int obs_ind,
        unsigned int nsteps) {

    std::vector<unsigned int> dims2 = dims;
    dims2[2] = domain_heights.size();

    // initial probability container
    LogARMap log_a0;

    double betaAR = 0;
    double log_self_tx = -.69;

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
    ZRN zrn2(dims2, surface_heights.data(), domain_heights.data());
    RD rd(a0.ncol());
    TXM txm(zrn, rd);
    txm.setBetaAR(betaAR);

    // initialize GPS, depth, and composite likelihoods
    GpsLikGridded gps_lik(gps_trunc_alpha, obs_lons, obs_lats, obs_semi_majors,
                          obs_semi_minors, obs_orientations, lon_gridvals,
                          lat_gridvals);
    DepthLik depth_lik(obs_depth_bins, zrn2);
    LocationDepthLik sattag_lik(gps_lik, depth_lik);


    /*
     * expand initial neighborhood
     */

    sattag_lik.setLikToObs(obs_ind);

    LogARMap cur, prev;
    prev = log_a0;

    unsigned int nbhd_size = prev.data.size();

    // diffuse mass across until neighborhood stabilizes
    for(IndexType step_cur = 0; step_cur < nsteps; ++step_cur) {

        Rcpp::checkUserInterrupt();

        // diffuse mass (i.e., update prediction distribution)
        diffuseMassSelfTxAR<IndexType, ZRN, TXM, LocationDepthLik>(
                &prev, &cur, &zrn, &txm, log_self_tx, &sattag_lik
        );

        // swap state, to prepare for next diffusion
        prev.data.swap(cur.data);
        cur.data.clear();

        if(prev.data.size() == nbhd_size) {
            break;
        } else {
            nbhd_size = prev.data.size();
        }
    }

    // package output
    NumericMatrix out = NumericMatrix(prev.data.size(), a0.ncol() * 2);
    int i=0;
    for(auto iter = prev.data.begin(); iter != prev.data.end(); ++iter) {
        // extract coordinate info from array's storage format
        VectorIPair c = iter->first;
        for(int j=0; j < a0.ncol(); ++j) {
            out(i,j) = c.first[j];
        }
        for(int j=0; j < a0.ncol(); ++j) {
            out(i, a0.ncol() + j) = c.second[j];
        }
        ++i;
    }

    return out;
}