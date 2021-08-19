//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>
#include "CTDS2DDomain.h"
#include "CTDS2DProbs.h"
#include "ffar.h"
#include "GpsLik.h"
#include "DepthLik.h"

using namespace Rcpp;

// composite likelihood class
class LocationDepthLik {

    private:

        GpsLik *gps_lik;
        DepthLik *depth_lik;

    public:

        LocationDepthLik(GpsLik &g, DepthLik &d) :
            gps_lik(&g), depth_lik(&d) { };

        double ll(const CTDS2DState& state) {
            double d = depth_lik->ll(state.surface_height);
            if(!std::isfinite(d)) {
                return d;
            }
            double g = gps_lik->ll(state.lon_to, state.lat_to);
            return g + d;
        };

        void setLikToObs(unsigned int ind) {
            gps_lik->setLikToObs(ind);
            depth_lik->setLikToObs(ind);
        };

};


// [[Rcpp::export]]
double SattagFilteredLL(
    NumericMatrix init_dsts, NumericMatrix init_srcs,
    std::vector<double> init_log_probs,
    double gps_trunc_alpha, std::vector<double> obs_lons,
    std::vector<double> obs_lats, std::vector<double> obs_semi_majors,
    std::vector<double> obs_semi_minors, std::vector<double> obs_orientations,
    std::vector<double> obs_depths,
    std::vector<double> lon_gridvals, std::vector<double> lat_gridvals,
    std::vector<double> surface_heights, double min_elevation,
    double max_elevation, double log_self_tx, double betaAR
) {

    // number of steps to integrate over in likelihood
    unsigned int nsteps = obs_lats.size();

    // initialize log-likelihood
    double ll = 0;

    // initialize CTMC state space probability vector
    CTDS2DDomain pvec(lon_gridvals, lat_gridvals, surface_heights);

    // filter state space
    CTDS2DStateElevationFilter height_filter(min_elevation, max_elevation);
    pvec.filterStates(height_filter);

    // fill initial probability container
    for(unsigned int ind = 0; ind < init_dsts.nrow(); ++ind) {
        pvec.set(init_srcs(ind, 0), init_srcs(ind, 1), init_dsts(ind, 0),
                 init_dsts(ind, 1), init_log_probs[ind]);
    }

    // finalize initial probability vector
    pvec.swapActive();

    // set transition parameters
    TxProbs txmod;
    txmod.setBetaAR(betaAR);


    // initialize GPS and composite likelihoods
    GpsLik gps_lik(gps_trunc_alpha, obs_lons, obs_lats, obs_semi_majors,
                   obs_semi_minors, obs_orientations);
    DepthLik depth_lik(obs_depths);
    LocationDepthLik sattag_lik(gps_lik, depth_lik);


    /*
     * diffuse initial probability while aggregating likelihood
     */

    bool firstError = true;

    unsigned int checkmark_interval = .1 * nsteps;

    // diffuse mass across likelihood steps
    for(unsigned int step_cur = 0; step_cur < nsteps; ++step_cur) {

        sattag_lik.setLikToObs(step_cur);


        // aggregate marginal likelihood over last diffused probability state
        bool finiteStateWeight = false;
        bool finiteLikObs = false;
        double ll_step = -std::numeric_limits<double>::infinity();
        auto state_end = pvec.end();
        for(auto state_it = pvec.begin(); state_it != state_end; ++state_it) {
            // extract state weight
            double w = pvec.logProbCached(*state_it);
            if(std::isfinite(w)) {
                finiteStateWeight = true;
                // get log-likelihood for diffused location
                double ll_state = sattag_lik.ll(*state_it);
                if(std::isfinite(ll_state)) {
                    finiteLikObs = true;
                    ll_step = log_add(ll_step, ll_state + w);
                }
            }
        }

        if(!std::isfinite(ll_step)) {
            if(firstError) {
                // print error message
                Rcpp::Rcout <<
                    "first infinite marginal likelihood in step_cur: " <<
                    step_cur << ", ll_step: " << ll_step <<
                    ", finiteStateWeight: " << finiteStateWeight <<
                    ", finiteLikObs: " << finiteLikObs <<
                    std::endl;
                firstError = false;
                // read out state
                pvec.unswapActive();
                NumericMatrix state = pvec.toNumericMatrix();
                pvec.swapActive();
                // dump state to disk
                Environment base("package:base");
                Function saveRDS = base["saveRDS"];
                std::string fname = std::string("pvec_state__step_cur_").
                    append(std::to_string(step_cur)).append("__log_self_tx_").
                    append(std::to_string(log_self_tx)).append("__betaAR_").
                    append(std::to_string(betaAR)).append(".rds");
                saveRDS(state, Named("file",fname));
            }
        }

        ll += ll_step;

        Rcpp::checkUserInterrupt();

        // diffuse mass (i.e., update prediction distribution)
        diffuseMassPred<LocationDepthLik>(&pvec, &txmod, log_self_tx, &sattag_lik);
        // swap state
        pvec.swapActive();
    }

    return ll;
}

// [[Rcpp::export]]
std::vector<NumericMatrix> SattagPredDist(
        NumericMatrix init_dsts, NumericMatrix init_srcs,
        std::vector<double> init_log_probs,
        double gps_trunc_alpha, std::vector<double> obs_lons,
        std::vector<double> obs_lats, std::vector<double> obs_semi_majors,
        std::vector<double> obs_semi_minors, std::vector<double> obs_orientations,
        std::vector<double> obs_depths,
        std::vector<double> lon_gridvals, std::vector<double> lat_gridvals,
        std::vector<double> surface_heights, double min_elevation,
        double max_elevation, double log_self_tx, double betaAR,
        std::vector<unsigned int> pred_steps) {

    // Return for prediction at pred_steps given data up to pred_steps - 1

    // number of steps to integrate over in likelihood
    unsigned int nsteps = obs_lats.size();

    // initialize container for predictive distributions
    std::vector<NumericMatrix> pred_distns;
    pred_distns.reserve(pred_steps.size());

    // initialize CTMC state space probability vector
    CTDS2DDomain pvec(lon_gridvals, lat_gridvals, surface_heights);

    // filter state space
    CTDS2DStateElevationFilter height_filter(min_elevation, max_elevation);
    pvec.filterStates(height_filter);

    // fill initial probability container
    for(unsigned int ind = 0; ind < init_dsts.nrow(); ++ind) {
        pvec.set(init_srcs(ind, 0), init_srcs(ind, 1), init_dsts(ind, 0),
                 init_dsts(ind, 1), init_log_probs[ind]);
    }

    // finalize initial probability vector
    pvec.swapActive();

    // set transition parameters
    TxProbs txmod;
    txmod.setBetaAR(betaAR);

    // initialize GPS and composite likelihoods
    GpsLik gps_lik(gps_trunc_alpha, obs_lons, obs_lats, obs_semi_majors,
                   obs_semi_minors, obs_orientations);
    DepthLik depth_lik(obs_depths);
    LocationDepthLik sattag_lik(gps_lik, depth_lik);


    /*
     * diffuse initial probability while aggregating likelihood
     */

    auto pred_step_iter = pred_steps.begin();
    auto pred_step_end = pred_steps.end();

    // diffuse mass across likelihood steps
    for(unsigned int step_cur = 0; step_cur < nsteps; ++step_cur) {

        sattag_lik.setLikToObs(step_cur);

        // save prediction distribution
        if(pred_step_iter != pred_step_end) {
            if(*pred_step_iter == step_cur) {
                pvec.unswapActive();
                pred_distns.push_back(pvec.toNumericMatrix());
                pred_step_iter++;
                pvec.swapActive();
            }
        }

        // check end condition
        if(pred_step_iter == pred_step_end) {
            break;
        }

        Rcpp::checkUserInterrupt();

        // diffuse mass (i.e., update prediction distribution)
        diffuseMassPred<LocationDepthLik>(&pvec, &txmod, log_self_tx, &sattag_lik);
        // swap state
        pvec.swapActive();
    }

    return pred_distns;
}

// [[Rcpp::export]]
std::vector<NumericMatrix> BackInfoFilteringDist(
        NumericMatrix init_dsts, NumericMatrix init_srcs,
        std::vector<double> init_log_probs,
        double gps_trunc_alpha, std::vector<double> obs_lons,
        std::vector<double> obs_lats, std::vector<double> obs_semi_majors,
        std::vector<double> obs_semi_minors, std::vector<double> obs_orientations,
        std::vector<double> obs_depths,
        std::vector<double> lon_gridvals, std::vector<double> lat_gridvals,
        std::vector<double> surface_heights, double min_elevation,
        double max_elevation, double log_self_tx, double betaAR,
        std::vector<unsigned int> pred_steps) {

    // number of steps to integrate over in likelihood
    unsigned int nsteps = obs_lats.size();

    // initialize container for predictive distributions
    std::vector<NumericMatrix> pred_distns;
    pred_distns.reserve(pred_steps.size());

    // initialize CTMC state space probability vector
    CTDS2DDomain pvec(lon_gridvals, lat_gridvals, surface_heights);

    // filter state space
    CTDS2DStateElevationFilter height_filter(min_elevation, max_elevation);
    pvec.filterStates(height_filter);

    // fill initial probability container
    for(unsigned int ind = 0; ind < init_dsts.nrow(); ++ind) {
        pvec.set(init_srcs(ind, 0), init_srcs(ind, 1), init_dsts(ind, 0),
                 init_dsts(ind, 1), init_log_probs[ind]);
    }

    // finalize initial probability vector
    pvec.swapActive();

    // set transition parameters
    TxProbs txmod;
    txmod.setBetaAR(betaAR);

    // initialize GPS and composite likelihoods
    GpsLik gps_lik(gps_trunc_alpha, obs_lons, obs_lats, obs_semi_majors,
                   obs_semi_minors, obs_orientations);
    DepthLik depth_lik(obs_depths);
    LocationDepthLik sattag_lik(gps_lik, depth_lik);


    /*
     * diffuse initial probability while aggregating likelihood
     */

    auto pred_step_iter = pred_steps.rbegin();
    auto pred_step_end = pred_steps.rend();

    // diffuse mass across likelihood steps
    for(unsigned int step_cur = nsteps; step_cur > 1; --step_cur) {

        // save prediction distribution
        if(pred_step_iter != pred_step_end) {
            if(*pred_step_iter == step_cur) {
                pvec.unswapActive();
                pred_distns.push_back(pvec.toNumericMatrix());
                pred_step_iter++;
                pvec.swapActive();
            }
        }

        // check end condition
        if(pred_step_iter == pred_step_end) {
            break;
        }

        // set observation likelihood, for diffusion
        sattag_lik.setLikToObs(step_cur - 1);

        Rcpp::checkUserInterrupt();

        // diffuse mass (i.e., update prediction distribution)
        backFilterMass<LocationDepthLik>(&pvec, &txmod, log_self_tx, &sattag_lik);
        // swap state
        pvec.swapActive();
    }

    return pred_distns;
}