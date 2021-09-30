//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>
#include "CTDS2DDomain.h"
#include "CTDS2DProbs.h"
#include "ffar.h"
#include "ExactLocationLik.h"

using namespace Rcpp;

// [[Rcpp::export]]
double ExactLocFilteredLL(
    NumericMatrix init_dsts, NumericMatrix init_srcs,
    std::vector<double> init_log_probs,
    std::vector<double> obs_x_coords, std::vector<double> obs_y_coords,
    std::vector<double> x_coords, std::vector<double> y_coords,
    std::vector<double> surface_heights,
    double log_self_tx, double betaAR, double lptrunc
) {

    // number of steps to integrate over in likelihood
    unsigned int nsteps = obs_x_coords.size();

    // initialize log-likelihood
    double ll = 0;

    // initialize CTMC state space probability vector
    CTDS2DDomain pvec(x_coords, y_coords, surface_heights);

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


    // initialize likelihood
    ExactLocationLik loc_lik(obs_x_coords, obs_y_coords);



    /*
     * diffuse initial probability while aggregating likelihood
     */

    double lpmax = 0;
    double lpthresh = lpmax - lptrunc;

    bool firstError = true;

    unsigned int checkmark_interval = .1 * nsteps;

    // diffuse mass across likelihood steps
    for(unsigned int step_cur = 0; step_cur < nsteps; ++step_cur) {

        loc_lik.setLikToObs(step_cur);

        // aggregate marginal likelihood over last diffused probability state
        bool finiteStateWeight = false;
        bool finiteLikObs = false;
        double ll_step = -std::numeric_limits<double>::infinity();
//        auto state_end = pvec.end();
//        for(auto state_it = pvec.begin(); state_it != state_end; ++state_it) {
//        auto state_end = pvec.states_written.end();
//        for(auto state_it = pvec.states_written.begin(); state_it != state_end;
//            ++state_it) {
        auto state_end = pvec.end();
        for(auto state_it = pvec.begin(); state_it != state_end; ++state_it) {
            // extract state weight
            double w = pvec.logProbCached(*state_it);
            if(std::isfinite(w)) {
                finiteStateWeight = true;
                // get log-likelihood for diffused location
                double ll_state = loc_lik.ll(*state_it);
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
                NumericMatrix state = pvec.toNumericMatrix(false);
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

//        Rcpp::checkUserInterrupt();

        // diffuse mass (i.e., update prediction distribution)
        lpmax = diffuseMassPred<ExactLocationLik>(&pvec, &txmod, log_self_tx,
                                          &loc_lik, lpthresh);

        lpthresh = lpmax - lptrunc;

        // swap state
        pvec.swapActive();
    }

    return ll;
}
