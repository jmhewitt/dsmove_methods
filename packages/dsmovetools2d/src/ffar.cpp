//
// Created by Joshua Hewitt on 7/14/21.
//

#include "ffar.h"
#include "CTDS2DProbs.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix FF_DTMC(
    std::vector<double> lons, std::vector<double> lats,
    std::vector<double> surface_heights, NumericMatrix init_dsts,
    NumericMatrix init_srcs, std::vector<double> init_log_probs,
    unsigned int steps, double log_self_tx, double betaAR
) {

    // initialize CTMC state space probability vector
    CTDS2DDomain pvec(lons, lats, surface_heights);

    // populate vector
    for(unsigned int ind = 0; ind < init_dsts.nrow(); ++ind) {
        pvec.set(init_srcs(ind, 0), init_srcs(ind, 1), init_dsts(ind, 0),
                 init_dsts(ind, 1), init_log_probs[ind]);
    }

    // finalize initial probability vector
    pvec.swapActive();

    class UniformLik {
    public:
        double ll(const int lon_to_ind, const int lat_to_ind) {
            return 0;
        }
    };

    UniformLik ulik;

    TxProbs txmod;
    txmod.setBetaAR(betaAR);

    // diffuse mass
    for(unsigned int step_cur = 0; step_cur < steps; ++step_cur) {
        // forward-filter all mass from the most recently diffused vector
        diffuseMass<UniformLik>(&pvec, &txmod, log_self_tx, &ulik);
        // swap state
        pvec.swapActive();
    }

    // extract diffused probabilities
    pvec.unswapActive();
    return pvec.toNumericMatrix();
}