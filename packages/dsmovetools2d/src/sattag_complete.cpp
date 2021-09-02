//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>
#include "CTDS2DProbs.h"

using namespace Rcpp;

// [[Rcpp::export]]
double completeObsLL(
    double betaAR, double speed, double cell_size,
    std::vector<double> durations,
    std::vector<unsigned int> direction_of_movement
) {

    // transition rate parameter (constant across space in 2 parameter model)
    double txscale = cell_size / speed;

    // initialize cell transition model (won't account for boundary effects)
    TxProbs txmod;
    txmod.setBetaAR(betaAR);

    // initialize likelihood
    double ll = 0;

    // aggregate likelihood of cell durations (no final duration observed)
    auto durations_iter = durations.begin();
    auto durations_end = (durations.end()) - 1;
    for(durations_iter; durations_iter != durations_end; ++durations_iter)
        ll += R::dexp(*durations_iter, txscale, true);

    // aggregate likelihood of transitions
    auto dom_iter = direction_of_movement.begin();
    auto dom_iter_end = direction_of_movement.end();
    // skip ahead to first known direction of movement
    dom_iter++;
    // update transition probabilities
    txmod.constructProbs(*dom_iter);
    // skip ahead to first fully observed transition
    for(++dom_iter; dom_iter != dom_iter_end; ++dom_iter) {
        // aggregate transition probability
        ll += txmod.logProb(*dom_iter);
        // update transition probabilities
        txmod.constructProbs(*dom_iter);
    }

    return ll;
}