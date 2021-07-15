//
// Created by Joshua Hewitt on 7/14/21.
//

#ifndef DSMOVETOOLS2D_FFAR_H
#define DSMOVETOOLS2D_FFAR_H

#include "CTDS2DDomain.h"
#include "CTDS2DProbs.h"
#include "log_complement.h"

template<typename SrcLik>
void diffuseMass(CTDS2DDomain *src, TxProbs *txmod,
                 double log_self_tx, SrcLik *srclik) {

    // compute probability of making a transition
    double log_tx = log_complement(log_self_tx);

    // log-mass of dst vector
    double dst_mass = -std::numeric_limits<double>::infinity();

    // forward-filter all mass from the src vector to the dst vector
    auto src_end = src->end();
    for(auto src_state = src->begin(); src_state != src_end; ++src_state) {
        double src_wt = src->logProbCached(*src_state);
        // diffuse from entries with finite mass
        if(std::isfinite(src_wt)) {
            // log-likelihood for current location
            double ll = srclik->ll(*src_state);
            // diffuse from entries with finite observation likelihood
            if(std::isfinite(ll)) {
                // add mass for self-transition to dst vector
                double m = src_wt + log_self_tx;
                src->add(*src_state, m);
                // aggregate total output mass
                dst_mass = log_add(dst_mass, m);
                // compute transition probabilities to neighbors
                txmod->constructProbs(*src_state);
                // weight stepwise tx. probs by non self-tx prob and atom's mass
                for(unsigned int i = 0; i < 4; ++i) {
                    if(src_state->nbr_to[i]) {
                        // mass of transitioning to dst state
                        double m = src_wt + log_tx + txmod->logProb(i);
                        src->add(*src_state->nbr_to[i], m);
                        // aggregate total output mass
                        dst_mass = log_add(dst_mass, m);
                    }
                }
            }
        }
    }

    // standardize distribution
    if(dst_mass != 0) {
        // invert mass so it can normalize diffused probabilities
        dst_mass *= -1;
        for(auto src_state = src->begin(); src_state != src_end; ++src_state) {
            src->scale(*src_state, dst_mass);
        }
    }
}

#endif //DSMOVETOOLS2D_FFAR_H
