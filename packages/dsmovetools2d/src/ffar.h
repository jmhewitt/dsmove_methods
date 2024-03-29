//
// Created by Joshua Hewitt on 7/14/21.
//

#ifndef DSMOVETOOLS2D_FFAR_H
#define DSMOVETOOLS2D_FFAR_H

#include "CTDS2DDomain.h"
#include "CTDS2DProbs.h"
#include "log_complement.h"

template<typename SrcLik>
void backFilterMass(CTDS2DDomain *src, TxProbs *txmod, double log_self_tx,
                    SrcLik *srclik) {

    // compute probability of making a transition
    double log_tx = log_complement(log_self_tx);

    // backward-filter all mass from the src vector to the dst vector
    auto src_end = src->end();
    for(auto src_state = src->begin(); src_state != src_end; ++src_state) {
        double src_wt = src->logProbCached(*src_state);
        // (backward) diffuse to entries with finite mass
        if(std::isfinite(src_wt)) {
            // add self-transition probability
            double ll = srclik->ll(*src_state);
            if(std::isfinite(ll)) {
                // aggregate back-filtered mass in destination
                src->add(*src_state, ll + log_self_tx + src_wt);
            }
            // loop over locations that transition to src_state
            for(unsigned int i = 0; i < 4; ++i) {
                if(src_state->nbr_from[i]) {
                    // log-likelihood for state
                    ll = srclik->ll(*(src_state->nbr_from[i]));
                    if(std::isfinite(ll)) {
                        // compute transition probabilities from state
                        txmod->constructProbs(*(src_state->nbr_from[i]));
                        // aggregate back-filtered mass in dst vector
                        double m = ll + log_tx + txmod->logProb(i) + src_wt;
                        src->add(*(src_state->nbr_from[i]), m);
                    }
                }
            }
        }
    }
}

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

    src_end = src->end(true);

    // standardize distribution
    if(dst_mass != 0) {
        // invert mass so it can normalize diffused probabilities
        dst_mass *= -1;
        for(auto src_state = src->begin(true); src_state != src_end; ++src_state) {
            src->scale(*src_state, dst_mass);
        }
    }
}

template<typename SrcLik>
double diffuseMassPred(CTDS2DDomain *src, TxProbs *txmod,
                 double log_self_tx, SrcLik *srclik,
                 double lpthresh = -std::numeric_limits<double>::infinity()) {

    // compute probability of making a transition
    double log_tx = log_complement(log_self_tx);

    // log-mass of dst vector
    double dst_mass = -std::numeric_limits<double>::infinity();

    // forward-filter all mass from the src vector to the dst vector
    auto src_end = src->end();
    for(auto src_state = src->begin(); src_state != src_end; ++src_state) {
//    auto src_end = src->states_written.end();
//    for(auto src_state = src->states_written.begin();
//    src_state != src_end; ++src_state) {
        double src_wt = src->logProbCached(*src_state);
        // diffuse from entries with finite mass
//        Rcpp::Rcout << "src_wt: " << src_wt << std::endl;
//        Rcpp::Rcout << "lpthresh: " << lpthresh << std::endl;
        if(std::isfinite(src_wt)) {
        if(src_wt > lpthresh) {
            // log-likelihood for current location
            double ll = srclik->ll(*src_state);
//            Rcpp::Rcout << "ll: " << ll << std::endl;
            // diffuse from entries with finite observation likelihood
            if(std::isfinite(ll)) {
                // add mass for self-transition to dst vector
                double m = src_wt + ll + log_self_tx;
                src->add(*src_state, m);
//                src->getNonzeroActive()->insert(*src_state);
                // aggregate total output mass
                dst_mass = log_add(dst_mass, m);
                // compute transition probabilities to neighbors
                txmod->constructProbs(*src_state);
                // weight stepwise tx. probs by non self-tx prob and atom's mass
                for(unsigned int i = 0; i < 4; ++i) {
                    if((src_state)->nbr_to[i]) {
                        // mass of transitioning to dst state
                        double m = src_wt + ll + log_tx + txmod->logProb(i);
                        src->add(*((src_state)->nbr_to[i]), m);
//                        src->getNonzeroActive()->insert((*src_state)->nbr_to[i]);
                        // aggregate total output mass
                        dst_mass = log_add(dst_mass, m);
                    }
                }
            }
        }}
    }

    double lpmax = -std::numeric_limits<double>::infinity();

//    src_end = src->getNonzeroActive()->end();

    src_end = src->end(true);

    // standardize distribution
    if(dst_mass != 0) {
        // invert mass so it can normalize diffused probabilities
        dst_mass *= -1;
        for(auto src_state = src->begin(true); src_state != src_end; ++src_state) {
//        for(auto src_state = src->states_written.begin();
//            src_state != src_end; ++src_state) {
            src->scale(*src_state, dst_mass);
            double m = src->logProb(*src_state);
            if(m > lpmax) {
                lpmax = m;
            }
        }
    }

    return lpmax;
}

#endif //DSMOVETOOLS2D_FFAR_H
