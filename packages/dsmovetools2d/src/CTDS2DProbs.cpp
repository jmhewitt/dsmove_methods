//
// Created by Joshua Hewitt on 7/14/21.
//

#include "CTDS2DProbs.h"
#include "log_add.h"

using namespace Rcpp;

TxProbs::TxProbs() {

    // specify direction of movement for transition probability containers
    for(unsigned int i = 0; i < 4; ++i)
        tx_probs[i].direction_of_movement = i;

    // set up linked list to allow cyclic access to probability containers
    tx_probs[0].next = &tx_probs[1];
    tx_probs[1].next = &tx_probs[2];
    tx_probs[2].next = &tx_probs[3];
    tx_probs[3].next = &tx_probs[0];

    // initialize transition parameters
    for(unsigned int i = 0; i < 4; ++i)
        prob_order[i] = 0;
}

void TxProbs::setBetaAR(double betaAR) {
    // update non-static transition parameters
    prob_order[0] = betaAR;
    prob_order[2] = -betaAR;
}

void TxProbs::constructProbs(const CTDS2DState &state) {

    // reset total transition mass
    log_total_mass = -std::numeric_limits<double>::infinity();

    // assign first probability to the transition that continues movement dir.
    ProbContainer *cur = &tx_probs[state.direction_of_movement];

    // cycle through all 4 neighbors CW from initial position
    for(unsigned int i = 0; i < 4; ++i) {
        if(state.nbr_to[cur->direction_of_movement]) {
            cur->log_prob = prob_order[i];
            log_total_mass = log_add(log_total_mass, prob_order[i]);
        } else {
            cur->log_prob = -std::numeric_limits<double>::infinity();
        }
        // increment direction
        cur = cur->next;
    }

}

double TxProbs::logProb(unsigned int direction_of_movement) {
    return tx_probs[direction_of_movement].log_prob - log_total_mass;
}

// [[Rcpp::export]]
NumericMatrix LogTxProbs(
        std::vector<double> lons, std::vector<double> lats,
        std::vector<double> surface_heights, int lon_from_ind, int lat_from_ind,
        int lon_to_ind, int lat_to_ind, double betaAR
) {

    // initialize CTMC state space
    CTDS2DDomain domain(lons, lats, surface_heights);

    // access state from which transitions are made
    CTDS2DState *tgt = domain.statePtr(
        lon_from_ind, lat_from_ind, lon_to_ind, lat_to_ind
    );

    // compute transition probabilities for state
    TxProbs probs;
    probs.setBetaAR(betaAR);
    probs.constructProbs(*tgt);

    // extract transition probabilities
    for(unsigned int i = 0; i < 4; ++i) {
        if(tgt->nbr_to[i]) {
            domain.set(*tgt->nbr_to[i], probs.logProb(i));
        }
    }

    return domain.toNumericMatrix();
}

// [[Rcpp::export]]
NumericMatrix LogTxProbsElevation(
        std::vector<double> lons, std::vector<double> lats,
        std::vector<double> surface_heights, int lon_from_ind, int lat_from_ind,
        int lon_to_ind, int lat_to_ind, double betaAR, double min_elevation,
        double max_elevation
) {

    // initialize CTMC state space
    CTDS2DDomain domain(lons, lats, surface_heights);

    // filter state space
    CTDS2DStateElevationFilter height_filter(min_elevation, max_elevation);
    domain.filterStates(height_filter);

    // access state from which transitions are made
    CTDS2DState *tgt = domain.statePtr(
            lon_from_ind, lat_from_ind, lon_to_ind, lat_to_ind
    );

    // compute transition probabilities for state
    TxProbs probs;
    probs.setBetaAR(betaAR);
    probs.constructProbs(*tgt);

    // extract transition probabilities
    for(unsigned int i = 0; i < 4; ++i) {
        if(tgt->nbr_to[i]) {
            domain.set(*tgt->nbr_to[i], probs.logProb(i));
        }
    }

    return domain.toNumericMatrix();
}