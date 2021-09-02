//
// Created by Joshua Hewitt on 7/14/21.
//

#ifndef DSMOVETOOLS2D_CTDS2DPROBS_H
#define DSMOVETOOLS2D_CTDS2DPROBS_H

#include "CTDS2DDomain.h"

class TxProbs {

private:

    // linked list element for probability container
    struct ProbContainer {
        ProbContainer *next;
        double log_prob;
        unsigned int direction_of_movement;
    };

    // probabilities for transitions to N,E,S,W neighbors
    double prob_order[4];
    ProbContainer tx_probs[4];

    // total mass available for transitions
    double log_total_mass;

public:

    // initialize probability container elements
    TxProbs();

    // set directional persistence parameter
    void setBetaAR(double betaAR);

    // compute transition probabilities for a CTMC state
    void constructProbs(const CTDS2DState &state);

    // compute transition probabilities for a known direction of movement
    void constructProbs(unsigned int direction_of_movement);

    // extract transition probability to neighboring N,E,S,W state
    double logProb(unsigned int direction_of_movement);

};

#endif //DSMOVETOOLS2D_CTDS2DPROBS_H
