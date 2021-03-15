//
// Created by Joshua Hewitt on 3/15/21.
//

#include "Rcpp.h"
#include "log_add.h"
#include "RookHeading.h"

#ifndef DSMOVETOOLS_TXMODEL_H
#define DSMOVETOOLS_TXMODEL_H

template<typename Neighborhood, typename Directions, typename Index,
        typename size_type>
class TxModel {

    private:

        Neighborhood *nbhd;
        Directions *dirs;

        // unstandardized log_probability for each neighbor
        std::vector<double> log_probs;
        // log of normalizing constant
        double log_mass;

        // auto-regressive model component
        double beta_ar;
        bool beta_ar_nonzero;

    public:

        TxModel(Neighborhood &n, Directions &d) : nbhd(&n), dirs(&d),
            beta_ar(0.0) { };

        void constructProbs(const Index&, const Index&);
        void setBetaAR(double v) { beta_ar = v; };

        std::vector<double> logProbs();


};

template<typename N, typename D, typename I, typename size_type>
void TxModel<N,D,I,size_type>::constructProbs(const I &cur_loc,
                                              const I & prev_loc) {
    // set neighborhood structure
    nbhd->setCenter(cur_loc);
    dirs->setOrientation(cur_loc, prev_loc);
    size_type nnbrs = nbhd->neighborhoodSize();

    // initialize probability vector
    log_probs.clear();
    log_probs.reserve(nnbrs);
    nbhd->nextNeighbor();

    // determine unnormalized probability of initial neighbor
    double lp = beta_ar * dirs->dotOrientation(nbhd->neighborHeading());
    log_probs.emplace_back(lp);
    log_mass = lp;

    // determine unnormalized probability of additional neighbors
    for(size_type i = 1; i < nnbrs; ++i) {
        nbhd->nextNeighbor();
        double lp = beta_ar * dirs->dotOrientation(nbhd->neighborHeading());
        log_probs.emplace_back(lp);
        log_mass = log_add(log_mass, lp);
    }
}

template<typename N, typename D, typename I, typename size_type>
std::vector<double> TxModel<N,D,I,size_type>::logProbs() {
    // initialize return
    std::vector<double> out(log_probs);
    // shift log probs by normalizing constant
    auto end = out.end();
    for (auto it = out.begin(); it != end; ++it)
        *it -= log_mass;
    // return normalized probabilities on log scale
    return out;
}

#endif //DSMOVETOOLS_TXMODEL_H
