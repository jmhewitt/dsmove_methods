//
// Created by Joshua Hewitt on 3/15/21.
//

#include <Rcpp.h>
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

        // labels for neighbors
        std::vector<Index> nbrs;

        // auto-regressive model component
        double beta_ar;
        bool beta_ar_nonzero;

    public:

        TxModel(Neighborhood &n, Directions &d) : nbhd(&n), dirs(&d),
            beta_ar(0.0) { };

        void constructProbs(const Index&, const Index&);
        void setBetaAR(double v) { beta_ar = v; };

        std::vector<double> logProbs();
        std::vector<double>* logProbs_ptr();

        std::vector<Index> neighbors();
        std::vector<Index>* neighbors_ptr();

        Index sampleNeighbor();

        // log density for a neighboring index
        double ld(const Index&);

};

template<typename N, typename D, typename I, typename size_type>
void TxModel<N,D,I,size_type>::constructProbs(const I &cur_loc,
                                              const I &prev_loc) {
    // set neighborhood structure
    nbhd->setCenter(cur_loc);
    dirs->setOrientation(cur_loc, prev_loc);
    size_type nnbrs = nbhd->neighborhoodSize();

    // initialize probability vector
    log_probs.clear();
    log_probs.reserve(nnbrs);

    // initialize storage for neighbor labels
    nbrs.clear();
    nbrs.reserve(nnbrs);

    // begin iterations
    nbrs.emplace_back(nbhd->nextNeighbor());

    // determine unnormalized probability of initial neighbor
    double lp = beta_ar * dirs->dotOrientation(nbhd->neighborHeading());
    log_probs.emplace_back(lp);
    log_mass = lp;

    // determine unnormalized probability of additional neighbors
    for(size_type i = 1; i < nnbrs; ++i) {
        nbrs.emplace_back(nbhd->nextNeighbor());
        double lp = beta_ar * dirs->dotOrientation(nbhd->neighborHeading());
        log_probs.emplace_back(lp);
        log_mass = log_add(log_mass, lp);
    }

    // standardize probabilities
    auto end = log_probs.end();
    for (auto it = log_probs.begin(); it != end; ++it)
        *it -= log_mass;
}

template<typename N, typename D, typename I, typename size_type>
std::vector<double> TxModel<N,D,I,size_type>::logProbs() {
    return log_probs;
}

template<typename N, typename D, typename I, typename size_type>
std::vector<double>* TxModel<N,D,I,size_type>::logProbs_ptr() {
    return &log_probs;
}

template<typename N, typename D, typename I, typename size_type>
std::vector<I> TxModel<N,D,I,size_type>::neighbors() {
    return nbrs;
}

template<typename N, typename D, typename I, typename size_type>
std::vector<I>* TxModel<N,D,I,size_type>::neighbors_ptr() {
    return &nbrs;
}

template<typename N, typename D, typename Index, typename size_type>
Index TxModel<N,D,Index,size_type>::sampleNeighbor() {

    // re-align neighborhood iterator with probability vector
    nbhd->resetNeighborIterator();

    size_type nnbrs = nbhd->neighborhoodSize();

    //
    // sample via cumulative CDF
    //

    double u = R::runif(0,1);
    auto prob_it = log_probs.begin();

    // first neighbor and cumulative probability
    Index res = nbhd->nextNeighbor();
    double p = std::exp(*prob_it);

    // aggregate mass of other neighbors
    for(size_type i = 1; i < nnbrs; ++i) {
        if(p >= u) break;
        res = nbhd->nextNeighbor();
        p += std::exp(*(++prob_it));
    }

    return res;
}

template<typename N, typename D, typename Index, typename size_type>
double TxModel<N,D,Index,size_type>::ld(const Index &dst) {

    auto nbr = nbrs.begin();
    auto end = log_probs.end();
    for(auto lp = log_probs.begin(); lp != end; ++lp) {
        if(dst == *(nbr++)) {
            return *lp;
        }
    }

    // return -Inf if index not found in neighborhood
    return -std::numeric_limits<double>::infinity();
}

#endif //DSMOVETOOLS_TXMODEL_H
