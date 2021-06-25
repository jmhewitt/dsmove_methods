//
// Created by Joshua Hewitt on 12/14/20.
//

#include "CachedSparseNdimArray.h"

#ifndef DSMOVETOOLS_FFAR_CACHED_H
#define DSMOVETOOLS_FFAR_CACHED_H


template <typename size_type, typename Neighborhood, typename TxMod,
        typename SrcLik>
void diffuseMassSelfTxARDual(DualSparseCoordVec *probArray, Neighborhood *nbhd,
                             TxMod *txmod, double log_self_tx, SrcLik *srclik) {
    // Parameters
    //  log_self_tx - log of self-transition probability

    // temporary storage for neighbor structures
    DualSparseCoordVec::CoordPair ctmc_nbr;

    // compute probability of making a transition
    double log_tx = log_complement(log_self_tx);

    // log-mass of dst vector
    double dst_mass = -std::numeric_limits<double>::infinity();

    // forward-filter all mass from the src vector to the dst vector
    auto src_mass_entry = probArray->data.begin();
    auto src_mass_end = probArray->data.end();
    for(src_mass_entry; src_mass_entry != src_mass_end; ++src_mass_entry) {
        double src_wt = probArray->inactiveValue(src_mass_entry->second);
        if(std::isfinite(src_wt)) {
            // get log-likelihood for current location
            double ll = srclik->ll(src_mass_entry->first);
            // forward-diffuse from location if observation likelihood is finite
            if(std::isfinite(ll)) {
                // find reachable nodes for current location
                nbhd->setCenter(src_mass_entry->first.first);
                size_type nnbrs = nbhd->neighborhoodSize();
                // add mass for self-transition to dst vector
                double m = src_wt + log_self_tx;
                probArray->addToActive(src_mass_entry->first, m);
                // aggregate total output mass
                dst_mass = log_add(dst_mass, m);
                // compute and extract stepwise transition probabilities
                txmod->constructProbs(src_mass_entry->first.first,
                                      src_mass_entry->first.second);
                std::vector<DualSparseCoordVec::Coord> *nbrs = txmod->neighbors_ptr();
                std::vector<double> *lp = txmod->logProbs_ptr();
                // weight stepwise tx. probs. by non self-tx prob and mass of atom
                for(size_type i = 0; i < nnbrs; ++i) {
                    // mass of transitioning to dst vector
                    double m = src_wt + (*lp)[i] + log_tx;
                    ctmc_nbr.first = (*nbrs)[i];
                    ctmc_nbr.second = src_mass_entry->first.first;
                    probArray->addToActive(ctmc_nbr,m);
                    // aggregate total output mass
                    dst_mass = log_add(dst_mass, m);
                }
            }
        }
    }
    // standardize distribution
    if(dst_mass != 0) {
        auto dst_mass_entry = probArray->data.begin();
        auto dst_mass_end = probArray->data.end();
        for(dst_mass_entry; dst_mass_entry != dst_mass_end; ++dst_mass_entry) {
            probArray->unprotectedSubtractFromActive(dst_mass_entry->second,
                                                     dst_mass);
        }
    }
}

#endif //DSMOVETOOLS_FFAR_CACHED_H