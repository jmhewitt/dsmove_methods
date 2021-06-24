//
// Created by Joshua Hewitt on 12/14/20.
//

#include "SparseNdimArray.h"

#ifndef DSMOVETOOLS_FFAR_H
#define DSMOVETOOLS_FFAR_H

// let's use ints to index coordinates
typedef unsigned int IndexType;
// for specifying integer coordinates with grid dimension assigned at runtime
typedef std::vector<IndexType> VectorI;
// for specifying pairs of integer coordinates
typedef std::pair<VectorI, VectorI> VectorIPair;
// sparse spatial storage implemented via trees
typedef SparseNdimArrayLog<VectorIPair,
                           double,
                           std::map<VectorIPair, double>> LogARMap;

template <typename size_type, typename Neighborhood, typename TxMod,
        typename SrcLik>
void backFilterMassAF(LogARMap *src, LogARMap *dst, Neighborhood *nbhd,
                      TxMod *txmod, double log_self_tx, SrcLik *srclik) {
    // Parameters
    //  log_self_tx - log of self-transition probability

    // compute probability of making a transition
    double log_tx = log_complement(log_self_tx);

    // backward-filter all mass from the src vector to the dst vector
    auto src_mass_entry = src->data.begin();
    auto src_mass_end = src->data.end();
    for(src_mass_entry; src_mass_entry != src_mass_end; ++src_mass_entry) {
        // get weight for current entry
        double lw = src_mass_entry->second;
        // add self-transition probability
        double ll = srclik->ll(src_mass_entry->first);
        if(std::isfinite(ll)) {
            // aggregate back-filtered mass in destination
            dst->add(src_mass_entry->first, ll + log_self_tx + lw);
        }
        // loop over CTMC neighborhood of current entry "(cur_loc, prev_loc)"
        nbhd->setCenter(src_mass_entry->first.second);
        size_type nnbrs = nbhd->neighborhoodSize();
        for(size_type i = 0; i < nnbrs; ++i) {
            // get log-likelihood for observation
            VectorIPair ctmc_nbr = VectorIPair(
                    src_mass_entry->first.second, nbhd->nextNeighbor()
            );
            // compute and extract forward stepwise transition probabilities
            txmod->constructProbs(src_mass_entry->first.second,
                                  ctmc_nbr.second);
            double ltx = txmod->ld(src_mass_entry->first.first);
            double ll = srclik->ll(ctmc_nbr);
            if(std::isfinite(ll)) {
                // aggregate back-filtered mass in destination
                dst->add(ctmc_nbr, ll + log_tx + ltx + lw);
            }
        }
    }
}

template <typename size_type, typename Neighborhood, typename TxMod,
        typename SrcLik>
void diffuseMassSelfTxAR(LogARMap *src, LogARMap *dst, Neighborhood *nbhd,
                         TxMod *txmod, double log_self_tx, SrcLik *srclik) {
    // Parameters
    //  log_self_tx - log of self-transition probability

    // temporary storage for neighbor structures
    std::pair<VectorI, VectorI> ctmc_nbr;

    // compute probability of making a transition
    double log_tx = log_complement(log_self_tx);

    // log-mass of dst vector
    double dst_mass = -std::numeric_limits<double>::infinity();

    // forward-filter all mass from the src vector to the dst vector
    auto src_mass_entry = src->data.begin();
    auto src_mass_end = src->data.end();
    for(src_mass_entry; src_mass_entry != src_mass_end; ++src_mass_entry) {
        // get log-likelihood for current location
        double ll = srclik->ll(src_mass_entry->first);
        // forward-diffuse from location if observation likelihood is finite
        if(std::isfinite(ll)) {
            // find reachable nodes for current location
            nbhd->setCenter(src_mass_entry->first.first);
            size_type nnbrs = nbhd->neighborhoodSize();
            // add mass for self-transition to dst vector
            double m = src_mass_entry->second + log_self_tx;
            dst->add(src_mass_entry->first, m);
            // aggregate total output mass
            dst_mass = log_add(dst_mass, m);
            // compute and extract stepwise transition probabilities
            txmod->constructProbs(src_mass_entry->first.first,
                                  src_mass_entry->first.second);
            std::vector<VectorI> *nbrs = txmod->neighbors_ptr();
            std::vector<double> *lp = txmod->logProbs_ptr();
            // weight stepwise tx. probs. by non self-tx prob and mass of atom
            for(size_type i = 0; i < nnbrs; ++i) {
                // mass of transitioning to dst vector
                double m = src_mass_entry->second + (*lp)[i] + log_tx;
                ctmc_nbr.first = (*nbrs)[i];
                ctmc_nbr.second = src_mass_entry->first.first;
                dst->add(ctmc_nbr,m);
                // aggregate total output mass
                dst_mass = log_add(dst_mass, m);
            }
        }
    }
    // standardize distribution
    if(dst_mass != 0) {
        auto dst_mass_entry = dst->data.begin();
        auto dst_mass_end = dst->data.end();
        for(dst_mass_entry; dst_mass_entry != dst_mass_end; ++dst_mass_entry) {
            dst_mass_entry->second -= dst_mass;
        }
    }
}

#endif //DSMOVETOOLS_FFAR_H