//
// Created by Joshua Hewitt on 12/14/20.
//

#ifndef DSMOVETOOLS_FFRW_H
#define DSMOVETOOLS_FFRW_H

 template<typename size_type, typename StorageArray>
 struct ReachableProbs {
     std::vector<StorageArray> probs;
     std::vector<size_type> reachable;
 };

/**
* Forward filtering a random walk along a grid, while storing all intermediate
* distributions.  Forward filtering continues until a destination location is
* reached and at least "steps" numbers of diffusions are conducted.
*
* @param dims specify the number of locations along each dimension in grid
* @param a0 initial probability mass vector
* @param steps minimum number of forward-diffusion steps to take
* @param max_steps maximum number of steps to try before aborting diffusion
* @param dst location to reach via forward filtering
* @return (sparse) diffused mass vectors
*/
template<typename size_type, typename Index, typename Neighborhood, typename StorageArray>
ReachableProbs<size_type, StorageArray> ffrw_dst_reachable(
    const Index &dims, const StorageArray &a0, const size_type steps,
    Neighborhood &nbhd, const Index &dst, const size_type max_steps
);

#endif // DSMOVETOOLS_FFRW_H