//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>
#include "ffrw.h"
#include "SparseNdimArray.h"
#include "ZConstrainedRookNeighborhood.h"
#include "log_add.h"

using namespace Rcpp;

// let's use ints to index coordinates
typedef unsigned int IndexType;
// for specifying integer coordinates with grid dimension assigned at runtime
typedef std::vector<IndexType> VectorI;
// sparse spatial storage implemented via trees
typedef SparseNdimArrayLog<VectorI, double, std::map<VectorI, double>> LogArrayMap;
// rook neighborhoods on integer grids
typedef ZConstrainedRookNeighborhood<IndexType, VectorI> ZRN;


template<typename size_type, typename Index, typename Neighborhood>
std::vector<Index> bsrw_dst(
    const Index &src, const Index &dst, const Index &dims,
    Neighborhood &nbhd, Neighborhood &nbhd_cpy, const size_type steps,
    const size_type max_steps
) {

    // initial diffusion probability container: all mass at src
    LogArrayMap a0;
    a0.set(src, 0);

    // forward filter to destination
    ReachableProbs<size_type, LogArrayMap> a_diffused =
        ffrw_dst_reachable<size_type, Index, Neighborhood, LogArrayMap>(
            dims, a0, steps, nbhd, dst, max_steps
    );

    // extract actual path length from diffusion
    size_type path_len = a_diffused.probs.size();

    // initialize bridged random path
    std::vector<Index> path;
    path.reserve(path_len);
    path.emplace_back(dst);

    // backward sample
    auto diffusion = a_diffused.probs.rbegin() + 1;
    auto diffusion_end = a_diffused.probs.rend();
    for(diffusion; diffusion != diffusion_end; ++diffusion) {

        // last-sampled location in path can be reached by its neighbors
        nbhd.setCenter(path.back());
        size_type nnbrs = nbhd.neighborhoodSize();

        // initialize sampling probabilities for neighbors
        LogArrayMap nbr_log_probs;
        double log_mass;
        bool finite_log_mass = false;

        // determine unnormalized sampling probabilities for neighbors
        for(size_type i = 0; i < nnbrs; ++i) {
            Index neighbor = nbhd.nextNeighbor();
            // only connect to nbrs that can be reached in current diffusion
            if(diffusion->notNull(neighbor)) {
                // weight filtering dist'n. by neighbor's RW forward probability
                nbhd_cpy.setCenter(neighbor);
                double lp = diffusion->data[neighbor] -
                        std::log(nbhd_cpy.neighborhoodSize());
                nbr_log_probs.data[neighbor] = lp;
                // aggregate sampling distribution's log mass
                if(finite_log_mass) {
                    log_mass = log_add(log_mass, lp);
                } else {
                    log_mass = lp;
                    finite_log_mass = true;
                }
            }
        }

        // standardize probabilities
        auto end = nbr_log_probs.data.end();
        for (auto it = nbr_log_probs.data.begin(); it != end; ++it)
            it->second -= log_mass;

        //
        // sample neighbor via cumulative CDF
        //

        double u = R::runif(0,1);
        auto prob_it = nbr_log_probs.data.begin();

        // first neighbor's cumulative probability
        double p = std::exp(prob_it->second);

        // aggregate mass of other neighbors
        for(size_type i = 1; i < nnbrs; ++i) {
            if(p >= u) break;
            p += std::exp((++prob_it)->second);
        }

        // save sampled neighbor
        path.emplace_back(prob_it->first);

    }

    // correct path order, then return sample
    std::reverse(path.begin(), path.end());
    return path;
};


// [[Rcpp::export]]
NumericMatrix SampleConstrainedBridgedRWPath(
    std::vector<unsigned int> a0coords, std::vector<unsigned int> dstcoords,
    std::vector<unsigned int> dims, unsigned int steps, unsigned int max_steps,
    std::vector<double> surface_heights, std::vector<double> domain_heights
) {

    // munge source and destination coordinates into array's storage format
    VectorI src_index(dims.size());
    VectorI dst_index(dims.size());
    for(IndexType j=0; j < dims.size(); ++j) {
        src_index[j] = a0coords[j];
        dst_index[j] = dstcoords[j];
    }

    // construct neighborhood
    ZRN zrn(dims, surface_heights.data(), domain_heights.data());
    ZRN zrn2(dims, surface_heights.data(), domain_heights.data());

    // sample path
    std::vector<VectorI> path = bsrw_dst<IndexType , VectorI, ZRN>(
            src_index, dst_index, dims, zrn, zrn2, steps, max_steps
    );

    // initialize container to store extracted sampled path
    NumericMatrix out = NumericMatrix(path.size(), dims.size());

    // extract path: loop over path steps
    for(unsigned int step = 0; step < path.size(); ++step) {
        // extract each coordinate for step
        for(unsigned int i = 0; i < dims.size(); ++i) {
            out(step,i) = path[step][i];
        }
    }

    return out;
}