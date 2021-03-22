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


template<typename Index>
struct BSRWFamily {
    std::vector<std::vector<Index>> path;
    std::vector<double> log_weight;
};

/**
 * Forward-filter backwards-sample a family of n paths from src to dst.  Path
 * lengths are uniformly sampled from among the collection of viable path
 * lengths.
 *
 * @tparam size_type
 * @tparam Index
 * @tparam Neighborhood
 * @param src
 * @param dst
 * @param dims
 * @param nbhd
 * @param nbhd_cpy
 * @param steps
 * @param max_steps
 * @param n
 * @return
 */
template<typename size_type, typename Index, typename Neighborhood>
BSRWFamily<Index> bsrw_dst_family(
        const Index &src, const Index &dst, const Index &dims,
        Neighborhood &nbhd, Neighborhood &nbhd_cpy, const size_type steps,
        const size_type max_steps, const size_type n
) {

    // initial diffusion probability container: all mass at src
    LogArrayMap a0;
    a0.set(src, 0);

    // forward filter to destination
    ReachableProbs<size_type, LogArrayMap> a_diffused =
            ffrw_dst_reachable<size_type, Index, Neighborhood, LogArrayMap>(
                    dims, a0, steps, nbhd, dst, max_steps
            );

    // initialize container for sampled paths
    BSRWFamily<Index> path_fam;
    path_fam.path.reserve(n);
    path_fam.log_weight.reserve(n);

    double base_weight = - std::log(a_diffused.reachable.size());

    for(size_type path_id = 0; path_id < n; ++path_id) {

        // sample a viable path length from diffusion
        size_type path_len = a_diffused.reachable[
            std::floor(R::runif(0,a_diffused.reachable.size()))
        ];

        // initialize path likelihood with likelihood for path length
        double log_weight = base_weight;

        // initialize bridged random path
        std::vector<Index> path;
        path.reserve(path_len);
        path.emplace_back(dst);

        // backward sample
        auto diffusion = a_diffused.probs.rbegin() +
            (a_diffused.probs.size() - path_len);
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
                    // weight filtering dist'n. by neighbor's RW forward prob.
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

            // aggregate sampling weight
            log_weight += prob_it->second;

        }

        // correct path order, then export the sample
        std::reverse(path.begin(), path.end());
        path_fam.path.emplace_back(path);
        path_fam.log_weight.emplace_back(log_weight);
    }

    return path_fam;
};


/**
 * Forward-filter backwards-sample a path of at least steps length long from
 * src to dst.
 *
 * @tparam size_type
 * @tparam Index
 * @tparam Neighborhood
 * @param src
 * @param dst
 * @param dims
 * @param nbhd
 * @param nbhd_cpy
 * @param steps
 * @param max_steps
 * @return
 */
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

// [[Rcpp::export]]
List SampleConstrainedBridgedRWPathFamily(
    std::vector<unsigned int> a0coords, std::vector<unsigned int> dstcoords,
    std::vector<unsigned int> dims, unsigned int steps, unsigned int max_steps,
    std::vector<double> surface_heights, std::vector<double> domain_heights,
    unsigned int n
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

    // sample family of paths
    BSRWFamily<VectorI> path_fam_raw =
        bsrw_dst_family<IndexType , VectorI, ZRN>(
            src_index, dst_index, dims, zrn, zrn2, steps, max_steps, n
    );

    // initialize container to store family of sampled paths
    Rcpp:List path_fam(n);

    std::vector<std::vector<VectorI>> path = path_fam_raw.path;
    for(unsigned int path_ind = 0; path_ind < n; ++path_ind) {
        // initialize container to store extracted sampled path
        NumericMatrix out = NumericMatrix(path[path_ind].size(), dims.size());
        // extract path: loop over path steps
        for(unsigned int step = 0; step < path[path_ind].size(); ++step) {
            // extract each coordinate for step
            for(unsigned int i = 0; i < dims.size(); ++i) {
                out(step,i) = path[path_ind][step][i];
            }
        }
        // save extracted path
        path_fam[path_ind] = out;
    }

    return List::create(Named("path") = path_fam,
                        Named("log_weights") = path_fam_raw.log_weight);
}