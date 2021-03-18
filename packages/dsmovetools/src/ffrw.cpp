//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>
#include "SparseNdimArray.h"
#include "RookNeighborhood.h"
#include "ZConstrainedRookNeighborhood.h"

using namespace Rcpp;

// let's use ints to index coordinates
typedef unsigned int IndexType;
// for specifying integer coordinates with grid dimension assigned at runtime
typedef std::vector<IndexType> VectorI;
// sparse spatial storage implemented via trees
typedef SparseNdimArray<VectorI, double, std::map<VectorI, double>> ArrayMap;
typedef SparseNdimArrayLog<VectorI, double, std::map<VectorI, double>> LogArrayMap;
// rook neighborhoods on integer grids
typedef RookNeighborhood<IndexType, VectorI> RN;
typedef ZConstrainedRookNeighborhood<IndexType, VectorI> ZRN;


template <typename I, typename V, typename S, typename size_type, 
          typename Neighborhood>
void diffuseMass(SparseNdimArrayBase<I,V,S> *src,
                 SparseNdimArrayBase<I,V,S> *dst,
                 Neighborhood *nbhd) {
    // forward-filter all mass from the src vector to the dst vector
    auto src_mass_entry = src->data.begin();
    auto src_mass_end = src->data.end();
    for(src_mass_entry; src_mass_entry != src_mass_end; ++src_mass_entry) {
        // find neighborhood for previous location
        nbhd->setCenter(src_mass_entry->first);
        size_type nnbrs = nbhd->neighborhoodSize();
        // diffuse mass, following a random walk along neighbors
        double scaledNbrs = dst->normalizeScale(nnbrs);
        for(size_type i = 0; i < nnbrs; ++i) {
            dst->addScaled(nbhd->nextNeighbor(), src_mass_entry->second, 
                           scaledNbrs);
        }
    }
}

/**
 * Forward filtering a random walk along a grid, without storing all
 * intermediate distributions.
 *
 * @param dims specify the number of locations along each dimension in grid
 * @param a0 initial probability mass vector
 * @param steps number of forward-diffusion steps to take
 * @param nbhd Class that defines the neighborhood for arbitrary locations
 * @return (sparse) diffused mass vectors
 */
template<typename size_type, typename Neighborhood, typename StorageArray>
StorageArray ffrw_light(const VectorI &dims, const StorageArray &a0,
                        const IndexType steps, Neighborhood &nbhd) {

    // initialize forward filtering vectors and initial mass
    StorageArray cur, prev;
    prev = a0;

    // diffuse mass
    for(IndexType step_cur = 0; step_cur < steps; ++step_cur) {
        // forward-filter all mass from the most recently diffused vector
        diffuseMass<VectorI, double, std::map<VectorI, double>, size_type, 
                   Neighborhood>(&prev, &cur, &nbhd);
        // swap state
        prev.data.swap(cur.data);
        cur.data.clear();
    }

    return prev;
};

/**
 * Forward filtering a random walk along a grid, while storing all intermediate
 * distributions.
 *
 * @param dims specify the number of locations along each dimension in grid
 * @param a0 initial probability mass vector
 * @param steps number of forward-diffusion steps to take
 * @return (sparse) diffused mass vectors
 */
 template<typename size_type, typename Neighborhood, typename StorageArray>
std::vector<StorageArray> ffrw(const VectorI &dims, const StorageArray &a0,
                               const IndexType steps, Neighborhood &nbhd) {

    // initialize forward filtering vectors and initial mass
    std::vector<StorageArray> ffprob(steps + 1);
    ffprob[0] = a0;

    // diffuse mass
    auto step_cur = ffprob.begin();
    auto step_prev = step_cur;
    auto step_end = ffprob.end();
    for(++step_cur; step_cur != step_end;) {
        // forward-filter all mass from the most recently diffused vector
        diffuseMass<VectorI, double, std::map<VectorI, double>, size_type, 
                   Neighborhood>(&(*step_prev), &(*step_cur), &nbhd);
        // increment iterators
        ++step_cur;
        ++step_prev;
    }

    return ffprob;
};

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
template<typename size_type, typename Neighborhood, typename StorageArray>
ReachableProbs<IndexType, StorageArray> ffrw_dst_reachable(
    const VectorI &dims, const StorageArray &a0, const IndexType steps,
    Neighborhood &nbhd, const VectorI &dst, const IndexType max_steps
) {

    // initialize forward filtering vectors and initial mass
    ReachableProbs<IndexType, StorageArray> ffprob;
    ffprob.probs.reserve(steps + 1);
    ffprob.probs.emplace_back(a0);

    // initialize storage of steps at which dst is reachable
    ffprob.reachable.reserve(steps + 1);

    // diffuse mass until number of steps and mass at dst is attained
    bool dst_reached = false;
    for(IndexType i = 1; i < max_steps; ++i) {
        // initialize next diffused probability vector
        ffprob.probs.emplace_back();
        StorageArray *cur = &ffprob.probs.back();
        // last diffused mass vector
        //   NOTE: Can't get pointer sooner since emplacement may move all data
        //         if vector size needs to be expanded.
        StorageArray *prev = cur - 1;
        // forward-filter all mass from the most recently diffused vector
        diffuseMass<VectorI, double, std::map<VectorI, double>, size_type,
                Neighborhood>(prev, cur, &nbhd);
        // updates if destination is reachable
        if(cur->notNull(dst)) {
            ffprob.reachable.emplace_back(i);
            dst_reached = true;
        }
        // check termination conditions
        if((i >= steps) && dst_reached) {
            break;
        }
    }

    ffprob.reachable.shrink_to_fit();

    return ffprob;
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
 template<typename size_type, typename Neighborhood, typename StorageArray>
std::vector<StorageArray> ffrw_dst(const VectorI &dims, const StorageArray &a0,
                                   const IndexType steps, Neighborhood &nbhd,
                                   const VectorI &dst, 
                                   const IndexType max_steps) {

    ReachableProbs<IndexType, StorageArray> res =
            ffrw_dst_reachable<size_type, Neighborhood, StorageArray>(
                dims, a0, steps, nbhd, dst, max_steps
            );

    return res.probs;
};


// [[Rcpp::export]]
NumericMatrix TestFFRW(NumericMatrix a0coords, NumericVector a0values,
                       std::vector<unsigned int> dims, int steps) {
    // Parameters:
    //   dims - number of locations along each dimension

    // initial probability container
     ArrayMap a0;

    // fill initial probability container
    for(int i=0; i < a0coords.nrow(); i++) {
        // munge R coordinate into array's storage format
        VectorI c(a0coords.ncol());
        for(int j=0; j< a0coords.ncol(); ++j)
            c[j] = a0coords(i,j);
        // insert coord/value pair into sparse array
        a0.set(c, a0values(i));
    }

    // diffuse initial probability
    RN rn(dims);
    std::vector<ArrayMap> ffprobs = ffrw<IndexType, RN, ArrayMap>(dims, a0, 
                                                                  steps, rn);

    // extract final, diffused probability
    ArrayMap af = ffprobs.back();
    NumericMatrix out = NumericMatrix(af.data.size(), a0coords.ncol() + 1);
    int i=0;
    for(auto iter = af.data.begin(); iter != af.data.end(); ++iter) {
        // extract coordinate info from array's storage format
        VectorI c = iter->first;
        for(int j=0; j < a0coords.ncol(); ++j) {
            out(i,j) = c[j];
        }
        // extract value information
        out(i++, out.ncol() - 1) = iter->second;
    }

    // return lexicographically sorted coord/val pairs
    return out;
}

// [[Rcpp::export]]
NumericMatrix TestFFRWLight(NumericMatrix a0coords, NumericVector a0values,
                            std::vector<unsigned int> dims, int steps) {
    // Parameters:
    //   dims - number of locations along each dimension

    // initial probability container
    ArrayMap a0;

    // fill initial probability container
    for(int i=0; i < a0coords.nrow(); i++) {
        // munge R coordinate into array's storage format
        VectorI c(a0coords.ncol());
        for(int j=0; j< a0coords.ncol(); ++j)
            c[j] = a0coords(i,j);
        // insert coord/value pair into sparse array
        a0.set(c, a0values(i));
    }

    // diffuse initial probability
    RN rn(dims);
    ArrayMap af = ffrw_light<IndexType, RN, ArrayMap>(dims, a0, steps, rn);

    // extract final, diffused probability
    NumericMatrix out = NumericMatrix(af.data.size(), a0coords.ncol() + 1);
    int i=0;
    for(auto iter = af.data.begin(); iter != af.data.end(); ++iter) {
        // extract coordinate info from array's storage format
        VectorI c = iter->first;
        for(int j=0; j < a0coords.ncol(); ++j) {
            out(i,j) = c[j];
        }
        // extract value information
        out(i++, out.ncol() - 1) = iter->second;
    }

    // return lexicographically sorted coord/val pairs
    return out;
}

// [[Rcpp::export]]
NumericMatrix TestFFRWLightLog(NumericMatrix a0coords,
                               NumericVector log_a0values,
                               std::vector<unsigned int> dims, int steps) {
    // Parameters:
    //   dims - number of locations along each dimension

    // initial probability container
    LogArrayMap log_a0;

    // fill initial probability container
    for(int i=0; i < a0coords.nrow(); i++) {
        // munge R coordinate into array's storage format
        VectorI c(a0coords.ncol());
        for(int j=0; j< a0coords.ncol(); ++j)
            c[j] = a0coords(i,j);
        // insert coord/value pair into sparse array
        log_a0.set(c, log_a0values(i));
    }

    // diffuse initial probability
    RN rn(dims);
    LogArrayMap log_af = ffrw_light<IndexType, RN, 
                                    LogArrayMap>(dims, log_a0, steps, rn);

    // extract final, diffused probability
    NumericMatrix out = NumericMatrix(log_af.data.size(), a0coords.ncol() + 1);
    int i=0;
    for(auto iter = log_af.data.begin(); iter != log_af.data.end(); ++iter) {
        // extract coordinate info from array's storage format
        VectorI c = iter->first;
        for(int j=0; j < a0coords.ncol(); ++j) {
            out(i,j) = c[j];
        }
        // extract value information
        out(i++, out.ncol() - 1) = iter->second;
    }

    // return lexicographically sorted coord/val pairs
    return out;
}

// [[Rcpp::export]]
NumericMatrix FFRWLightLogConstrained(
    NumericMatrix a0coords, NumericVector log_a0values,
    std::vector<unsigned int> dims, int steps,
    std::vector<double> surface_heights,
    std::vector<double> domain_heights) {

    // Parameters:
    //   dims - number of locations along each dimension

    // initial probability container
    LogArrayMap log_a0;

    // fill initial probability container
    for(int i=0; i < a0coords.nrow(); i++) {
        // munge R coordinate into array's storage format
        VectorI c(a0coords.ncol());
        for(int j=0; j< a0coords.ncol(); ++j)
            c[j] = a0coords(i,j);
        // insert coord/value pair into sparse array
        log_a0.set(c, log_a0values(i));
    }

    // diffuse initial probability
    ZRN zrn(dims, surface_heights.data(), domain_heights.data());
    LogArrayMap log_af = ffrw_light<IndexType, ZRN, LogArrayMap>(dims, log_a0, 
                                                                 steps, zrn);

    // extract final, diffused probability
    NumericMatrix out = NumericMatrix(log_af.data.size(), a0coords.ncol() + 1);
    int i=0;
    for(auto iter = log_af.data.begin(); iter != log_af.data.end(); ++iter) {
        // extract coordinate info from array's storage format
        VectorI c = iter->first;
        for(int j=0; j < a0coords.ncol(); ++j) {
            out(i,j) = c[j];
        }
        // extract value information
        out(i++, out.ncol() - 1) = iter->second;
    }

    // return lexicographically sorted coord/val pairs
    return out;
}

// [[Rcpp::export]]
std::vector<NumericMatrix> FFRWLogConstrained(
    NumericMatrix a0coords, NumericVector log_a0values,
    std::vector<unsigned int> dims, int steps,
    std::vector<double> surface_heights,
    std::vector<double> domain_heights) {
    // Parameters:
    //   dims - number of locations along each dimension

    // initial probability container
     LogArrayMap log_a0;

    // fill initial probability container
    for(IndexType i=0; i < a0coords.nrow(); i++) {
        // munge R coordinate into array's storage format
        VectorI c(a0coords.ncol());
        for(IndexType j=0; j< a0coords.ncol(); ++j)
            c[j] = a0coords(i,j);
        // insert coord/value pair into sparse array
        log_a0.set(c, log_a0values(i));
    }

    // diffuse initial probability
    ZRN zrn(dims, surface_heights.data(), domain_heights.data());
    std::vector<LogArrayMap> ffprobs = ffrw<IndexType, ZRN, LogArrayMap>(
        dims, log_a0, steps, zrn);

    // extract all diffused probability vectors
    std::vector<NumericMatrix> res;
    res.reserve(steps + 1);
    for(auto p = ffprobs.begin(); p != ffprobs.end(); ++p) {
        NumericMatrix out = NumericMatrix(p->data.size(), a0coords.ncol() + 1);
        IndexType i=0;
        for(auto iter = p->data.begin(); iter != p->data.end(); ++iter) {
            // extract coordinate info from array's storage format
            VectorI c = iter->first;
            for(IndexType j=0; j <a0coords.ncol(); ++j) {
                out(i,j) = c[j];
            }
            // extract value information
            out(i++, out.ncol() - 1) = iter->second;
        }
        res.emplace_back(out);
    }

    // return lexicographically sorted coord/val pairs
    return res;
}

// [[Rcpp::export]]
std::vector<NumericMatrix> FFRWLogConstrainedDst(
    NumericMatrix a0coords, NumericMatrix dstcoords, NumericVector log_a0values,
    std::vector<unsigned int> dims, unsigned int steps, unsigned int max_steps,
    std::vector<double> surface_heights, std::vector<double> domain_heights
) {
   
    // initial probability container
     LogArrayMap log_a0;

    // fill initial probability container
    for(IndexType i=0; i < a0coords.nrow(); i++) {
        // munge R coordinate into array's storage format
        VectorI c(a0coords.ncol());
        for(IndexType j=0; j< a0coords.ncol(); ++j)
            c[j] = a0coords(i,j);
        // insert coord/value pair into sparse array
        log_a0.set(c, log_a0values(i));
    }

    // munge destination coordinates into array's storage format
    VectorI dst_index(dstcoords.ncol());
    for(IndexType j=0; j < dstcoords.ncol(); ++j) {
        dst_index[j] = dstcoords(0,j);
    }

    // diffuse initial probability
    ZRN zrn(dims, surface_heights.data(), domain_heights.data());
    std::vector<LogArrayMap> ffprobs = ffrw_dst<IndexType, ZRN, LogArrayMap>(
        dims, log_a0, steps, zrn, dst_index, max_steps);

    // extract all diffused probability vectors
    std::vector<NumericMatrix> res;
    res.reserve(ffprobs.size());
    for(auto p = ffprobs.begin(); p != ffprobs.end(); ++p) {
        NumericMatrix out = NumericMatrix(p->data.size(), a0coords.ncol() + 1);
        IndexType i=0;
        for(auto iter = p->data.begin(); iter != p->data.end(); ++iter) {
            // extract coordinate info from array's storage format
            VectorI c = iter->first;
            for(IndexType j=0; j <a0coords.ncol(); ++j) {
                out(i,j) = c[j];
            }
            // extract value information
            out(i++, out.ncol() - 1) = iter->second;
        }
        res.emplace_back(out);
    }

    // return lexicographically sorted coord/val pairs
    return res; 
}

// [[Rcpp::export]]
List FFRWLogConstrainedDstReachable(
    NumericMatrix a0coords, NumericMatrix dstcoords, NumericVector log_a0values,
    std::vector<unsigned int> dims, unsigned int steps, unsigned int max_steps,
    std::vector<double> surface_heights, std::vector<double> domain_heights
) {

    // initial probability container
    LogArrayMap log_a0;

    // fill initial probability container
    for(IndexType i=0; i < a0coords.nrow(); i++) {
        // munge R coordinate into array's storage format
        VectorI c(a0coords.ncol());
        for(IndexType j=0; j< a0coords.ncol(); ++j)
            c[j] = a0coords(i,j);
        // insert coord/value pair into sparse array
        log_a0.set(c, log_a0values(i));
    }

    // munge destination coordinates into array's storage format
    VectorI dst_index(dstcoords.ncol());
    for(IndexType j=0; j < dstcoords.ncol(); ++j) {
        dst_index[j] = dstcoords(0,j);
    }

    // diffuse initial probability
    ZRN zrn(dims, surface_heights.data(), domain_heights.data());
    ReachableProbs<IndexType, LogArrayMap> ffprobs =
            ffrw_dst_reachable<IndexType, ZRN, LogArrayMap>(
                    dims, log_a0, steps, zrn, dst_index, max_steps
            );

    // extract all diffused probability vectors
    std::vector<NumericMatrix> res;
    res.reserve(ffprobs.probs.size());
    for(auto p = ffprobs.probs.begin(); p != ffprobs.probs.end(); ++p) {
        NumericMatrix out = NumericMatrix(p->data.size(), a0coords.ncol() + 1);
        IndexType i=0;
        for(auto iter = p->data.begin(); iter != p->data.end(); ++iter) {
            // extract coordinate info from array's storage format
            VectorI c = iter->first;
            for(IndexType j=0; j <a0coords.ncol(); ++j) {
                out(i,j) = c[j];
            }
            // extract value information
            out(i++, out.ncol() - 1) = iter->second;
        }
        res.emplace_back(out);
    }

    return List::create(res, ffprobs.reachable);
}