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


// diffuse mass
template <typename I, typename V, typename S, typename size_type, 
          typename Neighborhood>
void diffuseMass(SparseNdimArrayBase<I,V,S> *src,
                 SparseNdimArrayBase<I,V,S> *dst,
                 Neighborhood *nbhd) {
    // forward-filter all mass from the src vector to the dst vector
    auto src_mass_end = src->data.end();
    for(auto src_mass_entry = src->data.begin(); 
        src_mass_entry != src_mass_end;
        ++src_mass_entry) {
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
 * Forward filtering a random walk along a grid, while storing all intermediate
 * distributions.
 *
 * @param dims specify the number of locations along each dimension in grid
 * @param a0 initial probability mass vector
 * @param steps number of forward-diffusion steps to take
 * @return (sparse) diffused mass vectors
 */
 template<typename size_type, typename Neighborhood>
std::vector<LogArrayMap> ffrw_log(const VectorI &dims, const LogArrayMap &a0, 
                                  const unsigned int steps, 
                                  Neighborhood &nbhd) {

    // initialize forward filtering vectors and initial mass
    std::vector<LogArrayMap> ffprob(steps + 1);
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

/**
 * Forward filtering a random walk along a grid, without storing all
 * intermediate distributions.  Probabilities are stored and manipulated on the
 * log scale vs. the (natural) linear scale.
 *
 * @param dims specify the number of locations along each dimension in grid
 * @param a0 initial probability mass vector
 * @param steps number of forward-diffusion steps to take
 * @return (sparse) diffused mass vectors
 */
template<typename size_type, typename Neighborhood>
LogArrayMap ffrw_light_log(const VectorI &dims, const LogArrayMap &a0,
                           const unsigned int steps, Neighborhood &nbhd) {

    // initialize forward filtering vectors and initial mass
    LogArrayMap cur, prev;
    prev = a0;

   // diffuse mass
   for(unsigned int step_cur = 0; step_cur < steps; ++step_cur) {
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
 * Forward filtering a random walk along a grid, without storing all
 * intermediate distributions.
 *
 * @param dims specify the number of locations along each dimension in grid
 * @param a0 initial probability mass vector
 * @param steps number of forward-diffusion steps to take
 * @return (sparse) diffused mass vectors
 */
template<typename size_type, typename Neighborhood>
ArrayMap ffrw_light(const VectorI &dims, const ArrayMap &a0,
                    const unsigned int steps) {

    // initialize forward filtering vectors and initial mass
    ArrayMap cur, prev;
    prev = a0;

    // initialize object to iterate over neighborhoods
    Neighborhood nbhd(dims);

    // diffuse mass
    for(unsigned int step_cur = 0; step_cur < steps; ++step_cur) {
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
 template<typename size_type, typename Neighborhood>
std::vector<ArrayMap> ffrw(const VectorI &dims, const ArrayMap &a0,
                           const unsigned int steps) {

    // initialize forward filtering vectors and initial mass
    std::vector<ArrayMap> ffprob(steps + 1);
    ffprob[0] = a0;

    // initialize object to iterate over neighborhoods
    Neighborhood nbhd(dims);

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
    std::vector<ArrayMap> ffprobs = ffrw<IndexType, RN>(dims, a0, steps);

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
    ArrayMap af = ffrw_light<IndexType, RN>(dims, a0, steps);

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
    LogArrayMap log_af = ffrw_light_log<IndexType, RN>(dims, log_a0, steps, rn);

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
    LogArrayMap log_af = ffrw_light_log<IndexType, ZRN>(dims, log_a0, steps,
                                                        zrn);

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
    std::vector<LogArrayMap> ffprobs = ffrw_log<IndexType, ZRN>(dims, log_a0, 
                                                                steps, zrn);

    // extract all diffused probability vectors
    std::vector<NumericMatrix> res;
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
        res.push_back(out);
    }

    // return lexicographically sorted coord/val pairs
    return res;
}