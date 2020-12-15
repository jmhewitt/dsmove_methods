//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>
#include "SparseNdimArray.h"
#include "RookNeighborhood.h"

using namespace Rcpp;

// let's use ints to index coordinates
typedef int IndexType;
// for specifying integer coordinates with grid dimension assigned at runtime
typedef std::vector<IndexType> VectorI;
// sparse spatial storage implemented via trees
typedef SparseNdimArray<VectorI, double, std::map<VectorI, double>> ArrayMap;
// rook neighborhoods on integer grids
typedef RookNeighborhood<IndexType, VectorI> RN;

// iterator for neighborhood structure

/**
 *
 * neighborhood structure?
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
        auto prev_mass_end = step_prev->data.end();
        for(auto prev_mass_entry = step_prev->data.begin();
            prev_mass_entry != prev_mass_end;
            ++prev_mass_entry) {
            // find neighborhood for previous location
            nbhd.setCenter(prev_mass_entry->first);
            size_type nnbrs = nbhd.neighborhoodSize();
            // diffuse mass, following a random walk along neighbors
            double mass = prev_mass_entry->second / (double) nnbrs;
            for(size_type i = 0; i < nnbrs; ++i) {
                step_cur->add(nbhd.nextNeighbor(), mass);
            }
        }
        // increment iterators
        ++step_cur;
        ++step_prev;
    }

    return ffprob;
};


// [[Rcpp::export]]
NumericMatrix TestFFRW(NumericMatrix a0coords, NumericVector a0values,
                       std::vector<int> dims, int steps) {
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