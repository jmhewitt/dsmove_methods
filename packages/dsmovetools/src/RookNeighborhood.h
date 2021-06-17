//
// Created by Joshua Hewitt on 12/14/20.
//

#include "RookHeading.h"

#ifndef DSMOVETOOLS_ROOK_NEIGHBORHOOD_H
#define DSMOVETOOLS_ROOK_NEIGHBORHOOD_H

/**
 * Following C++ indexing conventions, coordinates for dimension i span
 * 0, ..., dimvec[i]-1.
 *
 * @tparam size_type
 * @tparam Index
 */
template<typename size_type, typename Index>
class RookNeighborhood {

    private: 

        Index center;                       // neighborhood center
        std::vector<bool> lwr_nbr, upr_nbr; // indicators for whether nbr exists
        size_type nnbrs;                    // number of neighbors
        size_type dim_iter;                 // current dimension iterator is on
        bool onLwr;                         // true when iterator is on lower
        size_type dim_cur;                  // dimension of current neighbor

        Index translateCoord(const Index&, const size_type, const int);

    protected:

        const size_type ndim;               // number of dimensions for grid
        const Index *dimvec;                // num. coords along each dimension

        virtual bool satisfiesDomainConstraints(const Index&) { return true; };

    public:

        RookNeighborhood(const Index &dims) : ndim(dims.size()), dimvec(&dims) {
            lwr_nbr.reserve(ndim);
            upr_nbr.reserve(ndim);
        }

        void setCenter(const Index&);
        size_type neighborhoodSize() { return nnbrs; }

        Index nextNeighbor();
        void resetNeighborIterator();

        bool inDomain(const Index&);

        RookHeading<size_type> neighborHeading();

};

template<typename size_type, typename Index>
RookHeading<size_type> RookNeighborhood<size_type, Index>::neighborHeading() {
    return RookHeading<size_type>(dim_cur, onLwr);
}


template<typename size_type, typename Index>
bool RookNeighborhood<size_type, Index>::inDomain(const Index& coord) {
    // verify coordinate is within bounds for each dimension
    for(size_type i=0; i < ndim; ++i) {
        if(coord[i] < 0) { return false; }
        if(coord[i] > (*dimvec)[i] - 1) { return false; }
    }
    // run additional checks
    return satisfiesDomainConstraints(coord);
}

template<typename size_type, typename Index>
Index RookNeighborhood<size_type, Index>::translateCoord(const Index& coord,
                                                         const size_type dim,
                                                         const int value) {
    // shift the coordinate by "value" units in dimension "dim"
    Index nbr = coord;
    nbr[dim] += value;
    return nbr;
}

template<typename size_type, typename Index>
void RookNeighborhood<size_type, Index>::setCenter(const Index& coord) {
    // keep a copy of the neighborhood center
    center = coord;

    // reset neighborhood size
    nnbrs = 0;

    // pre-explore the neighborhood
    if(inDomain(center)) {
        for(size_type i = 0; i < ndim; ++i) {
            if(inDomain(translateCoord(center, i, -1))) {
                lwr_nbr[i] = true;
                nnbrs++;
            } else {
                lwr_nbr[i] = false;
            }
            if(inDomain(translateCoord(center, i, 1))) {
                upr_nbr[i] = true;
                nnbrs++;
            } else {
                upr_nbr[i] = false;
            }
        }
    }

    // initialize neighborhood counters
    resetNeighborIterator();
}

template<typename size_type, typename Index>
Index RookNeighborhood<size_type, Index>::nextNeighbor() {
    if(nnbrs == 0) {
        // avoid infinite recursion if there is no neighborhood
        return center;
    } else if(dim_iter < ndim) {
        if(onLwr) {
            if(lwr_nbr[dim_iter]) {
                dim_cur = dim_iter;
                return translateCoord(center, dim_iter++, -1);
            } else {
                dim_iter++;
                return nextNeighbor();
            }
        } else {
            if(upr_nbr[dim_iter]) {
                dim_cur = dim_iter;
                return translateCoord(center, dim_iter++, 1);
            } else {
                dim_iter++;
                return nextNeighbor();
            }
        }
    } else {
        // finished this half-neighborhood; move to other half-neighborhood
        onLwr = !onLwr;
        dim_iter = 0;
        return nextNeighbor();
    }
}

template<typename size_type, typename Index>
void RookNeighborhood<size_type, Index>::resetNeighborIterator() {
    onLwr = true;
    dim_cur = 0;
    dim_iter = 0;
}

#endif //DSMOVETOOLS_ROOK_NEIGHBORHOOD_H
