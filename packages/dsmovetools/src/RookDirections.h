//
// Created by Joshua Hewitt on 12/14/20.
//

#include "RookHeading.h"

#ifndef DSMOVETOOLS_ROOK_DIRECTIONS_H
#define DSMOVETOOLS_ROOK_DIRECTIONS_H

/**
 * Following C++ indexing conventions, coordinates for dimension i span
 * 0, ..., dimvec[i]-1.  Given two locations that are Rook-adjacent, determine 
 * the orientation of the unit vector that points from the "tail" location to 
 * the "head" location.  Additionally, implement a dot product for two 
 * direction vectors.
 *
 * @tparam size_type
 * @tparam Index
 */
template<typename size_type, typename Index>
class RookDirections {

    private: 

        const size_type ndim;

    public:

        RookHeading<size_type> orientation;

        RookDirections(const size_type dims) : ndim(dims) { };

        void setOrientation(const Index&, const Index&);
        double dotOrientation(const RookHeading<size_type>&);

};

template<typename size_type, typename Index>
double RookDirections<size_type, Index>::dotOrientation(
    const RookHeading<size_type> &h
) {
    if(orientation.dim == h.dim) {
        // XOR operator "^" returns TRUE if the direction of the vectors align
        return orientation.lwr ^ h.lwr ? 1 : -1;
    } else {
        // Rook-adjacency direction vectors are orthogonal if the vectors are 
        // not oriented in the same dimension
        return 0;
    }
}

template<typename size_type, typename Index>
void RookDirections<size_type, Index>::setOrientation(
    const Index &head, const Index &tail
) {
    for(size_type i=0; i < ndim; ++i) {
        if(head[i] != tail[i]) {
            orientation.dim = i;
            orientation.lwr = head[i] > tail[i];
            break;
        }
    }
}

#endif //DSMOVETOOLS_ROOK_DIRECTIONS_H
