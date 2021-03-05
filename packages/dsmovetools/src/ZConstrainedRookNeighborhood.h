//
// Created by Joshua Hewitt on 12/14/20.
//

#include "RookNeighborhood.h"

#ifndef DSMOVETOOLS_ZCONSTRAINEDROOK_NEIGHBORHOOD_H
#define DSMOVETOOLS_ZCONSTRAINEDROOK_NEIGHBORHOOD_H

/**
 * Following C++ indexing conventions, coordinates for dimension i span
 * 0, ..., dimvec[i]-1.
 *
 * @tparam size_type
 * @tparam Index
 */
template<typename size_type, typename Index>
class ZConstrainedRookNeighborhood : public RookNeighborhood<size_type, Index> {

    private:

        // column-major format specification of spatially-varying range for z;
        // defines minimum z value for each (x,y) pair, max. value is Inf.
        const double * zconstraint;
        // specifies z value for each z index.  i.e., z(0) = a, z(1) = b, etc...
        const double * zdef;
        // index for the final dimension
        const size_type last_dim;

    public:

        ZConstrainedRookNeighborhood(const Index &dims,
                                     const double *zfield,
                                     const double *zvals) :
                                     RookNeighborhood<size_type, Index>(dims),
                                     zconstraint(zfield),
                                     zdef(zvals),
                                     last_dim(dims.size()-1) {};

        bool inDomain(const Index&);

};

template<typename size_type, typename Index>
bool ZConstrainedRookNeighborhood<size_type, Index>::inDomain(const Index& coord) {
    // column-major lookup index associated with coordinate's (x,y) pair
    size_type ind = coord[0] + coord[1] * (*(this->dimvec))[0];
    // verify the coordinate's z-value is larger than the field's minimum
    if(zdef[coord[last_dim]] < zconstraint[ind]) {
        return false;
    }
    // defer to default checks
    return RookNeighborhood<size_type, Index>::inDomain(coord);
}

#endif //DSMOVETOOLS_ZCONSTRAINED3DROOK_NEIGHBORHOOD_H
