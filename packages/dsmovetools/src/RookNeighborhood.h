//
// Created by Joshua Hewitt on 12/14/20.
//

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

        const size_type ndim;               // number of dimensions for grid
        const Index *dimvec;                // num. coords along each dimension
        Index center;                       // neighborhood center
        std::vector<bool> lwr_nbr, upr_nbr; // indicators for whether nbr exists
        size_type nnbrs;                    // number of neighbors
        double log_nnbrs;                   // log of number of neighbors
        size_type nbrs_visited;             // number of neighbors iterated over
        size_type dim_cur;                  // current dimension iterator is on
        bool onLwr;                         // true when iterator is on lower

    public:

        RookNeighborhood(const Index &dims) : ndim(dims.size()), dimvec(&dims) {
            lwr_nbr.reserve(ndim);
            upr_nbr.reserve(ndim);
        }
        void setCenter(const Index&);
        size_type neighborhoodSize() { return nnbrs; }
        double logNeighborhoodSize() { return log_nnbrs; }
        Index nextNeighbor();

};

template<typename size_type, typename Index>
void RookNeighborhood<size_type, Index>::setCenter(const Index& coord) {
    // keep a copy of the neighborhood center
    center = coord;
    // pre-explore the neighborhood
    nnbrs = 0;
    for(size_type i = 0; i < ndim; ++i) {
        if(center[i] == 0) {
            lwr_nbr[i] = false;
            upr_nbr[i] = true;
            nnbrs++;
        } else if(center[i] == (*dimvec)[i] - 1) {
            lwr_nbr[i] = true;
            upr_nbr[i] = false;
            nnbrs++;
        } else if(center[i] > 0 && center[i] < (*dimvec)[i] - 1) {
            lwr_nbr[i] = true;
            upr_nbr[i] = true;
            nnbrs += 2;
        } else {
            nnbrs = 0;
            break;
        }
    }
    // update constants
    log_nnbrs = log(nnbrs);
    // initialize neighborhood counters
    onLwr = true;
    dim_cur = 0;
    nbrs_visited = 0;
}

template<typename size_type, typename Index>
Index RookNeighborhood<size_type, Index>::nextNeighbor() {
    if(nnbrs == 0) {
        // avoid infinite recursive if there is no neighborhood
        return center;
    } else if(dim_cur < ndim) {
        if(onLwr) {
            if(lwr_nbr[dim_cur]) {
                Index nbr = center;
                nbr[dim_cur++] -= 1;
                return nbr;
            } else {
                dim_cur++;
                return nextNeighbor();
            }
        } else {
            if(upr_nbr[dim_cur]) {
                Index nbr = center;
                nbr[dim_cur++] += 1;
                return nbr;
            } else {
                dim_cur++;
                return nextNeighbor();
            }
        }
    } else {
        // finished this half-neighborhood; move to other half-neighborhood
        onLwr = !onLwr;
        dim_cur = 0;
        return nextNeighbor();
    }
}

#endif //DSMOVETOOLS_ROOK_NEIGHBORHOOD_H
