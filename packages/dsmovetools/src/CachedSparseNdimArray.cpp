//
// Created by Joshua Hewitt on 12/14/20.
//

#include "CachedSparseNdimArray.h"

// [[Rcpp::export]]
double TestDualSparseCoordVec(std::vector<unsigned int> x1,
                                     std::vector<unsigned int> x2,
                                     double v) {

    typedef std::vector<unsigned int> CoordVec;
    typedef std::pair<CoordVec, CoordVec> CoordVecPair;

    DualSparseCoordVec c;

    CoordVecPair cp(x1, x2);

    c.addToActive(cp, v);
    c.addToActive(cp, v);
    c.swapActive();

    return c.inactiveValue(cp);
}