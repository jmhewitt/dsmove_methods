#include "DepthLik.h"


void DepthLik::setLikToObs(unsigned int ind) {
    if(std::isnan((*depths)[ind])) {
        na_obs = true;
    } else {
        na_obs = false;
        depth = (*depths)[ind];
    }
}

double DepthLik::ll(std::vector<unsigned int> coord) {

    // uniform likelihood if no depth information available
    if(na_obs) {
        return -std::numeric_limits<double>::infinity();
    }

    // set depth in current location
    coord[2] = depth;

    return domain->inDomain(coord) ? 0 :
        -std::numeric_limits<double>::infinity();
}

// [[Rcpp::export]]
double DepthLikEval(std::vector<unsigned int> dims,
                    std::vector<unsigned int> coords,
                    std::vector<unsigned int> obs_depths,
                    std::vector<double> zfield,
                    std::vector<double> zvals,
                    int ind) {

    ZConstrainedRookNeighborhood<unsigned int, std::vector<unsigned int>>
        rn(dims, zfield.data(), zvals.data());

    DepthLik d(obs_depths, rn);

    d.setLikToObs(ind);

    return d.ll(coords);
}