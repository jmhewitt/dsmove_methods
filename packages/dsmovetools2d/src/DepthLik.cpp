#include "DepthLik.h"


/**
 * @param d Depth observation
 */
void DepthLikBase::setLik(double d) {
    if(std::isnan(d)) {
        na_obs = true;
    } else {
        na_obs = false;
        depth = d;
    }
}

double DepthLikBase::ll(double loc_depth) {

    // uniform likelihood if no depth information available
    if(na_obs) {
      return 0;
    }

    // uniform likelihood if current depth if above the surface at location
    return depth >= loc_depth ? 0 :
        -std::numeric_limits<double>::infinity();
}

void DepthLik::setLikToObs(unsigned int ind) {
    setLik((*depths)[ind]);
}