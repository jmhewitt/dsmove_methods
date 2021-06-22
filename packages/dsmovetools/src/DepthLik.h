//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>
#include "ZConstrainedRookNeighborhood.h"

#ifndef DSMOVETOOLS_DEPTH_LIK_H
#define DSMOVETOOLS_DEPTH_LIK_H

using namespace Rcpp;

class DepthLik {

    private:

        // quick switch to indicate that observations are NA
        bool na_obs;

        // depth bin to which likelihood is set
        unsigned int depth;

        // domain structure for bathymetry data
        ZConstrainedRookNeighborhood<unsigned int,
            std::vector<unsigned int>> *domain;

        // depth bin observations
        std::vector<unsigned int> *depths;

    public:

        DepthLik(
            std::vector<unsigned int> &obs_depths,
            ZConstrainedRookNeighborhood<unsigned int,
                                         std::vector<unsigned int>> &nbhd
        ) : depths(&obs_depths), domain(&nbhd) { };

        // parameterize log-likelihood using an observation, specified via index
        void setLikToObs(unsigned int);

        // log-likelihood for coordinates relative to observation uncertainty
        double ll(std::vector<unsigned int>);

};

#endif //DSMOVETOOLS_DEPTH_LIK_H
