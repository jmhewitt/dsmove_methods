//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>
#include "ZConstrainedRookNeighborhood.h"

#ifndef DSMOVETOOLS_DEPTH_LIK_H
#define DSMOVETOOLS_DEPTH_LIK_H

using namespace Rcpp;

class DepthLikBase {

    private:

        // domain structure for bathymetry data
        ZConstrainedRookNeighborhood<unsigned int,
                std::vector<unsigned int>> *domain;

        // quick switch to indicate that observations are NA
        bool na_obs;

        // depth bin to which likelihood is set
        unsigned int depth;

    public:

        DepthLikBase(ZConstrainedRookNeighborhood<unsigned int,
                     std::vector<unsigned int>> &nbhd) : domain(&nbhd) { };

        // log-likelihood for coordinates relative to observation uncertainty
        double ll(std::vector<unsigned int>);

        void setLik(unsigned int);

};

class DepthLik : public DepthLikBase {

    private:

        // depth bin observations
        std::vector<unsigned int> *depths;

    public:

        DepthLik(
            std::vector<unsigned int> &obs_depths,
            ZConstrainedRookNeighborhood<unsigned int,
                                         std::vector<unsigned int>> &nbhd
        ) : DepthLikBase(nbhd), depths(&obs_depths) { };

        // parameterize log-likelihood using an observation, specified via index
        void setLikToObs(unsigned int);

};

#endif //DSMOVETOOLS_DEPTH_LIK_H
