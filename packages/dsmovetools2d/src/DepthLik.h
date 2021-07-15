//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>

#ifndef DSMOVETOOLS2D_DEPTH_LIK_H
#define DSMOVETOOLS2D_DEPTH_LIK_H

using namespace Rcpp;

class DepthLikBase {

    private:

        // quick switch to indicate that observations are NA
        bool na_obs;

        // depth to which likelihood is set
        double depth;

    public:

        DepthLikBase() { }

        // log-likelihood for a depth value relative to observation
        double ll(double loc_depth);

        void setLik(double obs_depth);

};

class DepthLik : public DepthLikBase {

    private:

        // depth bin observations
        std::vector<double> *depths;

    public:

        DepthLik(std::vector<double> &obs_depths) : depths(&obs_depths) { };

        // parameterize log-likelihood using an observation, specified via index
        void setLikToObs(unsigned int);

};

#endif //DSMOVETOOLS2D_DEPTH_LIK_H
