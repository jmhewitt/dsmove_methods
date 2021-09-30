//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>
#include "CTDS2DDomain.h"

#ifndef EXACT_LOCATION_LIK_H
#define EXACT_LOCATION_LIK_H

using namespace Rcpp;

class ExactLocationLik {

    private:

        // quick switch to indicate that observations are NA
        bool na_obs;

        // current location information
        double x_coord_cur, y_coord_cur;

        // location observations
        std::vector<double> *obs_x_coord, *obs_y_coord;

    public:

        ExactLocationLik(std::vector<double> &obs_x,
                         std::vector<double> &obs_y) :
            obs_x_coord(&obs_x), obs_y_coord(&obs_y) { };

        // log-likelihood for coordinates relative to observation uncertainty
        double ll(double x_coord, double y_coord) {

            // flat likelihood when observations are NA (i.e., are uninformative)
            if(na_obs) {
                return 0;
            }

            if(x_coord != x_coord_cur) {
                return -std::numeric_limits<double>::infinity();
            }

            if(y_coord != y_coord_cur) {
                return -std::numeric_limits<double>::infinity();
            }

            return 0;
        };

        double ll(const CTDS2DState& state) {
            return ll(state.lon_to, state.lat_to);
        }

        void setLikToObs(unsigned int ind) {

            if(std::isnan((*obs_x_coord)[ind])) {
                na_obs = true;
                return;
            } else {
                na_obs = false;
            }

            x_coord_cur = (*obs_x_coord)[ind];
            y_coord_cur = (*obs_y_coord)[ind];

        };

};

#endif //EXACT_LOCATION_LIK_H
