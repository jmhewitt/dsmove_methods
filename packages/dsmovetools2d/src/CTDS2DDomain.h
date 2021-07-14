//
// Created by Joshua Hewitt on 7/14/21.
//

#include <Rcpp.h>

#ifndef DSMOVETOOLS2D_CTDS2DDOMAIN_H
#define DSMOVETOOLS2D_CTDS2DDOMAIN_H

struct CTDS2DState {
    double log_prob, surface_height, lon_to, lat_to;
    unsigned int prob_age, direction_of_movement, nnbrs;
    int lon_from_ind, lon_to_ind, lat_from_ind, lat_to_ind;
    bool well_defined;
    CTDS2DState *nbr_to[4], *nbr_from[4];
};

class CTDS2DDomain {

private:

    std::vector<CTDS2DState> states;
    unsigned int prob_age = 1;
    unsigned int nlons, nlats;

public:

    CTDS2DDomain(std::vector<double> &lons, std::vector<double> &lats,
                 std::vector<double> &surface_heights);

    // set log prob for a state in the domain
    void set(int lon_from_ind, int lat_from_ind, int lon_to_ind, int lat_to_ind,
             double log_prob);
    void set(CTDS2DState *state, double log_prob);

    // access a state in the domain
    CTDS2DState* statePtr(int lon_from_ind, int lat_from_ind, int lon_to_ind,
                          int lat_to_ind);

    // flatten non-zero entries to matrix format
    Rcpp::NumericMatrix toNumericMatrix();

};


#endif //DSMOVETOOLS2D_CTDS2DDOMAIN_H
