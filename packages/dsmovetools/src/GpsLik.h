//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>

#ifndef DSMOVETOOLS_GPS_LIK_H
#define DSMOVETOOLS_GPS_LIK_H

using namespace Rcpp;

class GpsLik {

    private:

        // quick switch to indicate that observations are NA
        bool na_obs;

        // bivariate normal parameters and constants
        double mu_lon, mu_lat, sd_lon, sd_lat, rho, rhosq_c, lcst;

        // truncation parameters
        double qform_thresh, ltrunc;

        // GPS lon/lat observations with error ellipse information
        std::vector<double> *lon_ctrs, *lat_ctrs, *semi_majors, *semi_minors,
            *orientations;

        // parameterize dist'n. using lon/lat/error information
        void parameterizeDistribution(double, double, double, double, double);

        // distance between two lon/lat coords, in meters
        double distance_m(double, double, double, double);

    public:

        GpsLik(double alpha, std::vector<double> &lons,
               std::vector<double> &lats, std::vector<double> &majors,
               std::vector<double> &minors, std::vector<double> &orients) :
                lon_ctrs(&lons), lat_ctrs(&lats), semi_majors(&majors),
                semi_minors(&minors), orientations(&orients),
                ltrunc(std::log(1-alpha)),
                qform_thresh(R::qchisq(1-alpha, 2, 1, 0)) { };

        // parameterize log-likelihood using an observation, specified via index
        void setLikToObs(unsigned int);

        // log-likelihood for coordinates relative to observation uncertainty
        double ll(double, double);

};

#endif //DSMOVETOOLS_GPS_LIK_H
