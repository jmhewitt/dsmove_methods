//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>

#ifndef DSMOVETOOLS_GPS_LIK_H
#define DSMOVETOOLS_GPS_LIK_H

using namespace Rcpp;

class GpsLikBase {

    private:

        // quick switch to indicate that observations are NA
        bool na_obs;

        // bivariate normal parameters and constants
        double mu_lon, mu_lat, sd_lon, sd_lat, rho, rhosq_c, lcst;

        // truncation parameters
        double qform_thresh, ltrunc;

        // distance between two lon/lat coords, in meters
        double distance_m(double, double, double, double);

    public:

        GpsLikBase(double alpha) :
            ltrunc(std::log(1-alpha)),
            qform_thresh(R::qchisq(1-alpha, 2, 1, 0)) { };

        // log-likelihood for coordinates relative to observation uncertainty
        double ll(double lon, double lat);

        // parameterize dist'n. using lon/lat/error information
        void parameterizeDistribution(
            double lon_ctr, double lat_ctr, double semi_major,
            double semi_minor, double orientation
        );

};

class GpsLik : public GpsLikBase {

    private:

        // GPS lon/lat observations with error ellipse information
        std::vector<double> *lon_ctrs, *lat_ctrs, *semi_majors, *semi_minors,
            *orientations;

    public:

        GpsLik(double alpha, std::vector<double> &lons,
               std::vector<double> &lats, std::vector<double> &majors,
               std::vector<double> &minors, std::vector<double> &orients) :
                lon_ctrs(&lons), lat_ctrs(&lats), semi_majors(&majors),
                semi_minors(&minors), orientations(&orients),
                GpsLikBase(alpha) { };

        // parameterize log-likelihood using an observation, specified via index
        void setLikToObs(unsigned int);

};

class GpsLikGridded : public GpsLik {

    private:

        // grid definitions
        std::vector<double> *lon_gridvals, *lat_gridvals;

    public:

        GpsLikGridded(double alpha, std::vector<double> &lons,
          std::vector<double> &lats, std::vector<double> &majors,
          std::vector<double> &minors, std::vector<double> &orients,
          std::vector<double> &lon_grid, std::vector<double> &lat_grid) :
                GpsLik(alpha, lons, lats, majors, minors, orients),
                lon_gridvals(&lon_grid), lat_gridvals(&lat_grid) { };

        double ll(unsigned int lon_ind, unsigned int lat_ind);

};

#endif //DSMOVETOOLS_GPS_LIK_H
