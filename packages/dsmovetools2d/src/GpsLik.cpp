#include "GpsLik.h"

void GpsLikBase::parameterizeDistribution(
            double lon_ctr, double lat_ctr, double semi_major, double semi_minor,
            double orientation) {

    if(std::isnan(lon_ctr)) {
        na_obs = true;
        return;
    } else {
        na_obs = false;
    }

    // Implements "Materials and Methods" from McClintock et al. (2015)

    double Mtsq_half = std::pow(semi_major, 2) / 2;
    double mtsq_half = std::pow(semi_minor, 2) / 2;

    double c = orientation / 180 * M_PI;
    double cos_c = std::cos(c);
    double cos2_c = std::pow(cos_c, 2);
    double sin_c = std::sin(c);
    double sin2_c = std::pow(sin_c, 2);

    obs_lon = lon_ctr;
    obs_lat = lat_ctr;

    sd_lon = std::sqrt(Mtsq_half * sin2_c + mtsq_half * cos2_c);
    sd_lat = std::sqrt(Mtsq_half * cos2_c + mtsq_half * sin2_c);
    rho = (Mtsq_half - mtsq_half) * cos_c * sin_c / sd_lat / sd_lon;

    rhosq_c = 1 - std::pow(rho, 2);

    // log-normalizing constant, including rescaling for truncation
    lcst = - std::log(2 * M_PI * sd_lon * sd_lat) - 0.5 * std::log(rhosq_c) -
            ltrunc;
}

double GpsLikBase::ll(double lon, double lat) {

    // flat likelihood when observations are NA (i.e., are uninformative)
    if(na_obs) {
        return 0;
    }

    // convert distances wrt. lon and lat to m; scale wrt. uncertainty
    double zx = distance_m(obs_lon, lat, lon, lat) / sd_lon;
    double zy = distance_m(lon, obs_lat, lon, lat) / sd_lat;

    // sign the distances
    if(obs_lon < lon) {
        zx *= -1;
    }
    if(obs_lat < lat) {
        zy *= -1;
    }

    // quadratic form
    double q = std::pow(zx, 2) - 2 * rho * zx * zy + std::pow(zy, 2);

    return - q / 2 / rhosq_c + lcst;
//    // return -Inf if lon/lat are outside the ellipse truncation boundary
//    if(q < qform_thresh) {
//        return - q / 2 / rhosq_c + lcst;
//    } else {
//        return -std::numeric_limits<double>::infinity();
//    }
}

double GpsLikGridded::ll(unsigned int lon_ind, unsigned int lat_ind) {
    return GpsLik::ll((*lon_gridvals)[lon_ind],
                      (*lat_gridvals)[lat_ind]);
}

double GpsLikBase::distance_m(double lon1, double lat1, double lon2, double lat2) {

    // convert units from degrees to radians
    lat1 = lat1 * M_PI / 180;
    lon1 = lon1 * M_PI / 180;
    lat2 = lat2 * M_PI / 180;
    lon2 = lon2 * M_PI / 180;

    // Haversine Formula
    double dlon = lon2 - lon1;
    double dlat = lat2 - lat1;

    double ans = pow(sin(dlat / 2), 2) +
                      cos(lat1) * cos(lat2) *
                      pow(sin(dlon / 2), 2);

    ans = 2 * asin(sqrt(ans));

    // radius of Earth in meters
    double R = 6378388;

    // Calculate the result
    return ans * R;
}

void GpsLik::setLikToObs(unsigned int ind) {
    parameterizeDistribution((*lon_ctrs)[ind], (*lat_ctrs)[ind],
                             (*semi_majors)[ind], (*semi_minors)[ind],
                             (*orientations)[ind]);
}

// [[Rcpp::export]]
double GpsLikEval(
        std::vector<double> obs_lons, std::vector<double> obs_lats,
        std::vector<double> semi_majors, std::vector<double> semi_minors,
        std::vector<double> orientations, double alpha, double test_lon,
        double test_lat, int ind) {

    GpsLik g(alpha, obs_lons, obs_lats, semi_majors, semi_minors, orientations);

    g.setLikToObs(ind);

    return g.ll(test_lon, test_lat);
}