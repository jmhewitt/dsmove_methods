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

    double c = (90 - orientation) / 180 * M_PI;
    double cos_c = std::cos(c);
    double cos2_c = std::pow(cos_c, 2);
    double sin_c = std::sin(c);
    double sin2_c = std::pow(sin_c, 2);

    mu_lon = lon_ctr;
    mu_lat = lat_ctr;

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
    double zx = distance_m(lon, mu_lat, mu_lon, mu_lat) / sd_lon;
    double zy = distance_m(mu_lon, lat, mu_lon, mu_lat) / sd_lat;

    // quadratic form
    double q = std::pow(zx, 2) - 2 * rho * zx * zy + std::pow(zy, 2);

    // return -Inf if lon/lat are outside the ellipse truncation boundary
    if(q < qform_thresh) {
        return - q / 2 / rhosq_c + lcst;
    } else {
        return -std::numeric_limits<double>::infinity();
    }
}

double GpsLikGridded::ll(unsigned int lon_ind, unsigned int lat_ind) {
    return GpsLik::ll((*lon_gridvals)[lon_ind],
                      (*lat_gridvals)[lat_ind]);
}

double GpsLikBase::distance_m(double lon1, double lat1, double lon2, double lat2) {

    double coslat1 = cos(lat1 * M_PI / 180);
    double coslon1 = cos(lon1 * M_PI / 180);
    double coslat2 = cos(lat2 * M_PI / 180);
    double coslon2 = cos(lon2 * M_PI / 180);
    double sinlon1 = sin(lon1 * M_PI / 180);
    double sinlon2 = sin(lon2 * M_PI / 180);
    double sinlat1 = sin(lat1 * M_PI / 180);
    double sinlat2 = sin(lat2 * M_PI / 180);
        
    double pp = coslat1 * coslon1 * coslat2 * coslon2 +
                coslat1 * sinlon1 * coslat2 * sinlon2 + 
                sinlat1 * sinlat2;
    
    // radius of Earth in meters
    double R = 6378388;
    
    if(abs(pp) > 1) {
        return R * acos(sign(pp));
    } else {
        return R * acos(pp);
    }
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