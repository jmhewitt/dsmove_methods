// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/dsmovetools2d.h"
#include <Rcpp.h>

using namespace Rcpp;

// TestCTDS2DDomainIO
NumericMatrix TestCTDS2DDomainIO(std::vector<double> lons, std::vector<double> lats, std::vector<double> surface_heights, NumericMatrix init_dsts, NumericMatrix init_srcs, std::vector<double> log_probs);
RcppExport SEXP _dsmovetools2d_TestCTDS2DDomainIO(SEXP lonsSEXP, SEXP latsSEXP, SEXP surface_heightsSEXP, SEXP init_dstsSEXP, SEXP init_srcsSEXP, SEXP log_probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type lons(lonsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type lats(latsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type surface_heights(surface_heightsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type init_dsts(init_dstsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type init_srcs(init_srcsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type log_probs(log_probsSEXP);
    rcpp_result_gen = Rcpp::wrap(TestCTDS2DDomainIO(lons, lats, surface_heights, init_dsts, init_srcs, log_probs));
    return rcpp_result_gen;
END_RCPP
}
// LogTxProbs
NumericMatrix LogTxProbs(std::vector<double> lons, std::vector<double> lats, std::vector<double> surface_heights, int lon_from_ind, int lat_from_ind, int lon_to_ind, int lat_to_ind, double betaAR);
RcppExport SEXP _dsmovetools2d_LogTxProbs(SEXP lonsSEXP, SEXP latsSEXP, SEXP surface_heightsSEXP, SEXP lon_from_indSEXP, SEXP lat_from_indSEXP, SEXP lon_to_indSEXP, SEXP lat_to_indSEXP, SEXP betaARSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type lons(lonsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type lats(latsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type surface_heights(surface_heightsSEXP);
    Rcpp::traits::input_parameter< int >::type lon_from_ind(lon_from_indSEXP);
    Rcpp::traits::input_parameter< int >::type lat_from_ind(lat_from_indSEXP);
    Rcpp::traits::input_parameter< int >::type lon_to_ind(lon_to_indSEXP);
    Rcpp::traits::input_parameter< int >::type lat_to_ind(lat_to_indSEXP);
    Rcpp::traits::input_parameter< double >::type betaAR(betaARSEXP);
    rcpp_result_gen = Rcpp::wrap(LogTxProbs(lons, lats, surface_heights, lon_from_ind, lat_from_ind, lon_to_ind, lat_to_ind, betaAR));
    return rcpp_result_gen;
END_RCPP
}
// LogTxProbsElevation
NumericMatrix LogTxProbsElevation(std::vector<double> lons, std::vector<double> lats, std::vector<double> surface_heights, int lon_from_ind, int lat_from_ind, int lon_to_ind, int lat_to_ind, double betaAR, double min_elevation, double max_elevation);
RcppExport SEXP _dsmovetools2d_LogTxProbsElevation(SEXP lonsSEXP, SEXP latsSEXP, SEXP surface_heightsSEXP, SEXP lon_from_indSEXP, SEXP lat_from_indSEXP, SEXP lon_to_indSEXP, SEXP lat_to_indSEXP, SEXP betaARSEXP, SEXP min_elevationSEXP, SEXP max_elevationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type lons(lonsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type lats(latsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type surface_heights(surface_heightsSEXP);
    Rcpp::traits::input_parameter< int >::type lon_from_ind(lon_from_indSEXP);
    Rcpp::traits::input_parameter< int >::type lat_from_ind(lat_from_indSEXP);
    Rcpp::traits::input_parameter< int >::type lon_to_ind(lon_to_indSEXP);
    Rcpp::traits::input_parameter< int >::type lat_to_ind(lat_to_indSEXP);
    Rcpp::traits::input_parameter< double >::type betaAR(betaARSEXP);
    Rcpp::traits::input_parameter< double >::type min_elevation(min_elevationSEXP);
    Rcpp::traits::input_parameter< double >::type max_elevation(max_elevationSEXP);
    rcpp_result_gen = Rcpp::wrap(LogTxProbsElevation(lons, lats, surface_heights, lon_from_ind, lat_from_ind, lon_to_ind, lat_to_ind, betaAR, min_elevation, max_elevation));
    return rcpp_result_gen;
END_RCPP
}
// GpsLikEval
double GpsLikEval(std::vector<double> obs_lons, std::vector<double> obs_lats, std::vector<double> semi_majors, std::vector<double> semi_minors, std::vector<double> orientations, double alpha, double test_lon, double test_lat, int ind);
RcppExport SEXP _dsmovetools2d_GpsLikEval(SEXP obs_lonsSEXP, SEXP obs_latsSEXP, SEXP semi_majorsSEXP, SEXP semi_minorsSEXP, SEXP orientationsSEXP, SEXP alphaSEXP, SEXP test_lonSEXP, SEXP test_latSEXP, SEXP indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_lons(obs_lonsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_lats(obs_latsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type semi_majors(semi_majorsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type semi_minors(semi_minorsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type orientations(orientationsSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type test_lon(test_lonSEXP);
    Rcpp::traits::input_parameter< double >::type test_lat(test_latSEXP);
    Rcpp::traits::input_parameter< int >::type ind(indSEXP);
    rcpp_result_gen = Rcpp::wrap(GpsLikEval(obs_lons, obs_lats, semi_majors, semi_minors, orientations, alpha, test_lon, test_lat, ind));
    return rcpp_result_gen;
END_RCPP
}
// FF_DTMC
NumericMatrix FF_DTMC(std::vector<double> lons, std::vector<double> lats, std::vector<double> surface_heights, NumericMatrix init_dsts, NumericMatrix init_srcs, std::vector<double> init_log_probs, unsigned int steps, double log_self_tx, double betaAR);
RcppExport SEXP _dsmovetools2d_FF_DTMC(SEXP lonsSEXP, SEXP latsSEXP, SEXP surface_heightsSEXP, SEXP init_dstsSEXP, SEXP init_srcsSEXP, SEXP init_log_probsSEXP, SEXP stepsSEXP, SEXP log_self_txSEXP, SEXP betaARSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type lons(lonsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type lats(latsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type surface_heights(surface_heightsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type init_dsts(init_dstsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type init_srcs(init_srcsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type init_log_probs(init_log_probsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< double >::type log_self_tx(log_self_txSEXP);
    Rcpp::traits::input_parameter< double >::type betaAR(betaARSEXP);
    rcpp_result_gen = Rcpp::wrap(FF_DTMC(lons, lats, surface_heights, init_dsts, init_srcs, init_log_probs, steps, log_self_tx, betaAR));
    return rcpp_result_gen;
END_RCPP
}
// SattagFilteredLL
double SattagFilteredLL(NumericMatrix init_dsts, NumericMatrix init_srcs, std::vector<double> init_log_probs, double gps_trunc_alpha, std::vector<double> obs_lons, std::vector<double> obs_lats, std::vector<double> obs_semi_majors, std::vector<double> obs_semi_minors, std::vector<double> obs_orientations, std::vector<double> obs_depths, std::vector<double> lon_gridvals, std::vector<double> lat_gridvals, std::vector<double> surface_heights, double min_elevation, double max_elevation, double log_self_tx, double betaAR);
RcppExport SEXP _dsmovetools2d_SattagFilteredLL(SEXP init_dstsSEXP, SEXP init_srcsSEXP, SEXP init_log_probsSEXP, SEXP gps_trunc_alphaSEXP, SEXP obs_lonsSEXP, SEXP obs_latsSEXP, SEXP obs_semi_majorsSEXP, SEXP obs_semi_minorsSEXP, SEXP obs_orientationsSEXP, SEXP obs_depthsSEXP, SEXP lon_gridvalsSEXP, SEXP lat_gridvalsSEXP, SEXP surface_heightsSEXP, SEXP min_elevationSEXP, SEXP max_elevationSEXP, SEXP log_self_txSEXP, SEXP betaARSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type init_dsts(init_dstsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type init_srcs(init_srcsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type init_log_probs(init_log_probsSEXP);
    Rcpp::traits::input_parameter< double >::type gps_trunc_alpha(gps_trunc_alphaSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_lons(obs_lonsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_lats(obs_latsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_semi_majors(obs_semi_majorsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_semi_minors(obs_semi_minorsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_orientations(obs_orientationsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_depths(obs_depthsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type lon_gridvals(lon_gridvalsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type lat_gridvals(lat_gridvalsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type surface_heights(surface_heightsSEXP);
    Rcpp::traits::input_parameter< double >::type min_elevation(min_elevationSEXP);
    Rcpp::traits::input_parameter< double >::type max_elevation(max_elevationSEXP);
    Rcpp::traits::input_parameter< double >::type log_self_tx(log_self_txSEXP);
    Rcpp::traits::input_parameter< double >::type betaAR(betaARSEXP);
    rcpp_result_gen = Rcpp::wrap(SattagFilteredLL(init_dsts, init_srcs, init_log_probs, gps_trunc_alpha, obs_lons, obs_lats, obs_semi_majors, obs_semi_minors, obs_orientations, obs_depths, lon_gridvals, lat_gridvals, surface_heights, min_elevation, max_elevation, log_self_tx, betaAR));
    return rcpp_result_gen;
END_RCPP
}
// SattagPredDist
std::vector<NumericMatrix> SattagPredDist(NumericMatrix init_dsts, NumericMatrix init_srcs, std::vector<double> init_log_probs, double gps_trunc_alpha, std::vector<double> obs_lons, std::vector<double> obs_lats, std::vector<double> obs_semi_majors, std::vector<double> obs_semi_minors, std::vector<double> obs_orientations, std::vector<double> obs_depths, std::vector<double> lon_gridvals, std::vector<double> lat_gridvals, std::vector<double> surface_heights, double min_elevation, double max_elevation, double log_self_tx, double betaAR, std::vector<unsigned int> pred_steps);
RcppExport SEXP _dsmovetools2d_SattagPredDist(SEXP init_dstsSEXP, SEXP init_srcsSEXP, SEXP init_log_probsSEXP, SEXP gps_trunc_alphaSEXP, SEXP obs_lonsSEXP, SEXP obs_latsSEXP, SEXP obs_semi_majorsSEXP, SEXP obs_semi_minorsSEXP, SEXP obs_orientationsSEXP, SEXP obs_depthsSEXP, SEXP lon_gridvalsSEXP, SEXP lat_gridvalsSEXP, SEXP surface_heightsSEXP, SEXP min_elevationSEXP, SEXP max_elevationSEXP, SEXP log_self_txSEXP, SEXP betaARSEXP, SEXP pred_stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type init_dsts(init_dstsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type init_srcs(init_srcsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type init_log_probs(init_log_probsSEXP);
    Rcpp::traits::input_parameter< double >::type gps_trunc_alpha(gps_trunc_alphaSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_lons(obs_lonsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_lats(obs_latsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_semi_majors(obs_semi_majorsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_semi_minors(obs_semi_minorsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_orientations(obs_orientationsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_depths(obs_depthsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type lon_gridvals(lon_gridvalsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type lat_gridvals(lat_gridvalsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type surface_heights(surface_heightsSEXP);
    Rcpp::traits::input_parameter< double >::type min_elevation(min_elevationSEXP);
    Rcpp::traits::input_parameter< double >::type max_elevation(max_elevationSEXP);
    Rcpp::traits::input_parameter< double >::type log_self_tx(log_self_txSEXP);
    Rcpp::traits::input_parameter< double >::type betaAR(betaARSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type pred_steps(pred_stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(SattagPredDist(init_dsts, init_srcs, init_log_probs, gps_trunc_alpha, obs_lons, obs_lats, obs_semi_majors, obs_semi_minors, obs_orientations, obs_depths, lon_gridvals, lat_gridvals, surface_heights, min_elevation, max_elevation, log_self_tx, betaAR, pred_steps));
    return rcpp_result_gen;
END_RCPP
}
// BackInfoFilteringDist
std::vector<NumericMatrix> BackInfoFilteringDist(NumericMatrix init_dsts, NumericMatrix init_srcs, std::vector<double> init_log_probs, double gps_trunc_alpha, std::vector<double> obs_lons, std::vector<double> obs_lats, std::vector<double> obs_semi_majors, std::vector<double> obs_semi_minors, std::vector<double> obs_orientations, std::vector<double> obs_depths, std::vector<double> lon_gridvals, std::vector<double> lat_gridvals, std::vector<double> surface_heights, double min_elevation, double max_elevation, double log_self_tx, double betaAR, std::vector<unsigned int> pred_steps);
RcppExport SEXP _dsmovetools2d_BackInfoFilteringDist(SEXP init_dstsSEXP, SEXP init_srcsSEXP, SEXP init_log_probsSEXP, SEXP gps_trunc_alphaSEXP, SEXP obs_lonsSEXP, SEXP obs_latsSEXP, SEXP obs_semi_majorsSEXP, SEXP obs_semi_minorsSEXP, SEXP obs_orientationsSEXP, SEXP obs_depthsSEXP, SEXP lon_gridvalsSEXP, SEXP lat_gridvalsSEXP, SEXP surface_heightsSEXP, SEXP min_elevationSEXP, SEXP max_elevationSEXP, SEXP log_self_txSEXP, SEXP betaARSEXP, SEXP pred_stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type init_dsts(init_dstsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type init_srcs(init_srcsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type init_log_probs(init_log_probsSEXP);
    Rcpp::traits::input_parameter< double >::type gps_trunc_alpha(gps_trunc_alphaSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_lons(obs_lonsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_lats(obs_latsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_semi_majors(obs_semi_majorsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_semi_minors(obs_semi_minorsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_orientations(obs_orientationsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type obs_depths(obs_depthsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type lon_gridvals(lon_gridvalsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type lat_gridvals(lat_gridvalsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type surface_heights(surface_heightsSEXP);
    Rcpp::traits::input_parameter< double >::type min_elevation(min_elevationSEXP);
    Rcpp::traits::input_parameter< double >::type max_elevation(max_elevationSEXP);
    Rcpp::traits::input_parameter< double >::type log_self_tx(log_self_txSEXP);
    Rcpp::traits::input_parameter< double >::type betaAR(betaARSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type pred_steps(pred_stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(BackInfoFilteringDist(init_dsts, init_srcs, init_log_probs, gps_trunc_alpha, obs_lons, obs_lats, obs_semi_majors, obs_semi_minors, obs_orientations, obs_depths, lon_gridvals, lat_gridvals, surface_heights, min_elevation, max_elevation, log_self_tx, betaAR, pred_steps));
    return rcpp_result_gen;
END_RCPP
}
// log_sum_c
double log_sum_c(std::vector<double> x);
RcppExport SEXP _dsmovetools2d_log_sum_c(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(log_sum_c(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dsmovetools2d_TestCTDS2DDomainIO", (DL_FUNC) &_dsmovetools2d_TestCTDS2DDomainIO, 6},
    {"_dsmovetools2d_LogTxProbs", (DL_FUNC) &_dsmovetools2d_LogTxProbs, 8},
    {"_dsmovetools2d_LogTxProbsElevation", (DL_FUNC) &_dsmovetools2d_LogTxProbsElevation, 10},
    {"_dsmovetools2d_GpsLikEval", (DL_FUNC) &_dsmovetools2d_GpsLikEval, 9},
    {"_dsmovetools2d_FF_DTMC", (DL_FUNC) &_dsmovetools2d_FF_DTMC, 9},
    {"_dsmovetools2d_SattagFilteredLL", (DL_FUNC) &_dsmovetools2d_SattagFilteredLL, 17},
    {"_dsmovetools2d_SattagPredDist", (DL_FUNC) &_dsmovetools2d_SattagPredDist, 18},
    {"_dsmovetools2d_BackInfoFilteringDist", (DL_FUNC) &_dsmovetools2d_BackInfoFilteringDist, 18},
    {"_dsmovetools2d_log_sum_c", (DL_FUNC) &_dsmovetools2d_log_sum_c, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_dsmovetools2d(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
