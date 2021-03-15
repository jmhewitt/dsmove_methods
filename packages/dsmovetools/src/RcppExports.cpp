// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/dsmovetools.h"
#include <Rcpp.h>

using namespace Rcpp;

// TestRookOrientation
std::vector<unsigned int> TestRookOrientation(std::vector<unsigned int> head, std::vector<unsigned int> tail);
RcppExport SEXP _dsmovetools_TestRookOrientation(SEXP headSEXP, SEXP tailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type head(headSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type tail(tailSEXP);
    rcpp_result_gen = Rcpp::wrap(TestRookOrientation(head, tail));
    return rcpp_result_gen;
END_RCPP
}
// TestRookDot
double TestRookDot(std::vector<unsigned int> head, std::vector<unsigned int> tail, std::vector<unsigned int> nextHead);
RcppExport SEXP _dsmovetools_TestRookDot(SEXP headSEXP, SEXP tailSEXP, SEXP nextHeadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type head(headSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type tail(tailSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type nextHead(nextHeadSEXP);
    rcpp_result_gen = Rcpp::wrap(TestRookDot(head, tail, nextHead));
    return rcpp_result_gen;
END_RCPP
}
// TestRookNeighborhood
NumericMatrix TestRookNeighborhood(std::vector<int> dims, std::vector<int> x);
RcppExport SEXP _dsmovetools_TestRookNeighborhood(SEXP dimsSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(TestRookNeighborhood(dims, x));
    return rcpp_result_gen;
END_RCPP
}
// TestSparseNdimArrayReadWrite
NumericMatrix TestSparseNdimArrayReadWrite(NumericMatrix coords, NumericVector values);
RcppExport SEXP _dsmovetools_TestSparseNdimArrayReadWrite(SEXP coordsSEXP, SEXP valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type values(valuesSEXP);
    rcpp_result_gen = Rcpp::wrap(TestSparseNdimArrayReadWrite(coords, values));
    return rcpp_result_gen;
END_RCPP
}
// TxModelParams
NumericMatrix TxModelParams(std::vector<unsigned int> cur_loc, std::vector<unsigned int> prev_loc, std::vector<unsigned int> dims, double betaAR);
RcppExport SEXP _dsmovetools_TxModelParams(SEXP cur_locSEXP, SEXP prev_locSEXP, SEXP dimsSEXP, SEXP betaARSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type cur_loc(cur_locSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type prev_loc(prev_locSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< double >::type betaAR(betaARSEXP);
    rcpp_result_gen = Rcpp::wrap(TxModelParams(cur_loc, prev_loc, dims, betaAR));
    return rcpp_result_gen;
END_RCPP
}
// TxModelSample
std::vector<unsigned int> TxModelSample(std::vector<unsigned int> cur_loc, std::vector<unsigned int> prev_loc, std::vector<unsigned int> dims, double betaAR);
RcppExport SEXP _dsmovetools_TxModelSample(SEXP cur_locSEXP, SEXP prev_locSEXP, SEXP dimsSEXP, SEXP betaARSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type cur_loc(cur_locSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type prev_loc(prev_locSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< double >::type betaAR(betaARSEXP);
    rcpp_result_gen = Rcpp::wrap(TxModelSample(cur_loc, prev_loc, dims, betaAR));
    return rcpp_result_gen;
END_RCPP
}
// TxModelLd
double TxModelLd(std::vector<unsigned int> cur_loc, std::vector<unsigned int> prev_loc, std::vector<unsigned int> dims, double betaAR, std::vector<unsigned int> dst_loc);
RcppExport SEXP _dsmovetools_TxModelLd(SEXP cur_locSEXP, SEXP prev_locSEXP, SEXP dimsSEXP, SEXP betaARSEXP, SEXP dst_locSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type cur_loc(cur_locSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type prev_loc(prev_locSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< double >::type betaAR(betaARSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type dst_loc(dst_locSEXP);
    rcpp_result_gen = Rcpp::wrap(TxModelLd(cur_loc, prev_loc, dims, betaAR, dst_loc));
    return rcpp_result_gen;
END_RCPP
}
// TestZConstrainedRookNeighborhood
NumericMatrix TestZConstrainedRookNeighborhood(std::vector<unsigned int> dims, std::vector<unsigned int> x, std::vector<double> zfield, std::vector<double> zvals);
RcppExport SEXP _dsmovetools_TestZConstrainedRookNeighborhood(SEXP dimsSEXP, SEXP xSEXP, SEXP zfieldSEXP, SEXP zvalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type zfield(zfieldSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type zvals(zvalsSEXP);
    rcpp_result_gen = Rcpp::wrap(TestZConstrainedRookNeighborhood(dims, x, zfield, zvals));
    return rcpp_result_gen;
END_RCPP
}
// TestFFRW
NumericMatrix TestFFRW(NumericMatrix a0coords, NumericVector a0values, std::vector<unsigned int> dims, int steps);
RcppExport SEXP _dsmovetools_TestFFRW(SEXP a0coordsSEXP, SEXP a0valuesSEXP, SEXP dimsSEXP, SEXP stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type a0coords(a0coordsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a0values(a0valuesSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< int >::type steps(stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(TestFFRW(a0coords, a0values, dims, steps));
    return rcpp_result_gen;
END_RCPP
}
// TestFFRWLight
NumericMatrix TestFFRWLight(NumericMatrix a0coords, NumericVector a0values, std::vector<unsigned int> dims, int steps);
RcppExport SEXP _dsmovetools_TestFFRWLight(SEXP a0coordsSEXP, SEXP a0valuesSEXP, SEXP dimsSEXP, SEXP stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type a0coords(a0coordsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a0values(a0valuesSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< int >::type steps(stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(TestFFRWLight(a0coords, a0values, dims, steps));
    return rcpp_result_gen;
END_RCPP
}
// TestFFRWLightLog
NumericMatrix TestFFRWLightLog(NumericMatrix a0coords, NumericVector log_a0values, std::vector<unsigned int> dims, int steps);
RcppExport SEXP _dsmovetools_TestFFRWLightLog(SEXP a0coordsSEXP, SEXP log_a0valuesSEXP, SEXP dimsSEXP, SEXP stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type a0coords(a0coordsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type log_a0values(log_a0valuesSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< int >::type steps(stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(TestFFRWLightLog(a0coords, log_a0values, dims, steps));
    return rcpp_result_gen;
END_RCPP
}
// FFRWLightLogConstrained
NumericMatrix FFRWLightLogConstrained(NumericMatrix a0coords, NumericVector log_a0values, std::vector<unsigned int> dims, int steps, std::vector<double> surface_heights, std::vector<double> domain_heights);
RcppExport SEXP _dsmovetools_FFRWLightLogConstrained(SEXP a0coordsSEXP, SEXP log_a0valuesSEXP, SEXP dimsSEXP, SEXP stepsSEXP, SEXP surface_heightsSEXP, SEXP domain_heightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type a0coords(a0coordsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type log_a0values(log_a0valuesSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< int >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type surface_heights(surface_heightsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type domain_heights(domain_heightsSEXP);
    rcpp_result_gen = Rcpp::wrap(FFRWLightLogConstrained(a0coords, log_a0values, dims, steps, surface_heights, domain_heights));
    return rcpp_result_gen;
END_RCPP
}
// FFRWLogConstrained
std::vector<NumericMatrix> FFRWLogConstrained(NumericMatrix a0coords, NumericVector log_a0values, std::vector<unsigned int> dims, int steps, std::vector<double> surface_heights, std::vector<double> domain_heights);
RcppExport SEXP _dsmovetools_FFRWLogConstrained(SEXP a0coordsSEXP, SEXP log_a0valuesSEXP, SEXP dimsSEXP, SEXP stepsSEXP, SEXP surface_heightsSEXP, SEXP domain_heightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type a0coords(a0coordsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type log_a0values(log_a0valuesSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< int >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type surface_heights(surface_heightsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type domain_heights(domain_heightsSEXP);
    rcpp_result_gen = Rcpp::wrap(FFRWLogConstrained(a0coords, log_a0values, dims, steps, surface_heights, domain_heights));
    return rcpp_result_gen;
END_RCPP
}
// FFRWLogConstrainedDst
std::vector<NumericMatrix> FFRWLogConstrainedDst(NumericMatrix a0coords, NumericMatrix dstcoords, NumericVector log_a0values, std::vector<unsigned int> dims, unsigned int steps, unsigned int max_steps, std::vector<double> surface_heights, std::vector<double> domain_heights);
RcppExport SEXP _dsmovetools_FFRWLogConstrainedDst(SEXP a0coordsSEXP, SEXP dstcoordsSEXP, SEXP log_a0valuesSEXP, SEXP dimsSEXP, SEXP stepsSEXP, SEXP max_stepsSEXP, SEXP surface_heightsSEXP, SEXP domain_heightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type a0coords(a0coordsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dstcoords(dstcoordsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type log_a0values(log_a0valuesSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_steps(max_stepsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type surface_heights(surface_heightsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type domain_heights(domain_heightsSEXP);
    rcpp_result_gen = Rcpp::wrap(FFRWLogConstrainedDst(a0coords, dstcoords, log_a0values, dims, steps, max_steps, surface_heights, domain_heights));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dsmovetools_TestRookOrientation", (DL_FUNC) &_dsmovetools_TestRookOrientation, 2},
    {"_dsmovetools_TestRookDot", (DL_FUNC) &_dsmovetools_TestRookDot, 3},
    {"_dsmovetools_TestRookNeighborhood", (DL_FUNC) &_dsmovetools_TestRookNeighborhood, 2},
    {"_dsmovetools_TestSparseNdimArrayReadWrite", (DL_FUNC) &_dsmovetools_TestSparseNdimArrayReadWrite, 2},
    {"_dsmovetools_TxModelParams", (DL_FUNC) &_dsmovetools_TxModelParams, 4},
    {"_dsmovetools_TxModelSample", (DL_FUNC) &_dsmovetools_TxModelSample, 4},
    {"_dsmovetools_TxModelLd", (DL_FUNC) &_dsmovetools_TxModelLd, 5},
    {"_dsmovetools_TestZConstrainedRookNeighborhood", (DL_FUNC) &_dsmovetools_TestZConstrainedRookNeighborhood, 4},
    {"_dsmovetools_TestFFRW", (DL_FUNC) &_dsmovetools_TestFFRW, 4},
    {"_dsmovetools_TestFFRWLight", (DL_FUNC) &_dsmovetools_TestFFRWLight, 4},
    {"_dsmovetools_TestFFRWLightLog", (DL_FUNC) &_dsmovetools_TestFFRWLightLog, 4},
    {"_dsmovetools_FFRWLightLogConstrained", (DL_FUNC) &_dsmovetools_FFRWLightLogConstrained, 6},
    {"_dsmovetools_FFRWLogConstrained", (DL_FUNC) &_dsmovetools_FFRWLogConstrained, 6},
    {"_dsmovetools_FFRWLogConstrainedDst", (DL_FUNC) &_dsmovetools_FFRWLogConstrainedDst, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_dsmovetools(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
