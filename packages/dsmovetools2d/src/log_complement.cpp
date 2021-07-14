//
// Created by Joshua Hewitt on 3/15/21.
//

#include "Rcpp.h"

/*
 * implement log(c) = log(1 - a) given log(a) and 0 < a < 1.
 */
double log_complement(double log_a) {
    return log(1 - exp(log_a));
}

