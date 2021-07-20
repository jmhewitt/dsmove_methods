//
// Created by Joshua Hewitt on 3/15/21.
//

#include "Rcpp.h"

/*
 * implement log(c) = log(a + b) given log(a) = v, log(b) = data[i],
 * and log(c) = data[i] (updated).  uses the identity:
 *   log(c) = log(a + b) = log(b) + log( 1 + exp(log(a) - log(b)) )
 *
 * corrections are made to account for the possibility that b = 0
 */
double log_add(double log_a, double log_b) {
    if(!std::isfinite(log_a)) {
        if(log_a < 0) {
            return log_b;
        }
    }
    if(!std::isfinite(log_b)) {
        if(log_b < 0) {
            return log_a;
        }
    }
    double x = log_a - log_b;
    double exp_x = exp(x);
    // evaluate log(c)
    if(exp_x == HUGE_VAL) { // exp_x == "Inf"
        return log_a;
    }
    if(exp_x == 0) { // a has negligible size relative to b
        return log_b;
    }
    return log_b + std::log1p(exp_x);
}