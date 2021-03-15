//
// Created by Joshua Hewitt on 3/15/21.
//

#ifndef ROOKHEADING_H_LOG_ADD_H
#define ROOKHEADING_H_LOG_ADD_H

/*
 * implement log(c) = log(a + b) given log(a) = v, log(b) = data[i],
 * and log(c) = data[i] (updated).  uses the identity:
 *   log(c) = log(a + b) = log(b) + log( 1 + exp(log(a) - log(b)) )
 *
 * corrections are made to account for the possibility that b = 0
 */
double log_add(double, double);

#endif //ROOKHEADING_H_LOG_ADD_H
