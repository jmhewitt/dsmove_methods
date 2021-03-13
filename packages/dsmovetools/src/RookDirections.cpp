//
// Created by Joshua Hewitt on 12/14/20.
//

#include <Rcpp.h>
#include "RookDirections.h"

using namespace Rcpp;

// [[Rcpp::export]]
std::vector<unsigned int> TestRookOrientation(
    std::vector<unsigned int> head, std::vector<unsigned int> tail
) {

    RookDirections<unsigned int, std::vector<unsigned int>> rd(head.size());

    rd.setOrientation(head, tail);

    std::vector<unsigned int> o(2);
    o[0] = rd.orientation.dim;
    o[1] = rd.orientation.lwr;

    return o;
}

// [[Rcpp::export]]
double TestRookDot(
    std::vector<unsigned int> head, std::vector<unsigned int> tail,
    std::vector<unsigned int> nextHead
) {

    RookDirections<unsigned int, std::vector<unsigned int>> rd(head.size());
    RookDirections<unsigned int, std::vector<unsigned int>> rd2(head.size());

    rd.setOrientation(head, tail);
    rd2.setOrientation(nextHead, head);

    return rd.dotOrientation(rd2.orientation);
}