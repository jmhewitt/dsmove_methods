//
// Created by Joshua Hewitt on 7/14/21.
//

#include <Rcpp.h>
#include "log_add.h"

#ifndef DSMOVETOOLS2D_CTDS2DDOMAIN_H
#define DSMOVETOOLS2D_CTDS2DDOMAIN_H

struct CTDS2DState {
    double log_prob_a, log_prob_b, surface_height, lon_to, lat_to;
    unsigned int prob_age_a, prob_age_b, direction_of_movement, nnbrs;
    int lon_from_ind, lon_to_ind, lat_from_ind, lat_to_ind;
    bool well_defined;
    CTDS2DState *nbr_to[4] = {nullptr, nullptr, nullptr, nullptr},
                *nbr_from[4] = {nullptr, nullptr, nullptr, nullptr};
};

class CTDS2DDomain {

private:

    // Generic handler to update specific log-probability values
    class ProbMods {
        protected:
            void set(double &val_cur, unsigned int &age_cur, double v,
                     unsigned int age) {
                val_cur = v;
                age_cur = age;
            }
            void add(double &val_cur, unsigned int &age_cur, double v,
                     unsigned int age) {
                if(age_cur != age) {
                    set(val_cur, age_cur, v, age);
                } else {
                    val_cur = log_add(v, val_cur);
                }
            }
            void scale(double &val_cur, unsigned int &age_cur, double v,
                     unsigned int age) {
                if(age_cur == age) {
                    val_cur += v;
                }
            }
            double get(const double &val_cur, const unsigned int &age_cur,
                       const unsigned int age) {
                if(age_cur != age) {
                    // value is outdated, so return nan
                    return std::numeric_limits<double>::quiet_NaN();
                } else {
                    return val_cur;
                }
            }
        public:
            virtual void set(CTDS2DState& x, double v, unsigned int age) { }
            virtual void add(CTDS2DState& x, double v, unsigned int age) { }
            virtual void scale(CTDS2DState& x, double v, unsigned int age) { }
            virtual double get(const CTDS2DState& x, unsigned int age) { }
    };

    // handler to update log_prob_a values in CTDS2DState objects
    class ModProbA : public ProbMods {
        public:
            void set(CTDS2DState& x, double v, unsigned int age) {
                ProbMods::set(x.log_prob_a, x.prob_age_a, v, age);
            }
            void add(CTDS2DState& x, double v, unsigned int age) {
                ProbMods::add(x.log_prob_a, x.prob_age_a, v, age);
            }
            void scale(CTDS2DState& x, double v, unsigned int age) {
                ProbMods::scale(x.log_prob_a, x.prob_age_a, v, age);
            }
            double get(const CTDS2DState& x, unsigned int age) {
                return ProbMods::get(x.log_prob_a, x.prob_age_a, age);
            }
    };

    // handler to update log_prob_b values in CTDS2DState objects
    class ModProbB : public ProbMods {
        public:
            void set(CTDS2DState& x, double v, unsigned int age) {
                ProbMods::set(x.log_prob_b, x.prob_age_b, v, age);
            }
            void add(CTDS2DState& x, double v, unsigned int age) {
                ProbMods::add(x.log_prob_b, x.prob_age_b, v, age);
            }
            void scale(CTDS2DState& x, double v, unsigned int age) {
                ProbMods::scale(x.log_prob_b, x.prob_age_b, v, age);
            }
            double get(const CTDS2DState& x, unsigned int age) {
                return ProbMods::get(x.log_prob_b, x.prob_age_b, age);
            }
    };

    // CTMC states
    std::vector<CTDS2DState> states;

    // CTDS domain size
    unsigned int nlons, nlats;

    // cache age and probability handlers
    unsigned int prob_age = 1;
    ModProbA a_probs;
    ModProbB b_probs;
    ProbMods *active_tgt = &a_probs, *inactive_tgt = &b_probs;

public:

    /**
     *
     * @param lons unique list of longitudes in grid
     * @param lats unique list of latitudes in grid
     * @param surface_heights column-major list of heights at each grid point,
     *   where latitudes form the rows, and longitudes form the columns
     */
    CTDS2DDomain(std::vector<double> &lons, std::vector<double> &lats,
                 std::vector<double> &surface_heights);

    // increment age and change storage for new probabilities
    void swapActive() {
        std::swap(active_tgt, inactive_tgt);
        ++prob_age;
    }

    // dencrement age and change storage for old probabilities
    void unswapActive() {
        std::swap(active_tgt, inactive_tgt);
        --prob_age;
    }

    // set active log prob for a state in the domain
    void set(int lon_from_ind, int lat_from_ind, int lon_to_ind, int lat_to_ind,
             double log_prob);
    void set(CTDS2DState &state, double log_prob);

    // add to active log prob for a state in the domain
    void add(CTDS2DState &state, double log_prob);

    // scale active log prob for a state in the domain
    // (i.e., add "scale" to the active log prob)
    void scale(CTDS2DState &state, double scale);

    // access a state in the domain
    CTDS2DState* statePtr(int lon_from_ind, int lat_from_ind, int lon_to_ind,
                          int lat_to_ind);

    // retrieve active log-probability from a state
    double logProb(const CTDS2DState &state) {
        return active_tgt->get(state, prob_age);
    }

    // retrieve previous (i.e., inactive or cached) log-probability from a state
    double logProbCached(const CTDS2DState &state) {
        return inactive_tgt->get(state, prob_age - 1);
    }

    // flatten active non-zero entries to matrix format
    Rcpp::NumericMatrix toNumericMatrix();

    using CTDS2DStateType = std::vector<CTDS2DState>;

    CTDS2DStateType::iterator begin() { return states.begin(); }
    CTDS2DStateType::iterator end() { return states.end(); }

};


#endif //DSMOVETOOLS2D_CTDS2DDOMAIN_H
