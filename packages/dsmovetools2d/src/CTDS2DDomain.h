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
    CTDS2DState // spatial structure
                *nbr_to[4] = {nullptr, nullptr, nullptr, nullptr},
                *nbr_from[4] = {nullptr, nullptr, nullptr, nullptr},
                // linked list to iterate over non-zero log_prob_a/b entries
                *prev_written_state_a = nullptr,
                *prev_written_state_b = nullptr;
};

struct CTDS2DStateFilter {
    virtual bool isValid(const CTDS2DState& state) { return true; }
};

/**
 * Filter to remove states that have surface_height values above or below the
 * range set by the min and max elevation arguments
 */
struct CTDS2DStateElevationFilter : public CTDS2DStateFilter {
private:
    double min_elevation, max_elevation;
public:
    CTDS2DStateElevationFilter(double min, double max) : min_elevation(min),
        max_elevation(max) { }
    bool isValid(const CTDS2DState& state) {
        if(state.surface_height > max_elevation) {
            return false;
        }
        if(state.surface_height < min_elevation) {
            return false;
        }
        return true;
    }
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
            CTDS2DState* last_state_written = nullptr;
            virtual void set(CTDS2DState& x, double v, unsigned int age) { }
            virtual void add(CTDS2DState& x, double v, unsigned int age) { }
            virtual void scale(CTDS2DState& x, double v, unsigned int age) { }
            virtual double get(const CTDS2DState& x, unsigned int age) { }
            virtual CTDS2DState* getPrevStateWritten(CTDS2DState* x) { }
    };

    // handler to update log_prob_a values in CTDS2DState objects
    class ModProbA : public ProbMods {
        public:
            void set(CTDS2DState& x, double v, unsigned int age) {
                if(age != x.prob_age_a) {
                    x.prev_written_state_a = last_state_written;
                    last_state_written = &x;
                }
                ProbMods::set(x.log_prob_a, x.prob_age_a, v, age);
            }
            void add(CTDS2DState& x, double v, unsigned int age) {
                if(age != x.prob_age_a) {
                    x.prev_written_state_a = last_state_written;
                    last_state_written = &x;
                }
                ProbMods::add(x.log_prob_a, x.prob_age_a, v, age);
            }
            void scale(CTDS2DState& x, double v, unsigned int age) {
                ProbMods::scale(x.log_prob_a, x.prob_age_a, v, age);
            }
            double get(const CTDS2DState& x, unsigned int age) {
                return ProbMods::get(x.log_prob_a, x.prob_age_a, age);
            }
            CTDS2DState* getPrevStateWritten(CTDS2DState* x) {
                return x->prev_written_state_a;
            }
    };

    // handler to update log_prob_b values in CTDS2DState objects
    class ModProbB : public ProbMods {
        public:
            void set(CTDS2DState& x, double v, unsigned int age) {
                if(age != x.prob_age_b) {
                    x.prev_written_state_b = last_state_written;
                    last_state_written = &x;
                }
                ProbMods::set(x.log_prob_b, x.prob_age_b, v, age);
            }
            void add(CTDS2DState& x, double v, unsigned int age) {
                if(age != x.prob_age_b) {
                    x.prev_written_state_b = last_state_written;
                    last_state_written = &x;
                }
                ProbMods::add(x.log_prob_b, x.prob_age_b, v, age);
            }
            void scale(CTDS2DState& x, double v, unsigned int age) {
                ProbMods::scale(x.log_prob_b, x.prob_age_b, v, age);
            }
            double get(const CTDS2DState& x, unsigned int age) {
                return ProbMods::get(x.log_prob_b, x.prob_age_b, age);
            }
            CTDS2DState* getPrevStateWritten(CTDS2DState* x) {
                return x->prev_written_state_b;
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
        // reset tail element of linked lists over nonzero sparse vector entries
        active_tgt->last_state_written = nullptr;
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

    // remove pointers to/from CTDS2DState objects with well_defined == false
    void updateConnections();

    // remove states from domain using rules in provided filtering class
    void filterStates(CTDS2DStateFilter &stateFilter);

    // flatten active non-zero entries to matrix format
    Rcpp::NumericMatrix toNumericMatrix();
    // flatten active (true) or inactive (false) non-zero entries to matrix fmt
    Rcpp::NumericMatrix toNumericMatrix(bool active);

    using CTDS2DStateType = std::vector<CTDS2DState>;

    struct Iterator {
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = CTDS2DState;
        using pointer = value_type*;
        using reference = value_type&;

        Iterator(pointer ptr, ProbMods *pmod) : m_ptr(ptr), pmod_ptr(pmod) { }

        reference operator*() const { return *m_ptr; }
        pointer operator->() { return m_ptr; }

        // Prefix increment
        Iterator& operator++() {
            m_ptr = pmod_ptr->getPrevStateWritten(m_ptr);
            return *this;
        }

        friend bool operator== (const Iterator& a, const Iterator& b) {
            return a.m_ptr == b.m_ptr;
        };
        friend bool operator!= (const Iterator& a, const Iterator& b) {
            return a.m_ptr != b.m_ptr;
        };

    private:
        pointer m_ptr;
        ProbMods* pmod_ptr;
    };

//    CTDS2DStateType::iterator begin() { return states.begin(); }
//    CTDS2DStateType::iterator end() { return states.end(); }

    Iterator begin() { return begin(false); }
    Iterator end() { return end(false); }

    Iterator begin(bool active) {
        if(active) {
            return Iterator(active_tgt->last_state_written, active_tgt);
        } else {
            return Iterator(inactive_tgt->last_state_written, inactive_tgt);
        }
    }

    Iterator end(bool active) {
        if(active) {
            return Iterator(nullptr, active_tgt);
        } else {
            return Iterator(nullptr, inactive_tgt);
        }
    }

};


#endif //DSMOVETOOLS2D_CTDS2DDOMAIN_H
