//
// Created by Joshua Hewitt on 7/14/21.
//

#include "CTDS2DDomain.h"

using namespace Rcpp;

CTDS2DDomain::CTDS2DDomain(
    std::vector<double> &lons, std::vector<double> &lats,
    std::vector<double> &surface_heights
) {

    // grid dimensions
    nlons = lons.size();
    nlats = lats.size();

    // initialize CTDS state space
    states = std::vector<CTDS2DState>( nlons * nlats * 4 );
    auto state_it = states.begin();
    auto surface_height_it = surface_heights.begin();

    // fill as matrix with lons as columns and lats as rows, as if we were
    // visualizing the spatial grid as a matrix, rather than (lon,lat) pairs.
    //
    // the grid is filled by traversing destination locations.
    for(int lon_ind = 0; lon_ind < nlons; ++lon_ind) {

        // current destination longitude
        double lon_to = lons[lon_ind];

        for(int lat_ind = 0; lat_ind < nlats; ++lat_ind) {

            // current destination latitude and surface height
            double lat_to = lats[lat_ind];
            double surface_height = *(surface_height_it++);

            // addresses to connecting states, which can be transitioned to
            CTDS2DState *nbr_to[] = {nullptr, nullptr, nullptr, nullptr};
            unsigned int nnbrs = 0;
            if(lat_ind + 1 < nlats) {
                ++nnbrs;
                nbr_to[0] = &states[ 4 * ( lat_ind + 1 + lon_ind * nlats) + 0 ];
            }
            if(lon_ind + 1 < nlons) {
                ++nnbrs;
                nbr_to[1] = &states[ 4 * ( lat_ind + (lon_ind + 1) * nlats) + 1 ];
            }
            if(lat_ind - 1 > 0) {
                ++nnbrs;
                nbr_to[2] = &states[ 4 * ( lat_ind - 1 + lon_ind * nlats) + 2 ];
            }
            if(lon_ind - 1 > 0) {
                ++nnbrs;
                nbr_to[3] = &states[ 4 * ( lat_ind + (lon_ind - 1) * nlats) + 3 ];
            }

            // define states with N,E,S,W directions of movement to current loc.
            for(unsigned int dom = 0; dom < 4; ++dom) {
                // initialize to 0 probability mass and reset cache counter
                state_it->log_prob_a = -std::numeric_limits<double>::infinity();
                state_it->log_prob_b = -std::numeric_limits<double>::infinity();
                state_it->prob_age_a = 0;
                state_it->prob_age_b = 0;
                // destination location's spatial information
                state_it->surface_height = surface_height;
                state_it->lon_to = lon_to;
                state_it->lat_to = lat_to;
                state_it->lon_to_ind = lon_ind;
                state_it->lat_to_ind = lat_ind;
                // state orientation
                state_it->direction_of_movement = dom;
                // source location coordinates
                switch(dom) {
                    case 0:
                        state_it->lon_from_ind = lon_ind;
                        state_it->lat_from_ind = lat_ind - 1;
                        break;
                    case 1:
                        state_it->lon_from_ind = lon_ind - 1;
                        state_it->lat_from_ind = lat_ind;
                        break;
                    case 2:
                        state_it->lon_from_ind = lon_ind;
                        state_it->lat_from_ind = lat_ind + 1;
                        break;
                    case 3:
                        state_it->lon_from_ind = lon_ind + 1;
                        state_it->lat_from_ind = lat_ind;
                        break;
                    default:
                        break;
                }
                // check state bounds
                state_it->well_defined = true;
                if(state_it->lon_from_ind < 0) {
                    state_it->well_defined = false;
                }
                if(state_it->lat_from_ind < 0) {
                    state_it->well_defined = false;
                }
                if(state_it->lon_from_ind >= nlons) {
                    state_it->well_defined = false;
                }
                if(state_it->lat_from_ind >= nlats) {
                    state_it->well_defined = false;
                }
                // link to connecting states, which can be transitioned to
                state_it->nnbrs = nnbrs;
                std::memcpy(
                    &(state_it->nbr_to), &nbr_to, 4 * sizeof(CTDS2DState*)
                );

                // increment to new state object
                ++state_it;
            }
        }
    }

    // prune connections
    updateConnections();
}

CTDS2DState * CTDS2DDomain::statePtr(
        int lon_from_ind, int lat_from_ind, int lon_to_ind, int lat_to_ind
) {
    // direction of movement to associate with transition
    unsigned int direction_of_movement = 0;
    if(lat_from_ind < lat_to_ind) {
        direction_of_movement = 0;
    } else if(lon_from_ind < lon_to_ind) {
        direction_of_movement = 1;
    } else if(lat_to_ind < lat_from_ind) {
        direction_of_movement = 2;
    }  else if(lon_to_ind < lon_from_ind) {
        direction_of_movement = 3;
    }

    return &(
        states[ 4 * (lat_to_ind + lon_to_ind * nlons) + direction_of_movement ]
    );
}

void CTDS2DDomain::set(
    int lon_from_ind, int lat_from_ind, int lon_to_ind, int lat_to_ind,
    double log_prob
) {
    // state pointer
    CTDS2DState *tgt = statePtr(
        lon_from_ind, lat_from_ind, lon_to_ind, lat_to_ind
    );
    // set active value in state
    set(*tgt, log_prob);
}

void CTDS2DDomain::set(CTDS2DState &state, double log_prob) {
    active_tgt->set(state, log_prob, prob_age);
}

void CTDS2DDomain::add(CTDS2DState &state, double log_prob) {
    active_tgt->add(state, log_prob, prob_age);
}

void CTDS2DDomain::scale(CTDS2DState &state, double log_prob) {
    active_tgt->scale(state, log_prob, prob_age);
}

void CTDS2DDomain::updateConnections() {

    CTDS2DState *all_null[4] = {nullptr, nullptr, nullptr, nullptr};

    // update connections from all states based on current well_defined value
    auto state_end = states.end();
    for(auto state_it = states.begin(); state_it != state_end; ++state_it) {
        if(!state_it->well_defined) {
            // remove all connections from an ill-defined state
            std::memcpy(
                &(state_it->nbr_to), &all_null, 4 * sizeof(CTDS2DState*)
            );
            // update neighbor count
            state_it->nnbrs = 0;
        } else {
            for(unsigned int i = 0; i < 4; ++i) {
                if(state_it->nbr_to[i]) {
                    if(!state_it->nbr_to[i]->well_defined) {
                        // remove connection to state that is ill-defined
                        state_it->nbr_to[i] = nullptr;
                        // update neighbor count
                        --state_it->nnbrs;
                    }
                }
            }
        }
    }
}

void CTDS2DDomain::filterStates(CTDS2DStateFilter &stateFilter) {

    // update well_defined flag on all states
    auto state_end = states.end();
    for(auto state_it = states.begin(); state_it != state_end; ++state_it) {
        state_it->well_defined = stateFilter.isValid(*state_it);
    }

    // update connections from all states
    updateConnections();
}

NumericMatrix CTDS2DDomain::toNumericMatrix() {

    //
    // extract probabilities
    //

    std::vector<int> lon_from_ind, lon_to_ind, lat_from_ind, lat_to_ind;
    std::vector<double> lp;

    auto state_end = states.end();

    for(auto state_it = states.begin(); state_it != state_end; ++state_it) {
        double state_lp = logProb(*state_it);
        if(std::isfinite(state_lp)) {
            lon_from_ind.emplace_back(state_it->lon_from_ind);
            lat_from_ind.emplace_back(state_it->lat_from_ind);
            lon_to_ind.emplace_back(state_it->lon_to_ind);
            lat_to_ind.emplace_back(state_it->lat_to_ind);
            lp.emplace_back(state_lp);
        }
    }

    //
    // package results
    //

    NumericMatrix out = NumericMatrix(lp.size(), 5);

    for(unsigned int ind = 0; ind < out.nrow(); ++ind) {
        out(ind,0) = lon_from_ind[ind];
        out(ind,1) = lat_from_ind[ind];
        out(ind,2) = lon_to_ind[ind];
        out(ind,3) = lat_to_ind[ind];
        out(ind,4) = lp[ind];
    }

    colnames(out) = CharacterVector({
        "lon_from_ind", "lat_from_ind", "lon_to_ind", "lat_to_ind", "log_prob"
    });

    return out;
}

// [[Rcpp::export]]
NumericMatrix TestCTDS2DDomainIO(
    std::vector<double> lons, std::vector<double> lats,
    std::vector<double> surface_heights, NumericMatrix init_dsts,
    NumericMatrix init_srcs, std::vector<double> log_probs
) {

    // initialize CTMC state space probability vector
    CTDS2DDomain pvec(lons, lats, surface_heights);

    // populate vector
    for(unsigned int ind = 0; ind < init_dsts.nrow(); ++ind) {
        pvec.set(init_srcs(ind, 0), init_srcs(ind, 1), init_dsts(ind, 0),
                   init_dsts(ind, 1), log_probs[ind]);
    }

    return pvec.toNumericMatrix();
}