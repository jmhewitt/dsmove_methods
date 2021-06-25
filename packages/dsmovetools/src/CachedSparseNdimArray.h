//
// Created by Joshua Hewitt on 12/14/20.
//

// // [[Rcpp::depends(BH)]]
#include <Rcpp.h>

//#include <boost/container_hash/hash.hpp>

#include "log_add.h"

#ifndef DSMOVETOOLS_CACHEDSPARSENDIMARRAY_H
#define DSMOVETOOLS_CACHEDSPARSENDIMARRAY_H

class DualSparseCoordVec {

    public:
        //
        // Data structures
        //

        using Coord = std::vector<uint>;
        using CoordPair = std::pair<Coord, Coord>;

        using uint = unsigned int;

        //
        // Data structures and storage helpers
        //

//        struct CoordPairHash {
//            std::size_t operator()(const CoordPair &V) const {
//                std::size_t hash = boost::hash_range(V.first.begin(),
//                                                     V.first.end());
//                boost::hash_range(hash, V.second.begin(), V.second.end());
//                return hash;
//            }
//        };

        struct DualVal {
            double val_a, val_b;
            unsigned int age_a = 0, age_b = 0;
        };

//        std::unordered_map<CoordPair, DualVal, CoordPairHash> data;

        using StorageType = std::map<CoordPair, DualVal>;
        StorageType data;

    private:

        //
        // Value modification helpers
        //

        class IOTgt {
            protected:
                void add(double &val_cur, uint &age_cur, double v, uint age) {
                    if(age_cur != age) {
                        // value is outdated, so update age and store value
                        // (i.e., add as if we were starting from 0)
                        age_cur = age;
                        val_cur = v;
                    } else {
                        val_cur = log_add(v, val_cur);
                    }
                };
                void unprotectedSubtract(double &val_cur, double v) {
                    val_cur -= v;
                };
                double get(double &val_cur, uint &age_cur, uint age) {
                    if(age_cur != age) {
                        // value is outdated, so return nan
                        return std::numeric_limits<double>::quiet_NaN();
                    } else {
                        return val_cur;
                    }
                }
            public:
                virtual void add(DualVal& x, double v, uint age) { };
                virtual void unprotectedSubtract(DualVal& x, double v) { };
                virtual double get(DualVal& x, uint age) { };
        };

        class IOA : public IOTgt {
            public:
                void add (DualVal& x, double v, uint age) {
                    IOTgt::add(x.val_a, x.age_a, v, age);
                }
                void unprotectedSubtract(DualVal& x, double v) {
                    IOTgt::unprotectedSubtract(x.val_a, v);
                }
                double get (DualVal& x, uint age) {
                    return IOTgt::get(x.val_a, x.age_a, age);
                }
        };

        class IOB : public IOTgt {
            public:
                void add (DualVal& x, double v, uint age) {
                    IOTgt::add(x.val_b, x.age_b, v, age);
                }
                void unprotectedSubtract(DualVal& x, double v) {
                    IOTgt::unprotectedSubtract(x.val_b, v);
                }
                double get (DualVal& x, uint age) {
                    return IOTgt::get(x.val_b, x.age_b, age);
                }
        };

        uint age = 1;
        IOA a_io;
        IOB b_io;
        IOTgt *active_tgt = &a_io, *inactive_tgt = &b_io;

    public:

        //
        // IO
        //

        void swapActive() {
            std::swap(active_tgt, inactive_tgt);
            ++age;
        };

        void addToActive(const CoordPair &x, double v) {
            active_tgt->add(data[x], v, age);
        };

        void unprotectedSubtractFromActive(DualVal &x, double v) {
            active_tgt->unprotectedSubtract(x, v);
        };

        void unprotectedSubtractFromActive(CoordPair &x, double v) {
            unprotectedSubtractFromActive(data[x], v);
        };

        double activeValue(DualVal &x) {
            return active_tgt->get(x, age);
        };

        double activeValue(CoordPair &x) {
            return activeValue(data[x]);
        };

        double inactiveValue(DualVal &x) {
            return inactive_tgt->get(x, age - 1);
        }

        double inactiveValue(CoordPair &x) {
            return inactiveValue(data[x]);
        }

};

#endif //DSMOVETOOLS_CACHEDSPARSENDIMARRAY_H
