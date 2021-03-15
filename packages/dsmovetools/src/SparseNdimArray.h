//
// Created by Joshua Hewitt on 12/14/20.
//

#include "log_add.h"

#ifndef DSMOVETOOLS_SPARSENDIMARRAY_H
#define DSMOVETOOLS_SPARSENDIMARRAY_H


template <typename Index, typename ValueType, typename Storage>
class SparseNdimArrayBase {

    public:
        Storage data;

        void set(const Index&, const ValueType&);
        bool notNull(const Index&);
        virtual void add(const Index&, const ValueType&) = 0;
        virtual void addScaled(const Index&, const ValueType&, 
                               const ValueType&) = 0;
        virtual ValueType normalizeScale(const ValueType&) = 0;

};

template <typename Index, typename ValueType, typename Storage>
void SparseNdimArrayBase<Index, ValueType, Storage>::set(const Index &i, 
                                                         const ValueType &v) {
    data[i] = v;
}

template <typename Index, typename ValueType, typename Storage>
bool SparseNdimArrayBase<Index, ValueType, Storage>::notNull(const Index &i) {
    // value is defined for index if its location is not the end
    return data.find(i) != data.end();
}

/*
 * SparseNdimArrayBase implementation where array values are stored and
 * manipulated on a linear scale.
 */
template <typename Index, typename ValueType, typename Storage>
class SparseNdimArray : public SparseNdimArrayBase<Index, ValueType, Storage> {
public:
    void add(const Index &i, const ValueType &v) {
        SparseNdimArrayBase<Index, ValueType, Storage>::data[i] += v;
    }
    void addScaled(const Index &i, const ValueType &v, const ValueType &sc) {
        add(i, sc * v);
    }
    // invert values to be used as scaling factors
    ValueType normalizeScale(const ValueType &v) { return 1.0/v; }
};

/*
 * SparseNdimArrayBase implementation where array values are stored and
 * manipulated on a log scale. so, addition is implemented such that
 * log(c) is stored where c = a + b and the numerical inputs are log(a) and
 * log(b), and log(c) is stored in the Storage object.
 */
template <typename Index, typename ValueType, typename Storage>
class SparseNdimArrayLog : public SparseNdimArrayBase<Index, ValueType, Storage> {
public:
    void add(const Index &i, const ValueType &v) {
        /*
         * implement log(c) = log(a + b) given log(a) = v, log(b) = data[i],
         * and log(c) = data[i] (updated).  uses the identity:
         *   log(c) = log(a + b) = log(b) + log( 1 + exp(log(a) - log(b)) )
         *
         * corrections are made to account for the possibility that b = 0
         */

        // get iterators pointing to target and storage end
        auto tgt = SparseNdimArrayBase<Index, ValueType, Storage>::data.find(i);
        auto end = SparseNdimArrayBase<Index, ValueType, Storage>::data.end();

        if(tgt == end) {
            // target has probability zero, so create and store new value
            SparseNdimArrayBase<Index, ValueType, Storage>::data[i] = v;
        } else {
            tgt->second = log_add(v, tgt->second);
        }
    }
    void addScaled(const Index &i, const ValueType &v, const ValueType &sc) {
        add(i, v + sc);
    }
    ValueType normalizeScale(const ValueType &v) { return - log(v); }
};

#endif //DSMOVETOOLS_SPARSENDIMARRAY_H
