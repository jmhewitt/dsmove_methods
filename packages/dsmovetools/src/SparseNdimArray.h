//
// Created by Joshua Hewitt on 12/14/20.
//

#ifndef DSMOVETOOLS_SPARSENDIMARRAY_H
#define DSMOVETOOLS_SPARSENDIMARRAY_H


template <typename Index, typename ValueType, typename Storage>
class SparseNdimArrayBase {

    public:
        Storage data;

        void set(const Index&, const ValueType&);
        virtual void add(const Index&, const ValueType&) = 0;
        ValueType get(const Index&);

};

template <typename Index, typename ValueType, typename Storage>
void SparseNdimArrayBase<Index, ValueType, Storage>::set(const Index &i, const ValueType &v) {
    data[i] = v;
}

template <typename Index, typename ValueType, typename Storage>
ValueType SparseNdimArrayBase<Index, ValueType, Storage>::get(const Index &i) {
    return data[i];
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
            // log(a) - log(b), and a/b
            double x = v - tgt->second;
            double exp_x = exp(x);
            // evaluate log(c)
            if(exp_x == HUGE_VAL) { // exp_x == "Inf"
                tgt->second = v;
            } else if(exp_x == 0) { // a has negligible size relative to b
            } else {
                tgt->second += log(1 + exp_x);
            }
        }
    }
};

#endif //DSMOVETOOLS_SPARSENDIMARRAY_H
