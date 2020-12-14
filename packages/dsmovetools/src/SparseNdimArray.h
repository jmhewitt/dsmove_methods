//
// Created by Joshua Hewitt on 12/14/20.
//

#ifndef DSMOVETOOLS_SPARSENDIMARRAY_H
#define DSMOVETOOLS_SPARSENDIMARRAY_H


template <typename Index, typename ValueType, typename Storage>
class SparseNdimArray {

    public:
        Storage data;

        void set(const Index&, const ValueType&);
        ValueType get(const Index&);

};

template <typename Index, typename ValueType, typename Storage>
void SparseNdimArray<Index, ValueType, Storage>::set(const Index &i, const ValueType &v) {
    data[i] = v;
}

template <typename Index, typename ValueType, typename Storage>
ValueType SparseNdimArray<Index, ValueType, Storage>::get(const Index &i) {
    return data[i];
}

#endif //DSMOVETOOLS_SPARSENDIMARRAY_H
