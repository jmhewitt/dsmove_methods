//
// Created by Joshua Hewitt on 12/14/20.
//

#ifndef DSMOVETOOLS_ROOK_HEADING_H
#define DSMOVETOOLS_ROOK_HEADING_H

/**
 * A unit vector that follows Rook adjacencies can be uniquely identified by 
 * the dimension and direction in which the vector points.  This struct
 * compactly represents such information.
 */
template<typename size_type>
struct RookHeading {
    
    size_type dim;
    bool lwr;
    
    RookHeading() { };
    RookHeading(size_type d, bool l) : dim(d), lwr(l) { };
    
};

#endif //DSMOVETOOLS_ROOK_HEADING_H
