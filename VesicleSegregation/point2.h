//
//  point2.h
//  SMBps1
//
//  Created by Matthew Bakalar on 1/27/13.
//  Copyright (c) 2013 Matthew Bakalar. All rights reserved.
//

#ifndef __SMBps1__point2__
#define __SMBps1__point2__

#include <iostream>
#include "point.h"

template <class T> class Point2 : public Point<T,2> {
public:
    Point2() : Point<T,2>() { };
    
    Point2(const Point2<T> &p) : Point<T,2>(p) {
        for (int i=0; i<2; i++)
            this->mat[i] = p.mat[i];
    }
    
    Point2( const T px, const T py)
    {
        this->mat[0]=px;
        this->mat[1]=py;
    }
    
    // Print a human readable lattice lattice to an output stream.
    friend std::ostream& operator<<(std::ostream& os, Point2<T>& pt) {
        os << pt.mat[0] << " " << pt.mat[1];
        return os;
    }
};

#endif /* defined(__SMBps1__point2__) */
