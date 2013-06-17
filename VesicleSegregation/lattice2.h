//
//  lattice2.h
//  Segregation
//
//  Created by Matthew Bakalar on 1/31/13.
//  Copyright (c) 2013 Matthew Bakalar. All rights reserved.
//

#ifndef __Segregation__lattice2__
#define __Segregation__lattice2__

#include <iostream>
#include <vector>
#include "point2.h"

template <class T>
class Lattice2 {
public:
    
    int edgelen;
    int components;
    std::vector<T> board;
    
    Lattice2(int iedgelen): edgelen(iedgelen), board(iedgelen*iedgelen,0) {
        // Pass
    }
    
    T& operator[] (const int nIndex) {
        return board[nIndex];
    }
    
    // Print a human readable lattice lattice to an output stream.
    friend std::ostream& operator<<(std::ostream& os, Lattice2<T>& lattice) {
        int c = 0;
        os << "[" << c << "]" << std::endl;
        for(int i=0; i<lattice.edgelen; i++) {
            for(int j=0; j<lattice.edgelen; j++) {
                os << lattice.board[i + j*lattice.edgelen] << " ";
            }
            os << std::endl;
        }
        c++;
        return os;
    }
            
    int InsideSquare(Point2<T> point, Point2<T> origin, int side) {
        T dx = point[0]-origin[0];
        T dy = point[1]-origin[1];
        return (dx >= 0) && (dx < side) && (dy >= 0) && (dy < side);
    }
        
    Point2<T> PeriodicWrap(Point2<T> pt) {
        if(pt[0] < 0) {
            pt[0] += edgelen;
        }
        else if(pt[0] >= edgelen) {
            pt[0] -= edgelen;
        }
        if(pt[1] < 0) {
            pt[1] += edgelen;
        }
        else if(pt[1] >= edgelen) {
            pt[1] -= edgelen;
        }
        return pt;
    }
        
};

#endif /* defined(__Segregation__lattice2__) */
