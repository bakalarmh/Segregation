//
//  Point.h
//  SMBps1
//
//  Created by Matthew Bakalar on 1/27/13.
//  Copyright (c) 2013 Matthew Bakalar. All rights reserved.
//

#ifndef __SMBps1__Point__
#define __SMBps1__Point__

#include <iostream>

template <class T, int n> class Point {

protected:
    T mat[n];
    void cp(const Point<T,n> &p)
    {
        memcpy(mat, p.mat, n*sizeof(T) );
    }
    
    void cp(const T p[n])
    {
        memcpy(mat, p, n*sizeof(T) );
    }
    
    void cp(const T &val)
    {
        for (int i=0; i<n; i++)
            mat[i] = val;
    }
    
public:
    // new constructor
    Point(const T pt[n])
    {
        cp(pt);
    }
    
    // copy constructor
    Point(const Point<T,n> &cp)
    {
        this->cp(cp);
    }
    
    // default constructor
    Point()
    {
        for (int i=0; i<n; i++)
            mat[i] = 0;
    }
    
    virtual ~Point(void) {};
    
    T& operator[] (const int nIndex) {
        return mat[nIndex];
    }
    
    //Point + some value
    Point<T,n> &operator+=(const T &val)
    {
        for (int i=0; i<n; i++)
            mat[i] += val;
        
        return (*this);
    }
    
    //sum of two Points
    Point<T,n> &operator+=(Point<T,n> &pt)
    {
        for (int i=0; i<n; i++)
            mat[i] += pt[i];
        
        return (*this);
    }
    
    //difference of two Points
    Point<T,n> &operator-=(Point<T,n> &pt)
    {
        for (int i=0; i<n; i++)
            mat[i] -= pt[i];
        
        return (*this);
    }
    
    void Print(std::ostream &output) const
    {
        for (int i=0; i<n-1; i++) {
            std::cout << mat[i] << ",";
        }
        std::cout << mat[n-1];
    }
};


#endif /* defined(__SMBps1__Point__) */
