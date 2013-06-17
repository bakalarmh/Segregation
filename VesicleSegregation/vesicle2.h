//
//  vesicle.h
//  Segregation
//
//  Created by Matthew Bakalar on 2/1/13.
//  Copyright (c) 2013 Matthew Bakalar. All rights reserved.
//

#ifndef __Segregation__vesicle__
#define __Segregation__vesicle__

#include <iostream>
#include <vector>
#include <gsl/gsl_rng.h>
#include "lattice2.h"

class Vesicle {
public:
    
    Lattice2<int> lattice;
    
    int size;
    
    // Scratch vector
    std::vector<Point2<int> > molecules;
    Point2<int> trial_source;
    Point2<int> trial_dest;
    int trial_idx;
    
    // GSL library random number generator components
    const gsl_rng_type * ranT;
    gsl_rng * ranr;
    
    Vesicle(int size, int contact_size);
    
    // Print to XYZ trajectory file
    std::ostream& PrintXYZ(std::ostream& os, std::vector<std::string> atoms);
    
    // Collects non-zero lattice positions into the scratch array. Returns the number of elements collected
    int CollectParticlesScratch();
    
    // Generate a random walk move for an individual particle
    void GenerateTrialStep();
    int TestTrialStepWithInterface(int radius, float U);
    int CountParticlesInside(int radius, int type);
    void AcceptTrialMove();
    
    int& operator[] (const int nIndex) {
        return lattice[nIndex];
    }
    
    void SeedRandParticles(std::vector<int> type, std::vector<int> number);
    int MoveRandomParticle();

};

#endif /* defined(__Segregation__vesicle__) */
