//
//  vesiclemanager.h
//  VesicleSegregation
//
//  Created by Matthew Bakalar on 6/2/13.
//  Copyright (c) 2013 Matthew Bakalar. All rights reserved.
//

#ifndef __VesicleSegregation__vesiclemanager__
#define __VesicleSegregation__vesiclemanager__

#include <iostream>
#include <vector>
#include "vesicledimer.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class VesicleManager {
public:
    
    // GSL library random number generator components
    const gsl_rng_type * ranT;
    gsl_rng * ranr;
    
    // Constructor. Set up random number generator
    VesicleManager();
    
    // Place particles down on the lattice. Returns 1 on success, -1 on failure. Method
    // will fail if there are not enough free sites remaining to place particles. Particles
    // are categorized by an integer type, where type=0 is reserved for empty lattice sites.
    int SeedParticles(VesicleDimer& dimer, int vesidx, int type, int number);
    int SeedParticlesCarefully(VesicleDimer& dimer, int vesidx, int type, int number);
    
    // Monte Carlo logic
    bool ShouldAccept(float dU);
    
    // Monte Carlo moves
    int AttemptMove(VesicleDimer& dimer);
    
    int AttemptDiffuseMove(VesicleDimer& dimer, int vesidx, int molnum, Point2<int> loc);
    int AttemptBindMove(VesicleDimer& dimer, int vesidx, int molnum, Point2<int> loc, std::vector<Point2<int>> partners);
    int AttemptUnbindMove(VesicleDimer& dimer, int vesidx, int molnum, Point2<int> loc, std::vector<Point2<int>> partners);
};

#endif /* defined(__VesicleSegregation__vesiclemanager__) */
