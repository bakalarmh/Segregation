//
//  vesicledimer.h
//  VesicleSegregation
//
//  Created by Matthew Bakalar on 6/2/13.
//  Copyright (c) 2013 Matthew Bakalar. All rights reserved.
//

#ifndef __VesicleSegregation__vesicledimer__
#define __VesicleSegregation__vesicledimer__

#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include "vesicle.h"

class VesicleDimer {
public:
    
    // Dimer parameters
    bool crowdAcross;
    // List of tuples, where each tuple defines a crowding relationship between (a,b)
    std::vector<std::tuple<int,int> > crowdPairs;
    
    // Binding types
    std::vector<bool> typeBinds;
    
    // Energy for each type to reside within contact
    std::vector<float> typeContactEnergy;
    // Energy for bound particle within the contact zone
    std::vector<float> typeBindingEnergy;
    
    // Pointer to the vesicles that form the dimer
    Vesicle& v0;
    Vesicle& v1;
    
    // Constructor
    VesicleDimer(Vesicle& v0, Vesicle& v1, std::vector<std::tuple<int,int> > crowdPairs, std::vector<bool> typeBinds, std::vector<float> typeContactEnergy, std::vector<float> typeBindingEnergy);
    
    // Energy of the vesicle system
    float systemEnergy();
    float testBrownianStepEnergy(int vesidx, Point2<int> loc, Point2<int> newloc);
    float testBindEnergy(int vesidx, Point2<int> loc);
    float testUnbindEnergy(int vesidx, Point2<int> loc);
    
    // Returns 1 if the site if open and available to type, -1 if the site is occupied or crowded
    int siteFree(int vesidx, Point2<int> idx, int type);
    
    // Returns the index of a valid binding partner on the opposite lattice
    std::vector<Point2<int>> bindingPartners(int vesidx, Point2<int> idx);
    std::vector<Point2<int>> unbindingPartners(int vesidx, Point2<int> idx);

    
    // Returns the number of sites available to particle of type
    int numFreeSites(int vesidx, int type);
    
    // Printing
    std::ostream& printContactMolecules(std::ostream& os, std::vector<std::string> atoms, float zoffset, float bulkoffset);
};

#endif /* defined(__VesicleSegregation__vesicledimer__) */