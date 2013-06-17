//
//  vesicle.h
//  VesicleSegregation
//
//  Created by Matthew Bakalar on 5/31/13.
//  Copyright (c) 2013 Matthew Bakalar. All rights reserved.
//

#ifndef __VesicleSegregation__vesicle__
#define __VesicleSegregation__vesicle__

#include <iostream>
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include "lattice2.h"

class Vesicle {
public:
    
    // Vesicle parameters
    int contactsz;
    int edgelen;
    int numtypes; // Maximum number of particle types that this vesicle supports
    
    // Vesicle state
    Lattice2<int> lattice;
    Lattice2<int> molnumlattice;
    
    // Molecules are named by their position in this molecule list
    std::vector<Point2<int>> molecules;
    // Bind list holds whether or not a particle is bound
    std::vector<int> bindstate;
    // Number of molecules in contact for a given type
    std::vector<int> numContact;
    std::vector<int> numBound;

    // Constructor
    Vesicle(int edgelen, int contactsz, int numtypes);
        
    // Returns true if the site at idx is a contact site
    bool contactSite(Point2<int> idx);
    
    // Transform the index of a point on the vesicle to an index wihtin the contact region. Should only be called
    // on points that are within the contact region
    Point2<int> contactIndex(Point2<int> idx);
    
    // Transform the index of a point within the contact region back to a global index on the vesicle lattice.
    Point2<int> globalIndex(Point2<int> idx);
    
    // Returns the number of free sites left on the vesicle lattice
    int freeSites();
    
    // Returns the total number of molecules on the veislce
    int numParticles();
    
    // Print vesicle grid
    std::ostream& printGrid(std::ostream& os);
    std::ostream& printMolnumGrid(std::ostream& os);
    
    // Print contact region grid
    std::ostream& printContactGrid(std::ostream& os);
    std::ostream& printMolecules(std::ostream& os, std::vector<std::string> atoms);

};

#endif /* defined(__VesicleSegregation__vesicle__) */
