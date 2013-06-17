//
//  vesiclemanager.cpp
//  VesicleSegregation
//
//  Created by Matthew Bakalar on 6/2/13.
//  Copyright (c) 2013 Matthew Bakalar. All rights reserved.
//

#include <cmath>
#include <vector>
#include "vesiclemanager.h"
#include "lattice2.h"

using namespace std;

VesicleManager::VesicleManager() {
    unsigned long seed;
    seed = time(NULL);
    gsl_rng_env_setup();
    ranT = gsl_rng_default;
    ranr = gsl_rng_alloc(ranT);
    gsl_rng_set(ranr, seed);
}

int VesicleManager::SeedParticles(VesicleDimer& dimer, int vesidx, int type, int number) {
    Vesicle* vesicle;
    if(vesidx == 0) {
        vesicle = &(dimer.v0);
    }
    else {
        vesicle = &(dimer.v1);
    }
    if(dimer.numFreeSites(vesidx, type) < number) {
        return -1;
    }
    else {        
        int size = vesicle->edgelen;
        for(int i=0; i < number; i++) {
            int rx = gsl_rng_uniform(ranr)*size;
            int ry = gsl_rng_uniform(ranr)*size;
            while(dimer.siteFree(vesidx, Point2<int>(rx,ry), type) != 1) {
                rx = gsl_rng_uniform(ranr)*size;
                ry = gsl_rng_uniform(ranr)*size;
            }
            // Place particles on vesicle lattice
            Point2<int> pt(rx,ry);
            if(vesicle->contactSite(pt) == 1) {
                vesicle->numContact[type] += 1;
                vesicle->numContact[0] -= 1;
                vesicle->numBound[type] = 0; // No effect, but a reminder that we init unbound particles
            }
            vesicle->lattice[rx + ry*size] = type;
            vesicle->molecules.push_back(Point2<int>(rx,ry));
            vesicle->bindstate.push_back(0);
            vesicle->molnumlattice[rx + ry*size] = (int)(vesicle->molecules.size() - 1);
        }
        return 1;
    }
}

int VesicleManager::SeedParticlesCarefully(VesicleDimer& dimer, int vesidx, int type, int number) {
    Vesicle* vesicle;
    if(vesidx == 0) {
        vesicle = &(dimer.v0);
    }
    else {
        vesicle = &(dimer.v1);
    }
    if(dimer.numFreeSites(vesidx, type) < number) {
        return -1;
    }
    else {
        // Prepare all free spots
        vector<Point2<int>> spots;
        for(int i=0; i<vesicle->edgelen; i++) {
            for(int j=0; j<vesicle->edgelen; j++) {
                Point2<int> loc(j,i);
                if(dimer.siteFree(vesidx, loc, type)) {
                    spots.push_back(Point2<int>(j,i));
                }
            }
        }
        
        for(int i=0; i < number; i++) {
            int ridx = gsl_rng_uniform(ranr)*spots.size();
            Point2<int> pt = spots[ridx];
            // Place particles on vesicle lattice
            if(vesicle->contactSite(pt) == 1) {
                vesicle->numContact[type] += 1;
                vesicle->numContact[0] -= 1;
                vesicle->numBound[type] = 0; // No effect, but a reminder that we init unbound particles
            }
            vesicle->lattice[pt[0] + pt[1]*vesicle->edgelen] = type;
            vesicle->molecules.push_back(pt);
            vesicle->bindstate.push_back(0);
            vesicle->molnumlattice[pt[0] + pt[1]*vesicle->edgelen] = (int)(vesicle->molecules.size() - 1);
            // Remove position from free spots list
            spots.erase(spots.begin()+ridx);
        }
        return 1;
    }
}

bool VesicleManager::ShouldAccept(float dU) {
    float r1 = gsl_rng_uniform(ranr);
    return exp(-dU) > r1;
}

int VesicleManager::AttemptMove(VesicleDimer& dimer) {
    Vesicle* vA;
    Vesicle* vB;
    int vesidx;
    // Select a random molecule from one of the two vesicles
    int totalmols = dimer.v0.numParticles() + dimer.v1.numParticles();
    int molnum = gsl_rng_uniform(ranr)*totalmols;
    if(molnum < dimer.v0.numParticles()) {
        vesidx = 0;
        vA = &(dimer.v0);
        vB = &(dimer.v1);
    }
    else {
        vesidx = 1;
        vA = &(dimer.v1);
        vB = &(dimer.v0);
        molnum -= dimer.v0.numParticles();
    }
    Point2<int> loc = vA->molecules[molnum];    
    int type = vA->lattice[loc[0] + vA->edgelen*loc[1]];
    
    // If molecule is bound, its only move is to unbind. If a molecule is unbound and in contact etc, it can bind.
    // If a molecule is not in contact or not adjacent to a partner it can only diffuse.
    bool iscontact = vA->contactSite(loc);
    bool isbound = vA->bindstate[molnum];
    
    int accepted = 0;
    
    // Outside of contact zone
    if(!iscontact) {
        accepted = AttemptDiffuseMove(dimer, vesidx, molnum, loc);
    }
    // Inside of contact zone
    else {
        // Bound across the interface
        if(isbound) {
            vector<Point2<int>> unbindpartners = dimer.unbindingPartners(vesidx, loc);
            if(unbindpartners.size() == 0) {
                cout << endl;
            }
            accepted = AttemptUnbindMove(dimer, vesidx, molnum, loc, unbindpartners);
        }
        else {
            vector<Point2<int>> partners = dimer.bindingPartners(vesidx, loc);
            // Move - bind or diffuse
            if((partners.size() > 0) && dimer.typeBinds[type]) {
                int movnum = gsl_rng_uniform(ranr)*2;
                // Move - bind
                if(movnum == 0) {
                    accepted = AttemptBindMove(dimer, vesidx, molnum, loc, partners);
                }
                // Move - diffuse
                else {
                    accepted = AttemptDiffuseMove(dimer, vesidx, molnum, loc);
                }
            }
            // Move - diffuse
            else {
                accepted = AttemptDiffuseMove(dimer, vesidx, molnum, loc);
            }
        }
    }
    
    return accepted;
}

int VesicleManager::AttemptDiffuseMove(VesicleDimer& dimer, int vesidx, int molnum, Point2<int> loc) {
    Vesicle* vesicle;
    if(vesidx == 0) {
        vesicle = &(dimer.v0);
    }
    else {
        vesicle = &(dimer.v1);
    }
    Point2<int> newloc(loc[0],loc[1]);
    
    // Select a random direction for the step
    int dir = floor(gsl_rng_uniform(ranr)*4);
    switch(dir)  {
        case 0:
            newloc[0] -= 1;
            break;
        case 1:
            newloc[0] += 1;
            break;
        case 2:
            newloc[1] -= 1;
            break;
        case 3:
            newloc[1] += 1;
            break;
        default:
            // pass
            break;
    }
    // Wrap the new location with the periodic boundary conditions
    newloc = vesicle->lattice.PeriodicWrap(newloc);
    
    // Test to see if the new location is available
    int type = vesicle->lattice[loc[0] + loc[1]*vesicle->edgelen];
    bool isfree = (dimer.siteFree(vesidx, newloc, type) == 1);

    if(isfree) {
        bool accept = false;
        // Compute the energy of the brownian step move
        float dE = dimer.testBrownianStepEnergy(vesidx, loc, newloc);
        if (dE < 0) {
            accept = true;
        }
        else {
            accept = ShouldAccept(dE);
        }
        if (accept) {
            // Complete the trial move!
            vesicle->lattice[loc[0] + loc[1]*vesicle->edgelen] = 0;
            vesicle->lattice[newloc[0] + newloc[1]*vesicle->edgelen] = type;
            
            vesicle->molnumlattice[loc[0] + loc[1]*vesicle->edgelen] = 0;
            vesicle->molnumlattice[newloc[0] + newloc[1]*vesicle->edgelen] = molnum;
            
            vesicle->molecules[molnum][0] = newloc[0];
            vesicle->molecules[molnum][1] = newloc[1];
            
            if(vesicle->bindstate[molnum] == 1) {
                cout << "trying to diffuse a bound particle!" << endl;
            }
            
            
            // Update the number within the contact zone if it changed
            if( (vesicle->contactSite(loc)) && !(vesicle->contactSite(newloc)) ) {
                vesicle->numContact[type] -= 1;
                vesicle->numContact[0] += 1;
            }
            else if( !(vesicle->contactSite(loc)) && (vesicle->contactSite(newloc)) ) {
                vesicle->numContact[type] += 1;
                vesicle->numContact[0] -= 1;
            }
            
            return 1;
        }
    }
    return 0;
}

int VesicleManager::AttemptBindMove(VesicleDimer& dimer, int vesidx, int molnum, Point2<int> loc, vector<Point2<int>> partners) {
    Vesicle* vA;
    Vesicle* vB;
    if(vesidx == 0) {
        vA = &(dimer.v0);
        vB = &(dimer.v1);
    }
    else {
        vA = &(dimer.v1);
        vB = &(dimer.v0);
    }
    
    // Debug
    if(!(vA->contactSite(loc))) {
        cout << "Trying to bind outside of the contact zone!" << endl;
    }
    else if(!(vA->bindstate[molnum] == 0)) {
        cout << "Trying to bind a bound particle!" << endl;
    }

    int type = vA->lattice[loc[0] + loc[1]*vA->edgelen];
    bool accept = false;
    // Compute the energy of the brownian step move
    float dE = dimer.testBindEnergy(vesidx, loc);
    if (dE < 0) {
        accept = true;
    }
    else {
        accept = ShouldAccept(dE);
    }
    if (accept) {
        // Complete the trial move! Pick one of the potential binding partners
        int randidx = gsl_rng_uniform(ranr) * (partners.size()-1);
        Point2<int> opploc = partners[randidx];
        int oppmolnum = vB->molnumlattice[opploc[0] + opploc[1]*vB->edgelen];
        
        // Update binding state information
        vA->bindstate[molnum] = 1;
        vB->bindstate[oppmolnum] = 1;
        vA->numBound[type] += 1;
        vB->numBound[type] += 1;
        
        return 1;
    }
    else {
        return 0;
    }
}

int VesicleManager::AttemptUnbindMove(VesicleDimer& dimer, int vesidx, int molnum, Point2<int> loc, vector<Point2<int>> partners) {
    Vesicle* vA;
    Vesicle* vB;
    if(vesidx == 0) {
        vA = &(dimer.v0);
        vB = &(dimer.v1);
    }
    else {
        vA = &(dimer.v1);
        vB = &(dimer.v0);
    }
    
    int type = vA->lattice[loc[0] + loc[1]*vA->edgelen];
    bool accept = false;
    // Compute the energy of the brownian step move
    float dE = dimer.testUnbindEnergy(vesidx, loc);
    if (dE < 0) {
        accept = true;
    }
    else {
        accept = ShouldAccept(dE);
    }
    if (accept) {
        // Complete the trial move! Pick one of the potential binding partners
        int randidx = gsl_rng_uniform(ranr) * (partners.size()-1);
        Point2<int> opploc = partners[randidx];
        int oppmolnum = vB->molnumlattice[opploc[0] + opploc[1]*vB->edgelen];
        
        // Update binding state information
        vA->bindstate[molnum] = 0;
        vB->bindstate[oppmolnum] = 0;
        vA->numBound[type] -= 1;
        vB->numBound[type] -= 1;

        return 1;
    }
    else {
        return 0;
    }
}
