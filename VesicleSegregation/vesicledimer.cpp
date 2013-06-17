//
//  vesicledimer.cpp
//  VesicleSegregation
//
//  Created by Matthew Bakalar on 6/2/13.
//  Copyright (c) 2013 Matthew Bakalar. All rights reserved.
//

#include <iomanip>
#include "vesicledimer.h"

using namespace std;

VesicleDimer::VesicleDimer(Vesicle& iv0, Vesicle& iv1, vector<tuple<int,int>> icrowdPairs, vector<bool> itypeBinds, vector<float> itypeContactEnergy, vector<float> itypeBindingEnergy) : v0(iv0), v1(iv1), crowdPairs(icrowdPairs), typeBinds(itypeBinds),typeContactEnergy(itypeContactEnergy), typeBindingEnergy(itypeBindingEnergy) {
    // Pass
}

float VesicleDimer::systemEnergy() {
    float totalenergy = 0.0;
    // Particles within the contact region
    for (int i=0; i<v0.numtypes; i++) {
        totalenergy += (typeContactEnergy[i] * (v0.numContact[i] - v0.numBound[i]));
        totalenergy += (typeBindingEnergy[i] * (v0.numBound[i]));
    }
    for (int i=0; i<v1.numtypes; i++) {
        totalenergy += (typeContactEnergy[i] * (v1.numContact[i] - v1.numBound[i]));
        totalenergy += (typeBindingEnergy[i] * (v1.numBound[i]));
    }
    
    return totalenergy;
}

float VesicleDimer::testBrownianStepEnergy(int vesidx, Point2<int> loc, Point2<int> newloc) {
    Vesicle* vesicle;
    if(vesidx == 0) {
        vesicle = &v0;
    }
    else {
        vesicle = &v1;
    }
    int type = vesicle->lattice[loc[0] + loc[1]*vesicle->edgelen];
    bool fromcontact = vesicle->contactSite(loc);
    bool tocontact = vesicle->contactSite(newloc);
    
    float dE;
    if(fromcontact && !tocontact) {
        dE = -typeContactEnergy[type];
    }
    else if(!fromcontact && tocontact) {
        dE = typeContactEnergy[type];
    }
    else {
        dE = 0.0;
    }
    return dE;
}

float VesicleDimer::testBindEnergy(int vesidx, Point2<int> loc) {
    Vesicle* vA;
    Vesicle* vB;
    if(vesidx == 0) {
        vA = &v0;
        vB = &v1;
    }
    else {
        vA = &v1;
        vB = &v0;
    }
    // Assume: point must be within contact site. Type of particle must be the same.
    Point2<int> contactloc = vA->contactIndex(loc);    
    int type = vA->lattice[loc[0] + loc[1]*vA->edgelen];
    
    return typeBindingEnergy[type];
}

float VesicleDimer::testUnbindEnergy(int vesidx, Point2<int> loc) {
    Vesicle* vA;
    Vesicle* vB;
    if(vesidx == 0) {
        vA = &v0;
        vB = &v1;
    }
    else {
        vA = &v1;
        vB = &v0;
    }
    // Assume: point must be within contact site. Type of particle must be the same.
    Point2<int> contactloc = vA->contactIndex(loc);
    int type = vA->lattice[loc[0] + loc[1]*vA->edgelen];
    
    return -1*typeBindingEnergy[type];
}

int VesicleDimer::numFreeSites(int vesidx, int type) {
    Vesicle* vA;
    Vesicle* vB;
    if(vesidx == 0) {
        vA = &v0;
        vB = &v1;
    }
    else {
        vA = &v1;
        vB = &v0;
    }
    
    int sitesfree = vA->freeSites();
    // Remove sites that are crowded by particles on the apposing vesicle
    for(tuple<int,int> crowdpair : crowdPairs) {
        if(get<1>(crowdpair) == type) {
            sitesfree -= vB->numContact[get<0>(crowdpair)];
        }
    }
    return sitesfree;
}

int VesicleDimer::siteFree(int vesidx, Point2<int> idx, int type) {
    Vesicle* vA;
    Vesicle* vB;
    if(vesidx == 0) {
        vA = &v0;
        vB = &v1;
    }
    else {
        vA = &v1;
        vB = &v0;
    }
    Point2<int> contactidx = vA->contactIndex(idx);
    Point2<int> bidx = vB->globalIndex(contactidx);
    int aid = vA->lattice[idx[0] + idx[1]*vA->edgelen];
    int bid = vB->lattice[bidx[0] + bidx[1]*vB->edgelen];
    
    // Assume the site is free, then check all ways it might not be available
    int sitefree = 1;
    // If A vesicle site is occupied, the site is not free
    if (aid != 0) {
        sitefree = -1;
    }
    else {
        // Does the particle on lattice B crowd type?
        for(tuple<int,int> crowdPair : crowdPairs) {
            if((get<0>(crowdPair) == bid) && (get<1>(crowdPair) == type)) {
                sitefree = -1;
            }
        }
    }
    return sitefree;
}

vector<Point2<int>> VesicleDimer::bindingPartners(int vesidx, Point2<int> loc) {
    Vesicle* vA;
    Vesicle* vB;
    if(vesidx == 0) {
        vA = &v0;
        vB = &v1;
    }
    else {
        vA = &v1;
        vB = &v0;
    }
    Point2<int> contactloc = vA->contactIndex(loc);
    Point2<int> opploc = vB->globalIndex(contactloc);
    int type = vA->lattice[loc[0] + loc[1]*vA->edgelen];
    int opptype = vB->lattice[opploc[0] + opploc[1]*vB->edgelen];
    vector<Point2<int>> options;
    
    if(type == opptype) {
        int molnum = vB->molnumlattice[opploc[0] + opploc[1]*vB->edgelen];
        int bound = vB->bindstate[molnum];
        if (bound == 0) {
            options.push_back(opploc);
        }
        else {
            cout << "Requested binding data on a bound particle!" << endl;
        }
    }
    
    return options;
}

vector<Point2<int>> VesicleDimer::unbindingPartners(int vesidx, Point2<int> loc) {
    Vesicle* vA;
    Vesicle* vB;
    if(vesidx == 0) {
        vA = &v0;
        vB = &v1;
    }
    else {
        vA = &v1;
        vB = &v0;
    }
    Point2<int> contactloc = vA->contactIndex(loc);
    Point2<int> opploc = vB->globalIndex(contactloc);
    int type = vA->lattice[loc[0] + loc[1]*vA->edgelen];
    int opptype = vB->lattice[opploc[0] + opploc[1]*vB->edgelen];
    vector<Point2<int>> options;

    if(type == opptype) {
        int molnum = vB->molnumlattice[opploc[0] + opploc[1]*vB->edgelen];
        int bound = vB->bindstate[molnum];
        if (bound == 1) {
            options.push_back(opploc);
        }
        else {
            cout << "Requested unbinding data on an unbound particle!" << endl;
        }
    }

    return options;
}

std::ostream& VesicleDimer::printContactMolecules(std::ostream& os, vector<string> atoms, float zoffset, float bulkoffset) {
    long count0 = v0.molecules.size();
    long count1 = v1.molecules.size();
    os << "\t" << (count0+count1) << "\n";
    os << "Title Line \n";
    for(int i=0; i<count0; i++) {
        float x = (float)(v0.molecules[i][0]);
        float y = (float)(v0.molecules[i][1]);
        float z;
        if(v0.contactSite(Point2<int>(x,y))) {
            z = 0.0;
        }
        else {
            z = bulkoffset;
        }
        int type = v0.lattice[x + y*v0.edgelen];
        if(type > 0) {
            os.setf(ios::fixed);
            os.precision(2);
            os.setf(ios::fixed);
            os.precision(10);
            os << atoms[type-1] << setw(15) << x << setw(15) << y << setw(15) << z << "\n";
        }
    }
    for(int i=0; i<count1; i++) {
        float x = (float)(v1.molecules[i][0]);
        float y = (float)(v1.molecules[i][1]);
        float z;
        if(v1.contactSite(Point2<int>(x,y))) {
            z = zoffset;
        }
        else {
            z = bulkoffset;
        }
        int type = v1.lattice[x + y*v1.edgelen];
        if(type > 0) {
            os.setf(ios::fixed);
            os.precision(2);
            os.setf(ios::fixed);
            os.precision(10);
            os << atoms[type-1] << setw(15) << x << setw(15) << y << setw(15) << z << "\n";
        }
    }
    return os;
}

