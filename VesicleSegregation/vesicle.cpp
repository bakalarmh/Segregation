//
//  vesicle.cpp
//  VesicleSegregation
//
//  Created by Matthew Bakalar on 5/31/13.
//  Copyright (c) 2013 Matthew Bakalar. All rights reserved.
//

#include "vesicle.h"
#include <iomanip>

using namespace std;

Vesicle::Vesicle(int iedgelen, int icontactsz, int inumtypes) : edgelen(iedgelen), contactsz(icontactsz), numtypes(inumtypes), molecules(vector<Point2<int>>()), lattice(Lattice2<int>(edgelen)), molnumlattice(Lattice2<int>(edgelen)), numContact(vector<int>(inumtypes+1)), numBound(vector<int>(inumtypes+1)) {
    
    numContact[0] = contactsz;

}

bool Vesicle::contactSite(Point2<int> idx) {
    return ( (idx[0] >= (edgelen/2 - contactsz/2)) && (idx[0] <= (edgelen/2 + contactsz/2)) &&
            (idx[1] >= (edgelen/2 - contactsz/2)) && (idx[1] <= (edgelen/2 + contactsz/2)) );
}

Point2<int> Vesicle::contactIndex(Point2<int> idx) {
    int x = idx[0] - ((edgelen/2) - (contactsz/2));
    int y = idx[1] - ((edgelen/2) - (contactsz/2));
    return Point2<int>(x,y);
}

Point2<int> Vesicle::globalIndex(Point2<int> idx) {
    int x = idx[0] + ((edgelen/2) - (contactsz/2));
    int y = idx[1] + ((edgelen/2) - (contactsz/2));
    return Point2<int>(x,y);
}

int Vesicle::freeSites() {
    int mols = (int)molecules.size();
    int sites = edgelen*edgelen;
    return sites - mols;
}

int Vesicle::numParticles() {
    return (int)molecules.size();
}

ostream& Vesicle::printGrid(ostream& os) {
    for(int y=0; y<edgelen; y++) {
        for(int x=0; x<edgelen; x++) {
            os << lattice[x + y*edgelen] << " ";
        }
        os << endl;
    }
    return os;
}

ostream& Vesicle::printMolnumGrid(ostream& os) {
    for(int y=0; y<edgelen; y++) {
        for(int x=0; x<edgelen; x++) {
            os << molnumlattice[x + y*edgelen] << " ";
        }
        os << endl;
    }
    return os;
}

ostream& Vesicle::printContactGrid(ostream& os) {
    for(int y=0; y<edgelen; y++) {
        for(int x=0; x<edgelen; x++) {
            if(contactSite(Point2<int>(x,y))) {
                os << lattice[x + y*edgelen] << " ";
            }
        }
        os << endl;
    }
    return os;
}

ostream& Vesicle::printMolecules(ostream& os, vector<string> atoms) {
    long count = molecules.size();
    os << "\t" << count << "\n";
    os << "Title Line \n";
    for(int i=0; i<count; i++) {
        float x = (float)molecules[i][0];
        float y = (float)molecules[i][1];
        int type = lattice[x + y*edgelen];
        if(type > 0) {
            os.setf(ios::fixed);
            os.precision(2);
            os.setf(ios::fixed);
            os.precision(10);
            os << atoms[type-1] << setw(15) << x << setw(15) << y << setw(15) << (float)0.0 << "\n";
        }
    }
    return os;
}