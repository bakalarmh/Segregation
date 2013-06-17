//
//  tests.cpp
//  VesicleSegregation
//
//  Created by Matthew Bakalar on 5/31/13.
//  Copyright (c) 2013 Matthew Bakalar. All rights reserved.
//

#include <fstream>
#include <sstream>
#include <string>
#include "tests.h"
#include "vesicle.h"
#include "vesicledimer.h"
#include "vesiclemanager.h"

using namespace std;

Tests::Tests() {
    // Pass
}

void Tests::VesicleTest() {
    int edgelen = 21;
    int contact = 5;
    int ntypes = 2;
    Vesicle vesicle(edgelen, contact, ntypes);
    for(int i=0; i<edgelen; i++) {
        for(int j=0; j<edgelen; j++) {
            if(vesicle.contactSite(Point2<int>(i,j))) {
                cout << "1 ";
            }
            else {
                cout << "0 ";
            }
        }
        cout << endl;
    }
}

void Tests::ContactTest() {
    int edgelen = 6;
    int contact = 3;
    int ntypes = 2;
    Vesicle vesicle(edgelen, contact, ntypes);
    for(int i=0; i<edgelen; i++) {
        for(int j=0; j<edgelen; j++) {
            Point2<int> pt = Point2<int>(i,j);
            if(vesicle.contactSite(pt)) {
                Point2<int> pos = vesicle.contactIndex(pt);
                cout << pos[1] << " ";
            }
            else {
                cout << ". ";
            }
        }
        cout << endl;
    }
    
    for(int i=0; i<contact; i++) {
        for(int j=0; j<contact; j++) {
            Point2<int> pt = Point2<int>(i,j);
            Point2<int> pos = vesicle.globalIndex(pt);
            cout << pos[0] << " ";
        }
        cout << endl;
    }
}

void Tests::DimerContactTest() {
    int edgelen = 9;
    int edgelen2 = 5;
    int contact = 3;
    int ntypes = 2;
    Vesicle v1(edgelen, contact, ntypes);
    Vesicle v2(edgelen2, contact, ntypes);
    
    // Define crowding relationships. 1 crowds 1, 1 crowds 2, 2 crowds 1
    vector<tuple<int,int>> crowdPairs;
    crowdPairs.push_back(tuple<int,int>(1,1));
    crowdPairs.push_back(tuple<int,int>(1,2));
    crowdPairs.push_back(tuple<int,int>(2,1));
    
    // Binding types
    const bool typebinds_arr[] = {false, true, false};
    vector<bool> typeBinds(typebinds_arr, end(typebinds_arr));
    
    // Define energy for particles in special states
    vector<float> contactEnergy(3);
    contactEnergy[0] = 0.0;
    contactEnergy[1] = 10.0;
    contactEnergy[2] = 0.0;
    vector<float> bindingEnergy(3);
    bindingEnergy[0] = 0.0;
    bindingEnergy[1] = 0.0;
    bindingEnergy[2] = 0.0;
    
    VesicleDimer dimer(v1, v2, crowdPairs, typeBinds, contactEnergy, bindingEnergy);
    VesicleManager manager;
    int seeded = manager.SeedParticles(dimer, 0, 1, 40);
    seeded = manager.SeedParticles(dimer, 1, 1, 10);
    
    v1.printGrid(cout);
    cout << endl;
    v2.printGrid(cout);
    
    v1.printContactGrid(cout);
    v2.printContactGrid(cout);
    
    cout << "N-contact = " << v1.numContact[1] << endl;
    cout << "N-bound = " << v1.numBound[1] << endl;
    cout << "N-contact = " << v2.numContact[1] << endl;
    cout << "N-bound = " << v2.numBound[1] << endl;
    
    cout << "Total system energy: " << dimer.systemEnergy() << endl;
}

void Tests::BrownianTest() {
    string outpath = "/Users/mhbakalar/Documents/Fletchlab/Data/Segsim/dual_lattice_2D/testoutput/brownmols.xyz";
    ofstream outfile;
    outfile.open(outpath);
    
    int edgelen = 1000;
    int edgelen2 = 1000;
    int contact = 100;
    int ntypes = 2;
    Vesicle v1(edgelen, contact, ntypes);
    Vesicle v2(edgelen2, contact, ntypes);
    
    // Name molecules for vmd printing
    const char* typenm_arr[] = {"C", "H"};
    const char* boundnm_arr[] = {"S", "P"};
    vector<string> atoms(typenm_arr, end(typenm_arr));
    vector<string> boundatoms(boundnm_arr, end(boundnm_arr));
    
    // Define crowding relationships. 1 crowds 1, 1 crowds 2, 2 crowds 1
    // No one will crowd across the interface today!
    vector<tuple<int,int>> crowdPairs;
    // crowdPairs.push_back(tuple<int,int>(1,1));
    // crowdPairs.push_back(tuple<int,int>(1,2));
    // crowdPairs.push_back(tuple<int,int>(2,1));
    
    // Binding types
    const bool typebinds_arr[] = {false, true, false};
    vector<bool> typeBinds(typebinds_arr, end(typebinds_arr));
    
    // Define energy for particles in special states
    vector<float> contactEnergy(3);
    contactEnergy[0] = 0.0;
    contactEnergy[1] = 0.0;
    contactEnergy[2] = 0.0;
    vector<float> bindingEnergy(3);
    bindingEnergy[0] = 0.0;
    bindingEnergy[1] = -10.0;
    bindingEnergy[2] = 0.0;
    
    VesicleDimer dimer(v1, v2, crowdPairs, typeBinds, contactEnergy, bindingEnergy);
    VesicleManager manager;
    int seeded = manager.SeedParticles(dimer, 0, 1, 1000);
    seeded = manager.SeedParticles(dimer, 1, 1, 1000);
    
    int acceptcount = 0;
    for(int sweep=0; sweep<1000; sweep++) {
        for(int mol=0; mol<100000; mol++) {
            int accepted = 0;
            while (accepted == 0) {
                accepted = manager.AttemptMove(dimer);
                acceptcount += accepted;
            }
        }
        // Print the grid
        dimer.v0.printMolecules(outfile, atoms);
    }
    outfile << endl;
    outfile.close();
}

void Tests::DimerTest() {
    int edgelen = 100;
    int contact = 25;
    int ntypes = 2;
    Vesicle v1(edgelen, contact, ntypes);
    Vesicle v2(edgelen, contact, ntypes);
    
    // Name molecules for vmd printing
    const char* typenm_arr[] = {"C", "H"};
    const char* boundnm_arr[] = {"S", "P"};
    vector<string> atoms(typenm_arr, end(typenm_arr));
    vector<string> boundatoms(boundnm_arr, end(boundnm_arr));
    
    // Define crowding relationships. 1 crowds 1, 1 crowds 2, 2 crowds 1
    vector<tuple<int,int>> crowdPairs;
    crowdPairs.push_back(tuple<int,int>(1,1));
    crowdPairs.push_back(tuple<int,int>(1,2));
    crowdPairs.push_back(tuple<int,int>(2,1));
    
    // Binding types
    const bool typebinds_arr[] = {false, true, false};
    vector<bool> typeBinds(typebinds_arr, end(typebinds_arr));
    
    // Define contact penalty
    vector<float> contactEnergy(3);
    contactEnergy[0] = 0.0;
    contactEnergy[1] = 0.0;
    contactEnergy[2] = 0.0;
    vector<float> bindingEnergy(3);
    bindingEnergy[0] = 0.0;
    bindingEnergy[1] = 0.0;
    bindingEnergy[2] = 0.0;
    
    VesicleDimer dimer(v1, v2, crowdPairs, typeBinds, contactEnergy, bindingEnergy);
    VesicleManager manager;
    int seeded = manager.SeedParticles(dimer, 0, 1, 100);
    cout << "seed worked: " << seeded << endl;
    seeded = manager.SeedParticles(dimer, 1, 1, 100);
    cout << "seed worked: " << seeded << endl;
    
    
    int count = 0;
    for(int i=0; i<edgelen; i++) {
        for(int j=0; j<edgelen; j++) {
            int spotsum = (v1.lattice[j + i*edgelen] + v2.lattice[j + i*edgelen]);
            count += spotsum;
            cout << spotsum << " ";
        }
        cout << endl;
    }
    
    cout << "count sum = " << count << endl;
    
}

void Tests::CarefulSeedTest() {
    int edgelen = 11;
    int contact = 5;
    int ntypes = 2;
    Vesicle v1(edgelen, contact, ntypes);
    Vesicle v2(edgelen, contact, ntypes);
    
    // Name molecules for vmd printing
    const char* typenm_arr[] = {"C", "H"};
    const char* boundnm_arr[] = {"S", "P"};
    vector<string> atoms(typenm_arr, end(typenm_arr));
    vector<string> boundatoms(boundnm_arr, end(boundnm_arr));
    
    // Define crowding relationships. 1 crowds 1, 1 crowds 2, 2 crowds 1
    vector<tuple<int,int>> crowdPairs;
    crowdPairs.push_back(tuple<int,int>(1,1));
    crowdPairs.push_back(tuple<int,int>(1,2));
    crowdPairs.push_back(tuple<int,int>(2,1));
    
    // Binding types
    const bool typebinds_arr[] = {false, true, false};
    vector<bool> typeBinds(typebinds_arr, end(typebinds_arr));
    
    // Define contact penalty
    vector<float> contactEnergy(3);
    contactEnergy[0] = 0.0;
    contactEnergy[1] = 0.0;
    contactEnergy[2] = 0.0;
    vector<float> bindingEnergy(3);
    bindingEnergy[0] = 0.0;
    bindingEnergy[1] = 0.0;
    bindingEnergy[2] = 0.0;
    
    VesicleDimer dimer(v1, v2, crowdPairs, typeBinds, contactEnergy, bindingEnergy);
    VesicleManager manager;
    int seeded = manager.SeedParticles(dimer, 0, 1, 50);
    cout << "seed worked: " << seeded << endl;
    seeded = manager.SeedParticles(dimer, 1, 1, 50);
    cout << "seed worked: " << seeded << endl;
    
    
    int count = 0;
    for(int i=0; i<edgelen; i++) {
        for(int j=0; j<edgelen; j++) {
            int spotsum = (v1.lattice[j + i*edgelen] + v2.lattice[j + i*edgelen]);
            count += spotsum;
            cout << spotsum << " ";
        }
        cout << endl;
    }
    
    cout << "count sum = " << count << endl;
    
}