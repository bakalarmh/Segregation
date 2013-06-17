//
//  simconfig.cpp
//  VesicleSegregation
//
//  Created by Matthew Bakalar on 6/10/13.
//  Copyright (c) 2013 Matthew Bakalar. All rights reserved.
//

#include "simconfig.h"
#include <fstream>
#include <sstream>
#include <string>
#include "tests.h"
#include "vesicle.h"
#include "vesicledimer.h"
#include "vesiclemanager.h"

using namespace std;

SimConfig::SimConfig() {
    // Pass
}

void SimConfig::SimulationA() {
    string outpath = "/Users/mhbakalar/Documents/Fletchlab/Data/Segsim/dual_lattice_2D/testoutput/state.txt";
    ofstream outfile;
    outfile.open(outpath);
    
    int edgelen = 1000;
    int contact = 100;
    float rho = 0.01;
    int ntypes = 2;
    
    Vesicle v1(edgelen, contact, ntypes);
    Vesicle v2(edgelen, contact, ntypes);
    
    // Name molecules for vmd printing
    const char* typenm_arr[] = {"C", "H"};
    vector<string> atoms(typenm_arr, end(typenm_arr));
    
    // No one will crowd across the interface today!
    vector<tuple<int,int>> crowdPairs;
    
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
    bindingEnergy[1] = -5.0;
    bindingEnergy[2] = 0.0;
    
    VesicleDimer dimer(v1, v2, crowdPairs, typeBinds, contactEnergy, bindingEnergy);
    VesicleManager manager;
    int n_particles = rho*(edgelen*edgelen);
    int seeded;
    seeded = manager.SeedParticles(dimer, 0, 1, n_particles);
    seeded = manager.SeedParticles(dimer, 1, 1, n_particles);
    
    cout << "Seeded: " << seeded << endl;
    
    // Simulation
    outfile << "Particles = " << n_particles << endl;
    
    int nsweeps = 1000;
    int steps = n_particles*1000;
    int acceptcount = 0;
    for(int sweep=0; sweep<nsweeps; sweep++) {
        for(int mol=0; mol<steps; mol++) {
            int accepted = 0;
            while (accepted == 0) {
                accepted = manager.AttemptMove(dimer);
                acceptcount += accepted;
            }
        }
        // Save the system state for vesicle v0
        // Record only particle type 1 for now
        outfile << v1.numContact[1] << " " << v1.numBound[1] << endl;
    }
    outfile << endl;
    outfile.close();
}

void SimConfig::SimulationB() {
    string basepath = "/Users/mhbakalar/Documents/Fletchlab/Data/Segsim/dual_lattice_2D/scan_U/";
    int edgelen = 200;
    int contact = 100;
    float rho = 0.01;
    int ntypes = 1;
    
    int nruns = 8;
    float U_max = -8.0;
    
    for(int run=0; run < nruns; run++) {
        
        float U = ((float)run/nruns) * U_max;
        
        ostringstream os;
        os << basepath << "U=" << U << "_p=" << rho << "_edge=" << edgelen << "_contact=" << contact << "_003.txt";
        string outpath = os.str();
        
        ofstream outfile;
        outfile.open(outpath);
        
        Vesicle v1(edgelen, contact, ntypes);
        Vesicle v2(edgelen, contact, ntypes);
        
        // Name molecules for vmd printing
        const char* typenm_arr[] = {"C"};
        vector<string> atoms(typenm_arr, end(typenm_arr));
        
        // No one will crowd across the interface today!
        vector<tuple<int,int>> crowdPairs;
        
        // Binding types
        const bool typebinds_arr[] = {false, true};
        vector<bool> typeBinds(typebinds_arr, end(typebinds_arr));
        
        // Define energy for particles in special states
        vector<float> contactEnergy(2);
        contactEnergy[0] = 0.0;
        contactEnergy[1] = 0.0;
        vector<float> bindingEnergy(2);
        bindingEnergy[0] = 0.0;
        bindingEnergy[1] = U;
        
        VesicleDimer dimer(v1, v2, crowdPairs, typeBinds, contactEnergy, bindingEnergy);
        VesicleManager manager;
        int n_particles = rho*(edgelen*edgelen);
        int seeded;
        seeded = manager.SeedParticles(dimer, 0, 1, n_particles);
        seeded = manager.SeedParticles(dimer, 1, 1, n_particles);
        
        cout << "Seeded: " << seeded << endl;
        
        // Simulation
        outfile << "Particles Edgelen Contact U" << n_particles << endl;
        outfile << n_particles << " " << edgelen << " " << contact << " " << U << endl;
        outfile << endl;
        
        int nsweeps = 1000;
        int steps = n_particles*100;
        int acceptcount = 0;
        for(int sweep=0; sweep<nsweeps; sweep++) {
            for(int mol=0; mol<steps; mol++) {
                int accepted = 0;
                while (accepted == 0) {
                    accepted = manager.AttemptMove(dimer);
                    acceptcount += accepted;
                }
            }
            // Save the system state for vesicle v0
            // Record only particle type 1 for now
            outfile << v1.numContact[1] << " " << v1.numBound[1] << endl;
        }
        outfile << endl;
        outfile.close();
    }
} 

void SimConfig::SimulationC() {
    string basepath = "/Users/mhbakalar/Documents/Fletchlab/Data/Segsim/dual_lattice_2D/rhoscan/";
    int edgelen = 200;
    int contact = 100;
    float U = -6.0;
    float rhovals[] = {0.15, 0.2, 0.25};
    int ntypes = 1;
    
    int nruns = 9;
    
    for(int run=0; run < nruns; run++) {
        
        float rho = rhovals[run];
        
        ostringstream os;
        os << basepath << "U=" << U << "_p=" << rho << "_edge=" << edgelen << "_contact=" << contact << "_001.txt";
        string outpath = os.str();
        
        ofstream outfile;
        outfile.open(outpath);
        
        Vesicle v1(edgelen, contact, ntypes);
        Vesicle v2(edgelen, contact, ntypes);
        
        // Name molecules for vmd printing
        const char* typenm_arr[] = {"C"};
        vector<string> atoms(typenm_arr, end(typenm_arr));
        
        // No one will crowd across the interface today!
        vector<tuple<int,int>> crowdPairs;
        
        // Binding types
        const bool typebinds_arr[] = {false, true};
        vector<bool> typeBinds(typebinds_arr, end(typebinds_arr));
        
        // Define energy for particles in special states
        vector<float> contactEnergy(2);
        contactEnergy[0] = 0.0;
        contactEnergy[1] = 0.0;
        vector<float> bindingEnergy(2);
        bindingEnergy[0] = 0.0;
        bindingEnergy[1] = U;
        
        VesicleDimer dimer(v1, v2, crowdPairs, typeBinds, contactEnergy, bindingEnergy);
        VesicleManager manager;
        int n_particles = rho*(edgelen*edgelen);
        int seeded;
        seeded = manager.SeedParticlesCarefully(dimer, 0, 1, n_particles);
        seeded = manager.SeedParticlesCarefully(dimer, 1, 1, n_particles);
        
        cout << "Seeded: " << seeded << endl;
        
        // Simulation
        outfile << "Particles Edgelen Contact U" << n_particles << endl;
        outfile << n_particles << " " << edgelen << " " << contact << " " << U << endl;
        outfile << endl;
        
        int nsweeps = 40000;
        int steps = n_particles*1000;
        int acceptcount = 0;
        for(int sweep=0; sweep<nsweeps; sweep++) {
            for(int mol=0; mol<steps; mol++) {
                int accepted = 0;
                accepted = manager.AttemptMove(dimer);
                acceptcount += accepted;
            }
            // Save the system state for vesicle v0
            // Record only particle type 1 for now
            outfile << v1.numContact[1] << " " << v1.numBound[1] << endl;
        }
        outfile << endl;
        outfile.close();
    }
}