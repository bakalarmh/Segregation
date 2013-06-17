//
//  Level0.cpp
//  Segregation
//
//  Created by Matthew Bakalar on 1/31/13.
//  Copyright (c) 2013 Matthew Bakalar. All rights reserved.
//

#include "level0.h"
#include "lattice2.h"
#include "Vesicle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include "point2.h"

using namespace std;

void Level0::AdhesionPartition() {
    float rhobase = 0.01;
    float rhoadd = 0.01;
    float rhob = 0.05;
    float Ubase = -2.0;
    float Uadd = -0.1;
    
    std::string basepath = "/Users/mhbakalar/Dropbox/Projects/Chem220B/Segregation/data/adhesion_rho_scan_U=2.0/";
    
    //for (int Umul=1; Umul<20; Umul++) {
    int Umul = 0.0;
        for(int rhomul=0; rhomul<10; rhomul++) {
            
            cout << "Hello world!" << endl;
            
            float rho = rhobase + rhoadd*rhomul;
            float U_energy = Ubase + Uadd*Umul;
            
            int size = 100;
            int mols = floor(((size*size)*rho));
            int mols_stable = floor(((size*size)*rhob));
            int contact_size = 25;
            
            ostringstream os1;
            ostringstream os2;
            os1 << basepath << "trajectory_U=" << U_energy << "_p=" << rho << "_sz=" << size << "_r=" << contact_size << ".xyz";
            string traj_path = os1.str();
            os2 << basepath << "counts_U=" << U_energy << "_p=" << rho << "_sz=" << size << "_r=" << contact_size << ".txt";
            string counts_path = os2.str();
            
            ofstream trajfile;
            ofstream countfile;
            
            // trajfile.open(traj_path);
            countfile.open(counts_path);
            
            Vesicle vesicle(size, contact_size);
            
            // Initialize vesicle with mols particles of types 1, mols particles of type 2
            static const int types_arr[] = {1,2};
            const char* typenm_arr[] = {"C", "H"};
            vector<int> types (types_arr, types_arr + sizeof(types_arr) / sizeof(types_arr[0]) );
            vector<int> counts(0);
            
            cout << mols << " a" << endl;
            cout << mols_stable << " b" << endl;
            counts.push_back(mols_stable);
            counts.push_back(mols);
            vector<string> atoms(typenm_arr, end(typenm_arr));
            
            cout << counts[0] << ", " << counts[1] << endl;
            
            vesicle.SeedRandParticles(types, counts);
            
            for(int sweep=0; sweep<1000; sweep++) {
                for(int run=0; run<mols*1000; run++) {
                    vesicle.GenerateTrialStep();
                    if(vesicle.TestTrialStepWithInterface(contact_size, U_energy) == 1) {
                        vesicle.AcceptTrialMove();
                    }
                }
                
                cout << "sweep: " << sweep << endl;
                if(sweep > 100) {
                    int in1 = vesicle.CountParticlesInside(contact_size, 1);
                    int in2 = vesicle.CountParticlesInside(contact_size, 2);
                    countfile << sweep << " " << in1 << " " << in2 << endl;
                    cout << in2 << endl;
                    // vesicle.PrintXYZ(trajfile, atoms);
                }
            }
            
            // trajfile << endl;
            countfile.close();
            // trajfile.close();
        }
    //}
}