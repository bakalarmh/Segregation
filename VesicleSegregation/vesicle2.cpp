//
//  vesicle.cpp
//  Segregation
//
//  Created by Matthew Bakalar on 2/1/13.
//  Copyright (c) 2013 Matthew Bakalar. All rights reserved.
//

#include "vesicle.h"
#include <cmath>
#include <iomanip>

using namespace std;

Vesicle::Vesicle(int isize, int contact_size): size(isize), lattice(isize), molecules(vector<Point2<int>>()) {
    gsl_rng_env_setup();
    ranT = gsl_rng_default;
    ranr = gsl_rng_alloc(ranT);
}

void Vesicle::SeedRandParticles(std::vector<int> type, std::vector<int> count) {
    int i=0;
    for(int t : type) {
        for(int j=0; j<count[i]; j++) {

            int rx = gsl_rng_uniform(ranr)*size;
            int ry = gsl_rng_uniform(ranr)*size;
            while(lattice[rx + ry*lattice.size] != 0) {
                rx = gsl_rng_uniform(ranr)*size;
                ry = gsl_rng_uniform(ranr)*size;
            }
            lattice[rx + ry*lattice.size] = t;
            molecules.push_back(Point2<int>(rx,ry));
        }
        i += 1;
    }
}

void Vesicle::GenerateTrialStep() {
    long count = molecules.size();
    float rn1 = gsl_rng_uniform(ranr);
    
    trial_idx = rn1*count;
    trial_source = molecules[trial_idx];
    trial_dest = trial_source;
    float rn2 = gsl_rng_uniform(ranr);
    int dir = floor(rn2*4);
    
    if(dir == 0) {
        int np = trial_source[0]-1;
        if(np < 0) {
            np += size;
        }
        trial_dest[0] = np;
    }
    else if(dir == 1) {
        int np = trial_source[0]+1;
        if(np >= size) {
            np -= size;
        }
        trial_dest[0] = np;
    }
    else if(dir == 2) {
        int np = trial_source[1]-1;
        if(np < 0) {
            np += size;
        }
        trial_dest[1] = np;
    }
    else if(dir == 3) {
        int np = trial_source[1]+1;
        if(np >= size) {
            np -= size;
        }
        trial_dest[1] = np;
    }
}

int Vesicle::TestTrialStepWithInterface(int radius, float U) {
    int dest_type = lattice[trial_dest[0] + trial_dest[1]*size];
    int source_type = lattice[trial_source[0] + trial_source[1]*size];
    if (dest_type == 0) {
        int center = floor(size/2);
        int dest_inside = (sqrt(pow(trial_dest[0]-center,2) + pow(trial_dest[1]-center,2))) < radius;
        int source_inside = (sqrt(pow(trial_source[0]-center,2) + pow(trial_source[1]-center,2))) < radius;
        if (source_type == 2) {
            if(source_inside == dest_inside) {
                return 1;
            }
            else if(source_inside == 0) {
                // Outside to inside
                return 1;
            }
            else {
                // Inside to outside
                float r1 = gsl_rng_uniform(ranr);
                float mc = exp(U);
                return mc > r1;
            }
        }
        else {
            return 1;
        }
    }
    else {
        return 0;
    }
}

int Vesicle::CountParticlesInside(int radius, int type) {
    int inside = 0;
    int center = floor(size/2);
    for(Point2<int> m : molecules) {
        if(lattice[m[0] + m[1]*size] == type) {
            inside += (sqrt(pow(m[0]-center,2) + pow(m[1]-center,2))) < radius;
        }
    }
    return inside;
}

void Vesicle::AcceptTrialMove() {
    lattice[trial_dest[0] + trial_dest[1]*size] = lattice[trial_source[0] + trial_source[1]*size];
    lattice[trial_source[0] + trial_source[1]*size] = 0;
    molecules[trial_idx][0] = trial_dest[0];
    molecules[trial_idx][1] = trial_dest[1];
}

int Vesicle::MoveRandomParticle() {
    long count = molecules.size();
    float rn1 = gsl_rng_uniform(ranr);
    
    Point2<int> frompt = molecules[floor(rn1*count)];
    Point2<int> topt = frompt;
    float rn2 = gsl_rng_uniform(ranr);
    int dir = floor(rn2*4);
    
    if(dir == 0) {
        int np = frompt[0]-1;
        if(np < 0) {
            np += size;
        }
        topt[0] = np;
    }
    else if(dir == 1) {
        int np = frompt[0]+1;
        if(np >= size) {
            np -= size;
        }
        topt[0] = np;
    }
    else if(dir == 2) {
        int np = frompt[1]-1;
        if(np < 0) {
            np += size;
        }
        topt[1] = np;
    }
    else if(dir == 3) {
        int np = frompt[1]+1;
        if(np >= size) {
            np -= size;
        }
        topt[1] = np;
    }
    
    if(lattice[topt[0] + topt[1]*size] == 0) {
        lattice[topt[0] + topt[1]*size] = lattice[frompt[0] + frompt[1]*size];
        lattice[frompt[0] + frompt[1]*size] = 0;
        molecules[floor(rn1*count)][0] = topt[0];
        molecules[floor(rn1*count)][1] = topt[1];
        return 1;
    }
    else {
        return 0;
    }
}

ostream& Vesicle::PrintXYZ(ostream& os, vector<string> atoms) {
    long count = molecules.size();
    os << "\t" << count << "\n";
    os << "Title Line \n";
    for(int i=0; i<count; i++) {
        float x = (float)molecules[i][0];
        float y = (float)molecules[i][1];
        int type = lattice[x + y*size];
        std::cout.setf(std::ios::fixed);
        std::cout.precision(2);
        os.setf(ios::fixed);
        os.precision(10);
        os << atoms[type-1] << setw(15) << x << setw(15) << y << setw(15) << (float)0.0 << "\n";
    }
    return os;
}
