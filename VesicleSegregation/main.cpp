//
//  main.cpp
//  VesicleSegregation
//
//  Created by Matthew Bakalar on 5/31/13.
//  Copyright (c) 2013 Matthew Bakalar. All rights reserved.
//

#include <iostream>
#include "tests.h"
#include "simconfig.h"

int main(int argc, const char * argv[])
{
    Tests tests;
    // tests.VesicleTest();
    // tests.ContactTest();
    // tests.DimerTest();
    // tests.DimerContactTest();
    // tests.BrownianTest();
    // tests.CarefulSeedTest();

    
    SimConfig sim;
    // sim.SimulationA();
    // sim.SimulationB();
    sim.SimulationC();
    
    return 0;
}

