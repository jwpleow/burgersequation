/**
* Code to solve General 2D Burger's Equation using an Explicit
* Forward Euler time-integration scheme.
*
*
*
*
*
*
*
*  Written by Joel Leow
*/
#include <chrono>
#include <string>
#include <iostream>

#include "Model.h"
#include "Burgers.h"

int main(int argc, char* argv[]) {





    Model m(argc, argv);
    m.PrintParameters();
    Burgers b(m);
    b.SetVelField(m);
    b.DisplayxVelField(m);

//    // Call code to initialise the problem here
//
//    typedef std::chrono::high_resolution_clock hrc;
//    typedef std::chrono::milliseconds ms;
//    hrc::time_point start = hrc::now();
//
//    // Call code to perform time integration here
//
//    hrc::time_point end = hrc::now();
//
//    // Calculate final energy and write output
//
    return 0;
}

