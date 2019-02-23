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

    // Initialising the problem
    Model m(argc, argv); ///< Initialise Model with input command line argument as parameters
    m.PrintParameters(); ///< Print Parameters to make sure they are correct
    Burgers b(m); ///< Initialise Burgers Class which will contain the data with our parameters from Model
    b.SetVelField(m); ///< Calculate the initial velocity field from initial conditions
//    b.DisplayuVelField(m);


    // Check the time taken for calculation
    std::cout << "Time integrating velocity field...\n";
    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();

    // Time integrate the velocity field
    b.TimeIntegrateVelField(m);

    hrc::time_point end = hrc::now();
    std::cout << "Time taken: " << std::chrono::duration_cast<ms>(end - start).count() << "ms\n";

    // Calculate final energy and write output

    std::cout << "Energy of velocity field: " << std::setprecision(5) << b.EnergyOfVelField(m) << std::endl;

    //  b.DisplayuVelField(m);

    b.PrintVelFields(m); ///< Print the final velocity field to VelocityFields.txt

    return 0;
}

