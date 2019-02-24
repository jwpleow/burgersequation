/**
* C++ Code (in serial) to solve General 2D Burger's Equation using
* an Explicit Forward Euler time-integration scheme.
*
*
* Use the Makefile included to run the program (e.g. make burg)
* Parameters of the problem can be edited in the Makefile
*
* For different initial conditions, edit Burgers::SetVelField in Burgers.cpp
*
*  Written by Joel Leow - 23/02/2019
*/
#include <chrono>
#include <iostream>

#include "Model.h"
#include "Burgers.h"

int main(int argc, char* argv[]) {


    // Initialising the problem
    Model m(argc, argv); ///< Initialise Model with input command line argument as parameters
    m.PrintParameters(); ///< Print Parameters to make sure they are correct
    Burgers b(m); ///< Initialise Burgers Class which will contain the data with our parameters from Model
    b.SetVelField(m); ///< Calculate the initial velocity field from initial conditions


    // Check the time taken for calculation
    std::cout << "Time integrating velocity field...\n";
    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();

    // Time integrate the velocity field
    b.TimeIntegrateVelField(m);

    // Print time taken to time integrate the field
    hrc::time_point end = hrc::now();
    std::cout << "Time taken: " << std::chrono::duration_cast<ms>(end - start).count() << "ms\n";

    // Calculate final energy
    std::cout << "Energy of velocity field: " << std::setprecision(5) << b.EnergyOfVelField(m) << std::endl;

    // Print the final velocity field to VelocityFields.txt
    b.PrintVelFields(m);


    return 0;
}

