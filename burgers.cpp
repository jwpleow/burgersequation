/**
* C++ Code (in serial) to solve General 2D Burger's Equation using
* an Explicit Forward Euler time-integration scheme.
*
*
* Use the Makefile included to run the program (e.g. make burg)
* Parameters of the problem can be edited in the Makefile
*
* For different initial conditions, edit Burgers::SetVelField in Burgers.cpp
* Arrays are stored in column-major format
*  Written by Joel Leow - 23/02/2019
*/
#include <chrono>
#include <iostream>

#include "Model.h"
#include "Burgers.h"

#include <thread>

int main(int argc, char* argv[]) {

    // Initialising the problem
    Model m(argc, argv); ///< Initialise Model with input command line argument as parameters
    if (m.world_rank == 0) m.PrintParameters(); ///< Print Parameters to make sure they are correct
    Burgers b(m); ///< Initialise the Burgers Class which will contain the data with our parameters from Model
    b.SetVelField(m); ///< Calculate the initial velocity field


    // Check the time taken for the time integration
    if (m.world_rank == 0) std::cout << "Time integrating velocity field...\n";
    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();

    // Time integrate the velocity field
    b.TimeIntegrateVelField(m);

    // Print time taken to time integrate the field
    hrc::time_point end = hrc::now();
    std::cout << "Rank: " << m.world_rank << ", Time taken: " << std::chrono::duration_cast<ms>(end - start).count()
              << "ms\n";

    // Display the Total Energy
    if (m.world_rank == 0)  std::cout << "Energy of velocity field: " << std::setprecision(5) << b.EnergyOfVelField(m) << std::endl;

    // Print the final velocity field to VelocityFields.txt
    b.FilePrintVelFields(m);


    return 0;
}