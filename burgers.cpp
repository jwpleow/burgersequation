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
    if(m.world_rank == 0) m.PrintParameters(); ///< Print Parameters to make sure they are correct
    Burgers b(m); ///< Initialise Burgers Class which will contain the data with our parameters from Model
    b.SetVelField(m); ///< Calculate the initial velocity field from initial conditions




    if (m.world_rank == 0) std::cout << "Time integrating velocity field...\n";
    // Check the time taken for calculation
    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();

//     Time integrate the velocity field
    b.TimeIntegrateVelField(m);

    // Print time taken to time integrate the field
    hrc::time_point end = hrc::now();
    std::cout << "Time taken: " << std::chrono::duration_cast<ms>(end - start).count() << "ms\n";



//    if (m.world_rank == 0) {
//        std::cout << " combined vel field: \n";
//        b.DisplayCombinedufield(m);
//
//    }
//
    // velocity field printers
//    if (m.world_rank == 0) {
//        std::chrono::seconds dura(1);
//        std::this_thread::sleep_for(dura);
//        std::cout << m.world_rank << " vel field: \n";
//        b.DisplayuVelField(m);
//
//    }
//    if (m.world_rank == 1) {
//        std::chrono::seconds dura(3);
//        std::this_thread::sleep_for(dura);
//        std::cout << m.world_rank << " vel field: \n";
//        b.DisplayuVelField(m);
//    }
//    if (m.world_rank == 2) {
//        std::chrono::seconds dura(4);
//        std::this_thread::sleep_for(dura);
//        std::cout << m.world_rank << " vel field: \n";
//        b.DisplayuVelField(m);
//    }
//
//    if (m.world_rank == 3) {
//        std::chrono::seconds dura(5);
//        std::this_thread::sleep_for(dura);
//        std::cout << m.world_rank << " vel field: \n";
//        b.DisplayuVelField(m);
//    }
//    if (m.world_rank == 4) {
//        std::chrono::seconds dura(6);
//        std::this_thread::sleep_for(dura);
//        std::cout << m.world_rank << " vel field: \n";
//        b.DisplayuVelField(m);
//    }
//    if (m.world_rank == 5) {
//        std::chrono::seconds dura(7);
//        std::this_thread::sleep_for(dura);
//        std::cout << m.world_rank << " vel field: \n";
//        b.DisplayuVelField(m);
//    }
//    if (m.world_rank == 6) {
//        std::chrono::seconds dura(8);
//        std::this_thread::sleep_for(dura);
//        std::cout << m.world_rank << " vel field: \n";
//        b.DisplayuVelField(m);
//    }
//
//    if (m.world_rank == 7) {
//        std::chrono::seconds dura(9);
//        std::this_thread::sleep_for(dura);
//        std::cout << m.world_rank << " vel field: \n";
//        b.DisplayuVelField(m);
//    }
//    if (m.world_rank == 8) {
//        std::chrono::seconds dura(10);
//        std::this_thread::sleep_for(dura);
//        std::cout << m.world_rank << " vel field: \n";
//        b.DisplayuVelField(m);
//    }
//    if (m.world_rank == 9) {
//        std::chrono::seconds dura(11);
//        std::this_thread::sleep_for(dura);
//        std::cout << m.world_rank << " vel field: \n";
//        b.DisplayuVelField(m);
//    }
//    if (m.world_rank == 10) {
//        std::chrono::seconds dura(12);
//        std::this_thread::sleep_for(dura);
//        std::cout << m.world_rank << " vel field: \n";
//        b.DisplayuVelField(m);
//    }
//    if (m.world_rank == 11) {
//        std::chrono::seconds dura(13);
//        std::this_thread::sleep_for(dura);
//        std::cout << m.world_rank << " vel field: \n";
//        b.DisplayuVelField(m);
//    }


//    // Calculate final energy
//    if (m.world_rank == 0)  std::cout << "Energy of velocity field: " << std::setprecision(5) << b.EnergyOfVelField(m) << std::endl;

    // Print the final velocity field to VelocityFields.txt
//    b.FilePrintVelFields(m);


    return 0;
}

