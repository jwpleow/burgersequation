// Created by jwl16 on 18/02/19.

#include <iostream>
#include <string>

#include "Model.h"

/// Constructor that takes in the required parameters from the argument
Model::Model(int argc, char *argv[]) {

    // Call functions to parse and validate the parameters provided through the command line
    ParseParameters(argc, argv);
    ValidateParameters(argc);

}

// Destructor
Model::~Model() {


}

void Model::ParseParameters(int argc, char *argv[]) {

    /// Parse argument char array into the relevant variables

    Lx = std::stod(argv[1]);
    Ly = std::stod(argv[2]);
    T = std::stod(argv[3]);
    Nx = std::stoi(argv[4]);
    Ny = std::stoi(argv[5]);
    Nt = std::stoi(argv[6]);
    ax = std::stod(argv[7]);
    ay = std::stod(argv[8]);
    b = std::stod(argv[9]);
    c = std::stod(argv[10]);
    dx = Lx / (Nx - 1);
    dy = Ly / (Ny - 1);
    dt = T / (Nt);
    x0 = -Lx / 2;
    y0 = -Ly / 2;
}

void Model::ValidateParameters(int argc) {
    // Check number of arguments are correct
    if (argc != 11) {
        throw std::invalid_argument("Incorrect number of arguments!");
    }
    // Check for negative final time input
    if (T <= 0) {
        throw std::invalid_argument("Final Time must be greater than 0!");
    }
    // Timestep & space discretisation validation
    if (Nx <= 0 || Ny <= 0 || Nt <= 0) {
        throw std::invalid_argument("Number of timesteps/grid points must be greater than 0!");
    }
    // Check domain given is OK
    if (Lx <= 0 || Ly <= 0) {
        throw std::invalid_argument("Domain given is invalid!");
    }

}



// Member function to print all parameters
void Model::PrintParameters() {
    std::cout << "x0: " << x0 << '\n';
    std::cout << "y0: " << y0 << '\n';
    std::cout << "Lx: " << Lx << '\n';
    std::cout << "Ly: " << Ly << '\n';
    std::cout << "T:  " << T << '\n';
    std::cout << "Nx: " << Nx << '\n';
    std::cout << "Ny: " << Ny << '\n';
    std::cout << "Nt: " << Nt << '\n';
    std::cout << "dx: " << dx << '\n';
    std::cout << "dy: " << dy << '\n';
    std::cout << "dt: " << dt << '\n';
    std::cout << "ax: " << ax << '\n';
    std::cout << "ay: " << ay << '\n';
    std::cout << "b:  " << b << '\n';
    std::cout << "c:  " << c << '\n';
}
















