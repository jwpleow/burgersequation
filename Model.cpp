// Created by jwl16 on 18/02/19.

#include <iostream>
#include <string>

#include "Model.h"

/// Constructor that takes in the required parameters from the argument
Model::Model(int argc, char *argv[]) {
    /// Check for correct number of arguments
    if (argc != 11) {
        throw std::invalid_argument("Incorrect number of input arguments!");
    }

    /// Parse argument char array into the relevant variables
    Lx = std::stod(argv[1]);
    Ly = std::stod(argv[2]);

    // Check for negative input
    if (std::stod(argv[3]) <= 0) {
        throw std::invalid_argument("Final Time must be greater than 0!");
    }
    T = std::stod(argv[3]);

    // Timestep & space discretisation validation and parse - check for negative input
    if (std::stoi(argv[4]) <= 0 || std::stoi(argv[5]) <= 0 || std::stoi(argv[6]) <= 0) {
        throw std::invalid_argument("Number of timesteps/grid points must be greater than 0!");
    }
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

// Destructor
Model::~Model() {}

// Member function that can reparse parameters
void Model::ParseParameters(int argc, char *argv[]) {
    /// Check for correct number of arguments
    if (argc != 11) {
        throw std::invalid_argument("Incorrect number of input arguments!");
    }
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
    if (Nx <= 0 || Ny <= 0 || Nt <= 0) {
        throw std::invalid_argument("Number of timesteps/grid points must be greater than 0!");
    }
    dx = Lx / (Nx - 1);
    dy = Ly / (Ny - 1);
    dt = T / (Nt);
    x0 = -Lx / 2;
    y0 = -Ly / 2;
}

void Model::ValidateParameters(){
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
















