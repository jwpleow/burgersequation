// Created by jwl16 on 18/02/19.

#include <iostream>
#include <string>

#include "Model.h"

/// Constructor that takes in the required parameters from the argument
Model::Model(int argc, char *argv[]) {

    // Initialise and check MPI
    retval = MPI_Init(&argc, &argv);
    if(retval != MPI_SUCCESS){
        throw std::runtime_error("An error occurred initialising MPI");
    }

    retval_rank = MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    retval_size = MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (retval_rank == MPI_ERR_COMM || retval_size == MPI_ERR_COMM){
        std::cout << "Invalid Communicator!" << std::endl;
    }


    // Check ranks are properly allocated
    std::cout << "My rank is: " << world_rank << '\n';

    ParseParameters(argc, argv);
    ValidateParameters(argc);

    // Arrange process ranks in column major format
    // Find the localNx and Ny sizes, as well as the start location of each array



    // allocate the last columns of processes to have the extra+1 in localNx (LHS of if is column no., RHS is the columns that will not need to be assigned the extra+1 localNx)
    if ( world_rank / nPy > nPx - 1 - (Nx - 2) % nPx) {
        localNx = (Nx - 2) / nPx + 2 + 1;
    } else {
        localNx = (Nx - 2) / nPx + 2;
    }
    // allocate the last rows of the processes to have the extra+1 in localNy (LHS of if is row no., RHS is the rows that will not need to be assigned the extra+1 localNx)
    if ( world_rank % nPy > nPy - 1 - (Ny - 2) % nPy) {
        localNy = (Ny - 2) / nPy + 2 + 1;
    } else {
        localNy = (Ny - 2) / nPy + 2;
    }


//    } else { ///< for other processes
//
//        if (nPx >  nPx - 1 - (Nx - 2) % nPx) { ///< distribute the extra array lengths needed equally among last few processes
//            localNx = (Nx - 2) / world_size + 2 + 1;
//        } else {
//            localNx = (Nx - 2) / world_size + 2;
//        }







//    // Check if arrays can be subdivided equally, and assign them their Nx sizes and starting location in the array
//    if (Nx % 2 == 0) {
//        localNx = (Nx / 2 + 1);
//        localstart = world_rank * ((localNx - 2) * localNy);
//    } else {
//        if (world_rank == 0) {
//            localNx = (Nx + 1) / 2;
//        } else{
//            localNx = (Nx + 1) / 2 + 1;
//        }
//        localstart = world_rank * (((Nx - 1) / world_size - 1) * localNy);
//    }

}

// Destructor
Model::~Model() {

    // Terminate Parallel Execution environment
    MPI_Finalize();

}

void Model::ParseParameters(int argc, char *argv[]) {

    /// Parse argument char array into the relevant variables (sto also checks for invalid datatype args)

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
    nPx = std::stoi(argv[11]);
    nPy = std::stoi(argv[12]);
    dx = Lx / (Nx - 1);
    dy = Ly / (Ny - 1);
    dt = T / (Nt);
    x0 = -Lx / 2;
    y0 = -Ly / 2;
}

void Model::ValidateParameters(int argc) {
    // Check Px * Py = number of processes
    if (nPx * nPy != world_size) {
        throw std::invalid_argument("Number of processes do not match -np initialisation!");
    }
    // Check number of arguments are correct
    if (argc != 13) {
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
















