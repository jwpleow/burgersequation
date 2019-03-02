// Created by jwl16 on 18/02/19.

#include "Burgers.h"

// Constructor for Burgers that accepts an argument of class Model
Burgers::Burgers(Model &A) {
    // Use dynamic memory for data arrays - initialise values to 0 using ()
    udata_ = new double[A.localNx * A.localNy]();
    vdata_ = new double[A.localNx * A.localNy]();
    // Initialise extra arrays for switching between during the time integration
    udata_2 = new double[A.localNx * A.localNy]();
    vdata_2 = new double[A.localNx * A.localNy]();
}

// Destructor
Burgers::~Burgers() {
    delete[] udata_;
    delete[] vdata_;
//    if (A.world_rank = 0) {
//        delete[] ucombineddata_;
//        delete[] vcombineddata_;
//    }
};

// Member function to set the initial velocity field and distribute them to the processes
// Calculates the initial velocity field from the input parameters for the initial condition:
// At t=0; r=sqrt(x^2+y^2); u & v = 2 * (1 - r) ^4 * (4 * r + 1) for r <= 1 and = 0 for r > 1
void Burgers::SetVelField(Model &A) {

    double r; ///< Initialise r value

    // Create the array to store combined velocity data
    ucombineddata_ = new double[A.Nx * A.Ny]();
    vcombineddata_ = new double[A.Nx * A.Ny]();

    // Calculate the initial velocity field in rank 0 process
    if (A.world_rank == 0) {
        for (int i = 0; i < A.Nx; ++i) { ///< i is the column index
            for (int j = 0; j < A.Ny; ++j) { ///< j is the row index
                r = sqrt((A.x0 + i * A.dx) * (A.x0 + i * A.dx) + (A.y0 + j * A.dy) * (A.y0 + j * A.dy));
                if (r <= 1.0) {
                    ucombineddata_[A.Ny * i + j] = 2.0 * pow(1.0 - r, 4) * (4.0 * r + 1);
                    vcombineddata_[A.Ny * i + j] = 2.0 * pow(1.0 - r, 4) * (4.0 * r + 1);
                }
            }
        }
    }

    // broadcast the initial velocity data
    MPI_Bcast(ucombineddata_, A.Nx * A.Ny , MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(vcombineddata_, A.Nx * A.Ny , MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Split the initialised velocity data
    for (int i = 0; i < A.localNx; ++i) {
        for (int j = 0; j < A.localNy; ++j) {
            udata_[i * A.localNy + j] = ucombineddata_[A.localstart + i * A.Ny + j];
            vdata_[i * A.localNy + j] = vcombineddata_[A.localstart + i * A.Ny + j];
        }
    }

    // Free the unnecessary dynamic arrays on other processes
    if (A.world_rank != 0) {
        delete[] ucombineddata_;
        delete[] vcombineddata_;
    }

}

void Burgers::DisplayuVelField(Model &A) {

    std::cout << "Current u Velocity Field: \n";

    for (int i = 0; i < A.Ny; ++i) { ///< i is the row index
        for (int j = 0; j < A.Nx; ++j) { ///< j is the column index
            std::cout << std::fixed << std::setprecision(3) << ucombineddata_[A.Ny * j + i] << ' ';
        }
        std::cout << '\n';
    }
}

void Burgers::DisplayvVelField(Model &A) {

    std::cout << "Current v Velocity Field: \n";

    for (int i = 0; i < A.Ny; ++i) { ///< i is the row index
        for (int j = 0; j < A.Nx; ++j) { ///< j is the column index
            std::cout << std::fixed << std::setprecision(3) << vcombineddata_[A.Ny * j + i] << ' ';
        }
        std::cout << '\n';
    }
}

// Member function that integrates the velocity fields for the parameters specified in Model &A and with initial conditions from SetVelField
void Burgers::TimeIntegrateVelField(Model &A) {

    // Initialise arrays to transfer the row edges - can make this faster by bringing this out of timeintegrate
    auto* utoparraysend = new double[A.localNx]();
    auto* utoparrayreceive = new double[A.localNx]();
    auto* ubtmarraysend = new double[A.localNx]();
    auto* ubtmarrayreceive = new double[A.localNx]();
    auto* vtoparraysend = new double[A.localNx]();
    auto* vtoparrayreceive = new double[A.localNx]();
    auto* vbtmarraysend = new double[A.localNx]();
    auto* vbtmarrayreceive = new double[A.localNx]();

    // Precalculate constants for use in the multiplication
    double cdt_dxsq, cdt_dysq, axdt_dx, aydt_dy, bdt_dx, bdt_dy, const1, const2, const3;

    cdt_dxsq = A.c * A.dt / (A.dx * A.dx);
    cdt_dysq = A.c * A.dt / (A.dy * A.dy);
    axdt_dx = A.ax * A.dt / A.dx;
    aydt_dy = A.ay * A.dt / A.dy;
    bdt_dx = A.b * A.dt / A.dx;
    bdt_dy = A.b * A.dt / A.dy;
    const1 = 1.0 - 2.0 * cdt_dxsq - 2.0 * cdt_dysq - axdt_dx - aydt_dy;
    const2 = (cdt_dxsq + axdt_dx);
    const3 = (cdt_dysq + aydt_dy);

    // * * * * * * * * * * * * * Time Integration Loop * * * * * * * * * * * * * * * //
    for (int t = 0; t < A.Nt; ++t) { ///< iterate over the timesteps
        // Switch to alternate between using data in udata_/vdata_ and udata_2/vdata_2
        switch (t % 2) {
            case 0 :
                for (int i = 1; i < A.localNx - 1; ++i) { ///< i is the column tracker
                    for (int j = 1; j < A.localNy - 1; ++j) { ///< j is the row tracker
                        udata_2[A.localNy * i + j] =
                                (const3 + bdt_dx * vdata_[A.localNy * i + j]) * udata_[A.localNy * i + (j - 1)] +
                                (const1 - bdt_dy * vdata_[A.localNy * i + j] +
                                 bdt_dx * (udata_[A.localNy * (i - 1) + j] - udata_[A.localNy * i + j])) *
                                udata_[A.localNy * i + j] + cdt_dysq * udata_[A.localNy * i + (j + 1)] +
                                const2 * udata_[A.localNy * (i - 1) + j] + cdt_dxsq * udata_[A.localNy * (i + 1) + j];

                        vdata_2[A.localNy * i + j] =
                                const3 * vdata_[A.localNy * i + (j - 1)] + (const1 - bdt_dx * udata_[A.localNy * i + j]
                                                                            + bdt_dy *
                                                                              (vdata_[A.localNy * i + (j - 1)] -
                                                                               vdata_[A.localNy * i + j])) *
                                                                           vdata_[A.localNy * i + j] +
                                cdt_dysq * vdata_[A.localNy * i + (j + 1)] +
                                (const2 + bdt_dx * udata_[A.localNy * i + j]) * vdata_[A.localNy * (i - 1) + j] +
                                cdt_dxsq * vdata_[A.localNy * (i + 1) + j];

                    }
                }

                // If not part of rightmost column - Send and receive data from the right
                if (A.world_rank / A.nPy != A.nPx - 1){
                    MPI_Send(&udata_2[(A.localNx - 2) * A.localNy], A.localNy, MPI_DOUBLE, A.world_rank + A.nPy, 4 * A.world_rank, MPI_COMM_WORLD); ///< Send rightside array to the rightside process
                    MPI_Recv(&udata_2[(A.localNx - 1) * A.localNy], A.localNy, MPI_DOUBLE, A.world_rank + A.nPy, 4 * A.world_rank + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); ///< Receive rightside array from the right
                    MPI_Send(&vdata_2[(A.localNx - 2) * A.localNy], A.localNy, MPI_DOUBLE, A.world_rank + A.nPy, 4 * A.world_rank + 2, MPI_COMM_WORLD); ///< Send rightside array to the rightside process
                    MPI_Recv(&vdata_2[(A.localNx - 1) * A.localNy], A.localNy, MPI_DOUBLE, A.world_rank + A.nPy, 4 * A.world_rank + 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE); ///< Receive rightside array from the right
                }
                // If not part of leftmost column - Send and receive data from the left
                if (A.world_rank / A.nPy != 0){
                    MPI_Recv(&udata_2[0], A.localNy, MPI_DOUBLE, A.world_rank - A.nPy, 4 * A.world_rank - (4 * A.nPy), MPI_COMM_WORLD, MPI_STATUS_IGNORE); ///< Receive leftside array from the leftside process
                    MPI_Send(&udata_2[A.localNy], A.localNy, MPI_DOUBLE, A.world_rank - A.nPy, 4 * A.world_rank - (4 * A.nPy) + 1, MPI_COMM_WORLD); ///< Send leftside array to the leftside process
                    MPI_Recv(&vdata_2[0], A.localNy, MPI_DOUBLE, A.world_rank - A.nPy, 4 * A.world_rank - (4 * A.nPy) + 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE); ///< Receive leftside array from the leftside process
                    MPI_Send(&vdata_2[A.localNy], A.localNy, MPI_DOUBLE, A.world_rank - A.nPy, 4 * A.world_rank - (4 * A.nPy) + 3, MPI_COMM_WORLD); ///< Send leftside array to the leftside process

                }
                 //If not part of topmost row - Send and receive data from the top
                if (A.world_rank % A.nPy != 0) {
                    // Prepare an array of values from the 2nd row from the top
                    for (int i = 0; i < A.localNx; ++i) {
                        utoparraysend[i] = udata_2[i * A.localNy + 1];
                        vtoparraysend[i] = vdata_2[i * A.localNy + 1];
                    }
                    MPI_Send(utoparraysend, A.localNx, MPI_DOUBLE, A.world_rank - 1, 4 * (A.world_rank - A.world_rank / A.nPy) + (A.nPx - 1) * 4 * A.nPy - 4, MPI_COMM_WORLD); ///< Send top array to the process above
                    MPI_Recv(utoparrayreceive, A.localNx, MPI_DOUBLE, A.world_rank - 1, 4 * (A.world_rank - A.world_rank / A.nPy) + (A.nPx - 1) * 4 * A.nPy - 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE); ///< Receive the top array from the process above
                    MPI_Send(vtoparraysend, A.localNx, MPI_DOUBLE, A.world_rank - 1, 4 * (A.world_rank - A.world_rank / A.nPy) + (A.nPx - 1) * 4 * A.nPy - 2, MPI_COMM_WORLD); ///< Send top array to the process above
                    MPI_Recv(vtoparrayreceive, A.localNx, MPI_DOUBLE, A.world_rank - 1, 4 * (A.world_rank - A.world_rank / A.nPy) + (A.nPx - 1) * 4 * A.nPy - 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); ///< Receive the top array from the process above
                    // After receiving, parse the data into the top row
                    for (int i = 0; i < A.localNx; ++i) {
                        udata_2[i * A.localNy] = utoparrayreceive[i];
                        vdata_2[i * A.localNy] = vtoparrayreceive[i];
                    }
                }
                // If not part of bottom row - Send and receive data from the bottom
                if (A.world_rank % A.nPy != A.nPy - 1) {
                    // Prepare an array of values from the 2nd row from the bottom
                    for (int i = 0; i < A.localNx; ++i) {
                        ubtmarraysend[i] = udata_2[i * A.localNy + (A.localNy - 2)];
                        vbtmarraysend[i] = vdata_2[i * A.localNy + (A.localNy - 2)];
                    }
                    MPI_Recv(ubtmarrayreceive, A.localNx, MPI_DOUBLE, A.world_rank + 1, 4 * (A.world_rank - A.world_rank / A.nPy) + (A.nPx - 1) * 4 * A.nPy, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(ubtmarraysend, A.localNx, MPI_DOUBLE, A.world_rank + 1, 4 * (A.world_rank - A.world_rank / A.nPy) + (A.nPx - 1) * 4 * A.nPy + 1, MPI_COMM_WORLD);
                    MPI_Recv(vbtmarrayreceive, A.localNx, MPI_DOUBLE, A.world_rank + 1, 4 * (A.world_rank - A.world_rank / A.nPy) + (A.nPx - 1) * 4 * A.nPy + 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(vbtmarraysend, A.localNx, MPI_DOUBLE, A.world_rank + 1, 4 * (A.world_rank - A.world_rank / A.nPy) + (A.nPx - 1) * 4 * A.nPy + 3, MPI_COMM_WORLD);
                    // After receiving, parse the data into the bottom row
                    for (int i = 0; i < A.localNx; ++i) {
                        udata_2[i * A.localNy + (A.localNy - 1)] = ubtmarrayreceive[i];
                        vdata_2[i * A.localNy + (A.localNy - 1)] = vbtmarrayreceive[i];
                    }
              }
                break;
            case 1 :
                for (int i = 1; i < A.localNx - 1; ++i) { ///< i is the column tracker
                    for (int j = 1; j < A.localNy - 1; ++j) { ///< j is the row tracker
                        udata_[A.localNy * i + j] =
                                (const3 + bdt_dx * vdata_2[A.localNy * i + j]) * udata_2[A.localNy * i + (j - 1)] +
                                (const1 - bdt_dy * vdata_2[A.localNy * i + j]
                                 + bdt_dx * (udata_2[A.localNy * (i - 1) + j] - udata_2[A.localNy * i + j])) *
                                udata_2[A.localNy * i + j] +
                                cdt_dysq * udata_2[A.localNy * i + (j + 1)] +
                                const2 * udata_2[A.localNy * (i - 1) + j] +
                                cdt_dxsq * udata_2[A.localNy * (i + 1) + j];

                        vdata_[A.localNy * i + j] = const3 * vdata_2[A.localNy * i + (j - 1)] +
                                                    (const1 - bdt_dx * udata_2[A.localNy * i + j]
                                                     + bdt_dy * (vdata_2[A.localNy * i + (j - 1)] -
                                                                 vdata_2[A.localNy * i + j])) *
                                                    vdata_2[A.localNy * i + j] +
                                                    cdt_dysq * vdata_2[A.localNy * i + (j + 1)] +
                                                    (const2 + bdt_dx * udata_2[A.localNy * i + j]) *
                                                    vdata_2[A.localNy * (i - 1) + j] +
                                                    cdt_dxsq * vdata_2[A.localNy * (i + 1) + j];


                    }
                }

                // If not part of rightmost column - Send and receive data from the right
                if (A.world_rank / A.nPy != A.nPx - 1){
                    MPI_Send(&udata_[(A.localNx - 2) * A.localNy], A.localNy, MPI_DOUBLE, A.world_rank + A.nPy, 4 * A.world_rank, MPI_COMM_WORLD); ///< Send rightside array to the rightside process
                    MPI_Recv(&udata_[(A.localNx - 1) * A.localNy], A.localNy, MPI_DOUBLE, A.world_rank + A.nPy, 4 * A.world_rank + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); ///< Receive rightside array from the right
                    MPI_Send(&vdata_[(A.localNx - 2) * A.localNy], A.localNy, MPI_DOUBLE, A.world_rank + A.nPy, 4 * A.world_rank + 2, MPI_COMM_WORLD); ///< Send rightside array to the rightside process
                    MPI_Recv(&vdata_[(A.localNx - 1) * A.localNy], A.localNy, MPI_DOUBLE, A.world_rank + A.nPy, 4 * A.world_rank + 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE); ///< Receive rightside array from the right
                }
                // If not part of leftmost column - Send and receive data from the left
                if (A.world_rank / A.nPy != 0){
                    MPI_Recv(&udata_[0], A.localNy, MPI_DOUBLE, A.world_rank - A.nPy, 4 * A.world_rank - (4 * A.nPy), MPI_COMM_WORLD, MPI_STATUS_IGNORE); ///< Receive leftside array from the leftside process
                    MPI_Send(&udata_[A.localNy], A.localNy, MPI_DOUBLE, A.world_rank - A.nPy, 4 * A.world_rank - (4 * A.nPy) + 1, MPI_COMM_WORLD); ///< Send leftside array to the leftside process
                    MPI_Recv(&vdata_[0], A.localNy, MPI_DOUBLE, A.world_rank - A.nPy, 4 * A.world_rank - (4 * A.nPy) + 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE); ///< Receive leftside array from the leftside process
                    MPI_Send(&vdata_[A.localNy], A.localNy, MPI_DOUBLE, A.world_rank - A.nPy, 4 * A.world_rank - (4 * A.nPy) + 3, MPI_COMM_WORLD); ///< Send leftside array to the leftside process

                }
                //If not part of topmost row - Send and receive data from the top
                if (A.world_rank % A.nPy != 0) {
                    // Prepare an array of values from the 2nd row from the top
                    for (int i = 0; i < A.localNx; ++i) {
                        utoparraysend[i] = udata_[i * A.localNy + 1];
                        vtoparraysend[i] = vdata_[i * A.localNy + 1];
                    }
                    MPI_Send(utoparraysend, A.localNx, MPI_DOUBLE, A.world_rank - 1, 4 * (A.world_rank - A.world_rank / A.nPy) + (A.nPx - 1) * 4 * A.nPy - 4, MPI_COMM_WORLD); ///< Send top array to the process above
                    MPI_Recv(utoparrayreceive, A.localNx, MPI_DOUBLE, A.world_rank - 1, 4 * (A.world_rank - A.world_rank / A.nPy) + (A.nPx - 1) * 4 * A.nPy - 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE); ///< Receive the top array from the process above
                    MPI_Send(vtoparraysend, A.localNx, MPI_DOUBLE, A.world_rank - 1, 4 * (A.world_rank - A.world_rank / A.nPy) + (A.nPx - 1) * 4 * A.nPy - 2, MPI_COMM_WORLD); ///< Send top array to the process above
                    MPI_Recv(vtoparrayreceive, A.localNx, MPI_DOUBLE, A.world_rank - 1, 4 * (A.world_rank - A.world_rank / A.nPy) + (A.nPx - 1) * 4 * A.nPy - 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); ///< Receive the top array from the process above
                    // After receiving, parse the data into the top row
                    for (int i = 0; i < A.localNx; ++i) {
                        udata_[i * A.localNy] = utoparrayreceive[i];
                        vdata_[i * A.localNy] = vtoparrayreceive[i];
                    }
                }

                // If not part of bottom row - Send and receive data from the bottom
                if (A.world_rank % A.nPy != A.nPy - 1) {
                    // Prepare an array of values from the 2nd row from the bottom
                    for (int i = 0; i < A.localNx; ++i) {
                        ubtmarraysend[i] = udata_[i * A.localNy + (A.localNy - 2)];
                        vbtmarraysend[i] = vdata_[i * A.localNy + (A.localNy - 2)];
                    }
                    MPI_Recv(ubtmarrayreceive, A.localNx, MPI_DOUBLE, A.world_rank + 1, 4 * (A.world_rank - A.world_rank / A.nPy) + (A.nPx - 1) * 4 * A.nPy, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(ubtmarraysend, A.localNx, MPI_DOUBLE, A.world_rank + 1, 4 * (A.world_rank - A.world_rank / A.nPy) + (A.nPx - 1) * 4 * A.nPy + 1, MPI_COMM_WORLD);
                    MPI_Recv(vbtmarrayreceive, A.localNx, MPI_DOUBLE, A.world_rank + 1, 4 * (A.world_rank - A.world_rank / A.nPy) + (A.nPx - 1) * 4 * A.nPy + 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(vbtmarraysend, A.localNx, MPI_DOUBLE, A.world_rank + 1, 4 * (A.world_rank - A.world_rank / A.nPy) + (A.nPx - 1) * 4 * A.nPy + 3, MPI_COMM_WORLD);
                    // After receiving, parse the data into the bottom row
                    for (int i = 0; i < A.localNx; ++i) {
                        udata_[i * A.localNy + (A.localNy - 1)] = ubtmarrayreceive[i];
                        vdata_[i * A.localNy + (A.localNy - 1)] = vbtmarrayreceive[i];
                    }
                }
                break;
            default : ///< Show an error if the switching somehow breaks
                std::cout << "Error in switching during time integration";
                break;
        }

    } ///< timestep close bracket

    // Make sure that final data is always in udata_ and vdata_, and frees the memory for the unused case
    if (A.Nt % 2 == 1) {
        delete[] udata_;
        delete[] vdata_;
        udata_ = udata_2;
        vdata_ = vdata_2;
    } else {
        delete[] udata_2;
        delete[] vdata_2;
    }


    // * * * * * * * * * * * * * Assemble the combined data * * * * * * * * * * * * * * //

    // First, remove the edge columns/rows so that there is no overlapping data

    // If not part of rightmost column - remove rightmost column
    if (A.world_rank / A.nPy != A.nPx - 1){
        A.localNx -= 1;
    }
    // If not part of leftmost column - remove leftmost column
    if (A.world_rank / A.nPy != 0){
        A.localNx -= 1;
        for (int i = 0; i < A.localNx; ++i) {
            for (int j = 0; j < A.localNy; ++j) {
                udata_[j + A.localNy * i] = udata_[j + (i + 1) * A.localNy];
                vdata_[j + A.localNy * i] = vdata_[j + (i + 1) * A.localNy];
            }
        }
    }
    //If not part of topmost row - remove topmost row
    if (A.world_rank % A.nPy != 0) {
        for (int i = 0; i < A.localNx; ++i) {
            for (int j = 0; j < A.localNy - 1; ++j) {
                udata_[j + (A.localNy - 1) * i] = udata_[j + 1 + i * A.localNy];
                vdata_[j + (A.localNy - 1) * i] = vdata_[j + 1 + i * A.localNy];
            }
        }
        A.localNy -= 1;
    }
    // If not part of bottom row - remove bottommost row
    if (A.world_rank % A.nPy != A.nPy - 1) {
        A.localNy -= 1;
        for (int i = 0; i < A.localNx; ++i) {
            for (int j = 0; j < A.localNy; ++j) {
                udata_[j + A.localNy * i] = udata_[j + i * (A.localNy + 1)];
                vdata_[j + A.localNy * i] = vdata_[j + i * (A.localNy + 1)];
            }
        }
    }


    // Create an array to receive the data from all processes
    double* ureceivebuffer;
    double* vreceivebuffer;

    // Create the receive counts, displacements, and local size vectors required for MPI_Gatherv and rearranging
    int myLen = A.localNx * A.localNy;
    int *lengths = nullptr;
    int *displs = nullptr;
    int *localNys = nullptr;
    int *localNxs = nullptr;

    // Allocate the dynamic memory to hold the information required for assembly in the rank 0 process
    if(A.world_rank == 0) {
        lengths = new int[A.world_size]();
        displs = new int[A.world_size]();
        localNys = new int[A.world_size]();
        localNxs = new int[A.world_size]();
        ureceivebuffer = new double[A.Nx * A.Ny];
        vreceivebuffer = new double[A.Nx * A.Ny];
    }
    // Gather the lengths and local sizes of all the local arrays onto the root process
    MPI_Gather(&myLen, 1, MPI_INT, lengths, 1, MPI_INT, 0 , MPI_COMM_WORLD);
    MPI_Gather(&A.localNx, 1, MPI_INT, localNxs, 1, MPI_INT, 0 , MPI_COMM_WORLD);
    MPI_Gather(&A.localNy, 1, MPI_INT, localNys, 1, MPI_INT, 0 , MPI_COMM_WORLD);

    // Calculate the displacement vector from the lengths
    if(A.world_rank == 0) {
        for (int i = 1; i < A.world_size; ++i) {
            displs[i] = displs[i - 1] + lengths[i - 1];
        }
    }

    // Gather all the data into receivebuffer
    MPI_Gatherv(udata_, A.localNx * A.localNy, MPI_DOUBLE, ureceivebuffer, lengths, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(vdata_, A.localNx * A.localNy, MPI_DOUBLE, vreceivebuffer, lengths, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    // Arrange the received data on rank 0
    if (A.world_rank == 0) {
        int yskip = 0;
        int columndisp = 0;
        int totaldisp = 0;
        for (int p = 0; p < A.nPx; p++) { ///< skipping to the next process to the right
            for (int k = 0; k < localNxs[p * A.nPy]; k++) { ///< skipping to the next column
                columndisp = 0; ///< reset the column displacement counter when on a new column
                for (int j = 0; j < A.nPy; j++) { ///< iterating over the processes below to assemble the entire column
                    for (int i = 0; i < localNys[j]; i++) { ///< retrieving a column from a local array
                        ucombineddata_[i + yskip] = ureceivebuffer[i + columndisp + k * localNys[j] + totaldisp];
                        vcombineddata_[i + yskip] = vreceivebuffer[i + columndisp + k * localNys[j] + totaldisp];
                    }
                    yskip += localNys[j]; // increase the yskip counter
                    columndisp += localNys[j] * localNxs[p * A.nPy]; ///< increase the column displacement counter
                }
            }
            // Generate the displacement for switching processes in the x-direction
            for (int l = 0; l < A.nPy ; ++l) {
                totaldisp += localNys[l] * localNxs[p * A.nPy];
            }
        }
    }


    // Free the dynamic memory arrays used
    delete[] utoparraysend;
    delete[] utoparrayreceive;
    delete[] ubtmarraysend;
    delete[] ubtmarrayreceive;
    delete[] vtoparraysend;
    delete[] vtoparrayreceive;
    delete[] vbtmarraysend;
    delete[] vbtmarrayreceive;

    if(A.world_rank == 0) {
        delete[] ureceivebuffer;
        delete[] vreceivebuffer;
        delete[] lengths;
        delete[] displs;
        delete[] localNxs;
        delete[] localNys;
    }

}


// Member function that prints the current velocity fields (in udata_ and vdata_) to the file VelocityFields.txt
void Burgers::FilePrintVelFields(Model &A) {
    // open file
    std::ofstream vMyFile("VelocityFields.txt");

    // check if filestream is OK
    if (vMyFile.good()) {
        vMyFile << "Current u Velocity Field: \n";   // think about flipping this to correct y axis

        for (int i = 0; i < A.Ny; ++i) { ///< i is the row index
            for (int j = 0; j < A.Nx; ++j) { ///< j is the column index
                vMyFile << std::fixed << std::setprecision(4) << std::scientific << ucombineddata_[A.Ny * j + i] << ' ';
            }
            vMyFile << '\n';
        }

        vMyFile << "Current v Velocity Field: \n";

        for (int i = 0; i < A.Ny; ++i) { ///< i is the row index
            for (int j = 0; j < A.Nx; ++j) { ///< j is the column index
                vMyFile << std::fixed << std::setprecision(4) << std::scientific << vcombineddata_[A.Ny * j + i] << ' ';
            }
            vMyFile << '\n';
        }

        std::cout << "Velocity Field printed to VelocityFields.txt. \n";
    } else std::cout << ("File could not be opened for writing.");

    vMyFile.close();
}

// Member function that returns the energy of the velocity field
double Burgers::EnergyOfVelField(Model &A) {

    return 0.5 * (F77NAME(ddot)(A.Nx * A.Ny, ucombineddata_, 1, ucombineddata_, 1) + F77NAME(ddot)(A.Nx * A.Ny, vcombineddata_, 1, vcombineddata_, 1)) * A.dx * A.dy;

}








