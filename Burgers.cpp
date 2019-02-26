// Created by jwl16 on 18/02/19.

#include "Burgers.h"

// Constructor for Burgers that accepts an argument of class Model
Burgers::Burgers(Model &A) {

    xsize_ = A.Nx;
    ysize_ = A.Ny;

    // Use dynamic memory for data arrays - initialise values to 0 using ()
    udata_ = new double[A.localNx * A.localNy]();
    vdata_ = new double[A.localNx * A.localNy]();
    // Initialise two extra arrays to pass between during the time integration
    udata_2 = new double[A.localNx * A.localNy]();
    vdata_2 = new double[A.localNx * A.localNy]();
}

// Destructor
Burgers::~Burgers() {
    delete[] udata_;
    delete[] vdata_;
};

/**
 * Set the velocity field
 * @param pVal number of something
 */
void Burgers::SetVelField(Model &A) {

    double r; ///< Initialise r value

    // Create new placeholder array to store initial velocity data
    auto* uinitveldata_ = new double[A.Nx * A.Ny]();
    auto* vinitveldata_ = new double[A.Nx * A.Ny]();

    // Calculate the initial velocity field
    if (A.world_rank == 0) {
        for (int i = 0; i < A.Nx; ++i) { ///< i is the column index
            for (int j = 0; j < A.Ny; ++j) { ///< j is the row index
                r = sqrt((A.x0 + i * A.dx) * (A.x0 + i * A.dx) + (A.y0 + j * A.dy) * (A.y0 + j * A.dy));
                if (r <= 1.0) {
                    uinitveldata_[A.Ny * i + j] = 2.0 * pow(1.0 - r, 4) * (4.0 * r + 1);
                    vinitveldata_[A.Ny * i + j] = 2.0 * pow(1.0 - r, 4) * (4.0 * r + 1);
                }
            }
        }
    }

//    // check velocity field
//    if (A.world_rank==0) {
//        std::cout << "Current u Velocity Field: \n";
//
//        for (int i = 0; i < A.Ny; ++i) { ///< i is the row index
//            for (int j = 0; j < A.Nx; ++j) { ///< j is the column index
//                std::cout << std::fixed << std::setprecision(2) << uinitveldata_[A.Ny * j + i] << ' ';
//            }
//            std::cout << '\n';
//        }
//    }


    // broadcast the initial velocity data
    MPI_Bcast(uinitveldata_, A.Nx * A.Ny , MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(vinitveldata_, A.Nx * A.Ny , MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Split the initialised velocity data
    for (int i = 0; i < A.localNx; ++i){
        for (int j = 0; j < A.localNy; ++j){
            udata_[ i * A.localNy + j] = uinitveldata_[A.localstart + i * A.localNy + j];
            vdata_[ i * A.localNy + j] = vinitveldata_[A.localstart + i * A.localNy + j];
        }
    }



    // Delete the unneeded initial velocity arrays
    delete[] uinitveldata_;
    delete[] vinitveldata_;
}

void Burgers::DisplayuVelField(Model &A) {

    std::cout << "Current u Velocity Field: \n";

    for (int i = 0; i < A.localNy; ++i) { ///< i is the row index
        for (int j = 0; j < A.localNx; ++j) { ///< j is the column index
            std::cout << std::fixed << std::setprecision(2) << udata_[A.localNy * j + i] << ' ';
        }
        std::cout << '\n';
    }
}

void Burgers::DisplayvVelField(Model &A) {

    std::cout << "Current v Velocity Field: \n";

    for (int i = 0; i < A.localNy; ++i) { ///< i is the row index
        for (int j = 0; j < A.localNx; ++j) { ///< j is the column index
            std::cout << std::fixed << std::setprecision(2) << vdata_[A.localNy * j + i] << ' ';
        }
        std::cout << '\n';
    }
}

// Member function that integrates the velocity fields for the parameters specified in Model &A and with initial conditions from SetVelField
void Burgers::TimeIntegrateVelField(Model &A) {

//    MPI_Status status;
//    MPI_Request request;

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

    for (int t = 0; t < A.Nt; ++t) { ///< iterate over the timesteps
        for (int i = 1; i < A.localNx - 1; ++i) { ///< i is the column tracker
            for (int j = 1; j < A.localNy - 1; ++j) { ///< j is the row tracker


                // Alternate between using data in udata_/vdata_ and udata_2/vdata_2
                switch (t % 2) {
                    case 0 :
                        udata_2[A.localNy * i + j] = const3 * udata_[A.localNy * i + (j - 1)] + const1 * udata_[A.localNy * i + j] +
                                                cdt_dysq * udata_[A.localNy * i + (j + 1)] +
                                                bdt_dx * (udata_[A.localNy * (i - 1) + j] - udata_[A.localNy * i + j]) *
                                                udata_[A.localNy * i + j] +
                                                bdt_dy * (udata_[A.localNy * i + (j - 1)] - udata_[A.localNy * i + j]) *
                                                vdata_[A.localNy * i + j] + const2 * udata_[A.localNy * (i - 1) + j] +
                                                cdt_dxsq * udata_[A.localNy * (i + 1) + j];

                        vdata_2[A.localNy * i + j] = const3 * vdata_[A.localNy * i + (j - 1)] + const1 * vdata_[A.localNy * i + j] +
                                                cdt_dysq * vdata_[A.localNy * i + (j + 1)] +
                                                bdt_dx * (vdata_[A.localNy * (i - 1) + j] - vdata_[A.localNy * i + j]) *
                                                udata_[A.localNy * i + j] +
                                                bdt_dy * (vdata_[A.localNy * i + (j - 1)] - vdata_[A.localNy * i + j]) *
                                                vdata_[A.localNy * i + j] + const2 * vdata_[A.localNy * (i - 1) + j] +
                                                cdt_dxsq * vdata_[A.localNy * (i + 1) + j];



                        break;
                    case 1 :
                        udata_[A.localNy * i + j] = const3 * udata_2[A.localNy * i + (j - 1)] + const1 * udata_2[A.localNy * i + j] +
                                               cdt_dysq * udata_2[A.localNy * i + (j + 1)] +
                                               bdt_dx * (udata_2[A.localNy * (i - 1) + j] - udata_2[A.localNy * i + j]) *
                                               udata_2[A.localNy * i + j] +
                                               bdt_dy * (udata_2[A.localNy * i + (j - 1)] - udata_2[A.localNy * i + j]) *
                                               vdata_2[A.localNy * i + j] + const2 * udata_2[A.localNy * (i - 1) + j] +
                                               cdt_dxsq * udata_2[A.localNy * (i + 1) + j];

                        vdata_[A.localNy * i + j] = const3 * vdata_2[A.localNy * i + (j - 1)] + const1 * vdata_2[A.localNy * i + j] +
                                               cdt_dysq * vdata_2[A.localNy * i + (j + 1)] +
                                               bdt_dx * (vdata_2[A.localNy * (i - 1) + j] - vdata_2[A.localNy * i + j]) *
                                               udata_2[A.localNy * i + j] +
                                               bdt_dy * (vdata_2[A.localNy * i + (j - 1)] - vdata_2[A.localNy * i + j]) *
                                               vdata_2[A.localNy * i + j] + const2 * vdata_2[A.localNy * (i - 1) + j] +
                                               cdt_dxsq * vdata_2[A.localNy * (i + 1) + j];

                        break;
                    default :
                        std::cout << "Error in switching during time integration";
                        break;
                }
            }
        }

        switch (t % 2) {
            case 0 :
                if (A.world_rank == 0) {
                    MPI_Send(&udata_2[(A.Nx / 2 - 1) * A.localNy], A.localNy, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
                    MPI_Recv(&udata_2[(A.Nx / 2) * A.localNy], A.localNy, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                    MPI_Send(&vdata_2[(A.Nx / 2 - 1) * A.localNy], A.localNy, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD);
                    MPI_Recv(&vdata_2[(A.Nx / 2) * A.localNy], A.localNy, MPI_DOUBLE, 1, 3, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                    //MPI_Wait();
                } else {
                    MPI_Recv(&udata_2[0], A.localNy, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&udata_2[A.localNy], A.localNy, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
                    MPI_Recv(&vdata_2[0], A.localNy, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&vdata_2[A.localNy], A.localNy, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
                }
                break;
            case 1 :
                if (A.world_rank == 0) {
                    MPI_Send(&udata_[(A.Nx / 2 - 1) * A.localNy], A.localNy, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
                    MPI_Recv(&udata_[(A.Nx / 2) * A.localNy], A.localNy, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                    MPI_Send(&vdata_[(A.Nx / 2 - 1) * A.localNy], A.localNy, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD);
                    MPI_Recv(&vdata_[(A.Nx / 2) * A.localNy], A.localNy, MPI_DOUBLE, 1, 3, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                    //MPI_Wait();
                } else {
                    MPI_Recv(&udata_[0], A.localNy, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&udata_[A.localNy], A.localNy, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
                    MPI_Recv(&vdata_[0], A.localNy, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&vdata_[A.localNy], A.localNy, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
                }
                break;
            default :
                std::cout << "Error in switching during time integration - MPI communication";
                break;
        }

    } ///< timestep close bracket

    // Make sure that final data is always in udata_ and vdata_, and clear the memory for the unused case
    if (A.Nt % 2 == 1) {
        delete[] udata_;
        delete[] vdata_;
        udata_ = udata_2;
        vdata_ = vdata_2;
    } else {
        delete[] udata_2;
        delete[] vdata_2;
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
                vMyFile << std::fixed << std::setprecision(4) << std::scientific << udata_[A.Ny * j + i] << ' ';
            }
            vMyFile << '\n';
        }

        vMyFile << "Current v Velocity Field: \n";

        for (int i = 0; i < A.Ny; ++i) { ///< i is the row index
            for (int j = 0; j < A.Nx; ++j) { ///< j is the column index
                vMyFile << std::fixed << std::setprecision(4) << std::scientific << vdata_[A.Ny * j + i] << ' ';
            }
            vMyFile << '\n';
        }

        std::cout << "Velocity Field printed to VelocityFields.txt. \n";
    } else std::cout << ("File could not be opened for writing.");

    vMyFile.close();
}

// Member function that returns the energy of the velocity field
double Burgers::EnergyOfVelField(Model &A) {

    return 0.5 * (F77NAME(ddot)(A.Nx * A.Ny, udata_, 1, udata_, 1) + F77NAME(ddot)(A.Nx * A.Ny, vdata_, 1, vdata_, 1)) * A.dx * A.dy;

}








