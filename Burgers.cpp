// Created by jwl16 on 18/02/19.

#include "Burgers.h"

// Constructor for Burgers that accepts an argument of class Model
Burgers::Burgers(Model &B) {
    A = &B;
    // Use dynamic memory for data arrays - initialise values to 0 using ()
    udata_ = new double[A->GetLocalNx() * A->GetLocalNy()]();
    vdata_ = new double[A->GetLocalNx() * A->GetLocalNy()]();
    // Initialise extra arrays for switching between during the time integration
    udata_2 = new double[A->GetLocalNx() * A->GetLocalNy()]();
    vdata_2 = new double[A->GetLocalNx() * A->GetLocalNy()]();
}

// Destructor
Burgers::~Burgers() {
    // Free dynamic memory of the velocity field data on rank 0 process
    if (A->GetWorldRank() == 0) {
        delete[] ucombineddata_;
        delete[] vcombineddata_;
    }
}

// Member function to set the initial velocity field and distribute them to the processes
// Calculates the initial velocity field from the input parameters for the initial condition:
// At t=0; r=sqrt(x^2+y^2); u & v = 2 * (1 - r) ^4 * (4 * r + 1) for r <= 1 and = 0 for r > 1
void Burgers::SetVelField() {

    double r; ///< Initialise r value

    // Calculate the initial velocity field for the local process from its offset position in the global system
    double c1 = A->x0 + (A->GetLocalStart() / A->GetNy()) * A->GetDx(); ///< Precalculate some constants in the loop
    double c2 = A->y0 + (A->GetLocalStart() % A->GetNy()) * A->GetDy();
    for (int i = 0; i < A->GetLocalNx(); ++i) { ///< i is the x direction index
            for (int j = 0; j < A->GetLocalNy(); ++j) { ///< j is the y direction index
                r = sqrt((c1 + i * A->GetDx()) * (c1 + i * A->GetDx()) + (c2 + j * A->GetDy()) * (c2 + j * A->GetDy()));
                if (r <= 1.0) {
                    udata_[A->GetLocalNy() * i + j] = 2.0 * pow(1.0 - r, 4) * (4.0 * r + 1);
                    vdata_[A->GetLocalNy() * i + j] = 2.0 * pow(1.0 - r, 4) * (4.0 * r + 1);
                }
            }
        }

}


// Member functions to print to screen the current u/v velocity fields (use only for checking small discretisations!)
void Burgers::DisplayuVelField() {

    std::cout << "Current u Velocity Field: \n";

    for (int i = 0; i < A->GetNy(); ++i) { ///< i is the row index
        for (int j = 0; j < A->GetNx(); ++j) { ///< j is the column index
            std::cout << std::fixed << std::setprecision(3) << ucombineddata_[A->GetNy() * j + i] << ' ';
        }
        std::cout << '\n';
    }
}

void Burgers::DisplayvVelField() {

    std::cout << "Current v Velocity Field: \n";

    for (int i = 0; i < A->GetNy(); ++i) { ///< i is the row index
        for (int j = 0; j < A->GetNx(); ++j) { ///< j is the column index
            std::cout << std::fixed << std::setprecision(3) << vcombineddata_[A->GetNy() * j + i] << ' ';
        }
        std::cout << '\n';
    }
}

// Member function that integrates the velocity fields for the parameters specified in Model &A and with initial conditions from SetVelField
void Burgers::TimeIntegrateVelField() {

    //  * * * * * * * * * * Prepare arrays and constants * * * * * * * * * * * //

    // Initialise arrays to transfer the row edges for MPI
    auto *utoparraysend = new double[A->GetLocalNx()]();
    auto *utoparrayreceive = new double[A->GetLocalNx()]();
    auto *ubtmarraysend = new double[A->GetLocalNx()]();
    auto *ubtmarrayreceive = new double[A->GetLocalNx()]();
    auto *vtoparraysend = new double[A->GetLocalNx()]();
    auto *vtoparrayreceive = new double[A->GetLocalNx()]();
    auto *vbtmarraysend = new double[A->GetLocalNx()]();
    auto *vbtmarrayreceive = new double[A->GetLocalNx()]();

    // Make a temporary pointer for switching of data during time iteration
    double *utempptr = nullptr;
    double *vtempptr = nullptr;

    // Initialise precalculatable constants used in the for loops
    double cdt_dxsq, cdt_dysq, axdt_dx, aydt_dy, bdt_dx, bdt_dy, const1, const2, const3;
    int MPIconst1, MPIconst2;

    // Precalculate constants used in the explicit Forward Euler method
    cdt_dxsq = A->GetC() * A->GetDt() / (A->GetDx() * A->GetDx());
    cdt_dysq = A->GetC() * A->GetDt() / (A->GetDy() * A->GetDy());
    axdt_dx = A->GetAx() * A->GetDt() / A->GetDx();
    aydt_dy = A->GetAy() * A->GetDt() / A->GetDy();
    bdt_dx = A->GetB() * A->GetDt() / A->GetDx();
    bdt_dy = A->GetB() * A->GetDt() / A->GetDy();
    const1 = 1.0 - 2.0 * cdt_dxsq - 2.0 * cdt_dysq - axdt_dx - aydt_dy;
    const2 = (cdt_dxsq + axdt_dx);
    const3 = (cdt_dysq + aydt_dy);

    // Constants used for tags in MPI Send/Receives
    MPIconst1 = 4 * A->GetWorldRank() - (4 * A->GetPy());
    MPIconst2 = 4 * (A->GetWorldRank() - A->GetWorldRank() / A->GetPy()) + (A->GetPx() - 1) * 4 * A->GetPy();


    // * * * * * * * * * * * * * Time Integration Loop * * * * * * * * * * * * * * * //

    for (int t = 0; t < A->Nt; ++t) { ///< iterate over the timesteps
        for (int i = 1; i < A->GetLocalNx() - 1; ++i) { ///< i is the column tracker
            for (int j = 1; j < A->GetLocalNy() - 1; ++j) { ///< j is the row tracker
                udata_2[A->GetLocalNy() * i + j] =
                        (const3 + bdt_dx * vdata_[A->GetLocalNy() * i + j]) * udata_[A->GetLocalNy() * i + (j - 1)]
                        + (const1 - bdt_dy * vdata_[A->GetLocalNy() * i + j] +
                           bdt_dx * (udata_[A->GetLocalNy() * (i - 1) + j]
                                     - udata_[A->GetLocalNy() * i + j])) * udata_[A->GetLocalNy() * i + j]
                        + cdt_dysq * udata_[A->GetLocalNy() * i + (j + 1)] +
                        const2 * udata_[A->GetLocalNy() * (i - 1) + j] +
                        cdt_dxsq * udata_[A->GetLocalNy() * (i + 1) + j];

                vdata_2[A->GetLocalNy() * i + j] =
                        const3 * vdata_[A->GetLocalNy() * i + (j - 1)] +
                        (const1 - bdt_dx * udata_[A->GetLocalNy() * i + j]
                         + bdt_dy *
                           (vdata_[A->GetLocalNy() * i + (j - 1)] -
                            vdata_[A->GetLocalNy() * i + j])) *
                        vdata_[A->GetLocalNy() * i + j] +
                        cdt_dysq * vdata_[A->GetLocalNy() * i + (j + 1)] +
                        (const2 + bdt_dx * udata_[A->GetLocalNy() * i + j]) * vdata_[A->GetLocalNy() * (i - 1) + j] +
                        cdt_dxsq * vdata_[A->GetLocalNy() * (i + 1) + j];

            }
        }

        // Passing the edges to neighbouring processes
        // If not part of rightmost column - Send and receive data from the right
        // FIX: tags are unnecessarily complicated
        if (A->GetWorldRank() / A->GetPy() != A->GetPx() - 1) {
            MPI_Send(&udata_2[(A->GetLocalNx() - 2) * A->GetLocalNy()], A->GetLocalNy(), MPI_DOUBLE,
                     A->GetWorldRank() + A->GetPy(), 4 * A->GetWorldRank(),
                     MPI_COMM_WORLD); ///< Send rightside array to the rightside process
            MPI_Recv(&udata_2[(A->GetLocalNx() - 1) * A->GetLocalNy()], A->GetLocalNy(), MPI_DOUBLE,
                     A->GetWorldRank() + A->GetPy(), 4 * A->GetWorldRank() + 1, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE); ///< Receive rightside array from the right
            MPI_Send(&vdata_2[(A->GetLocalNx() - 2) * A->GetLocalNy()], A->GetLocalNy(), MPI_DOUBLE,
                     A->GetWorldRank() + A->GetPy(), 4 * A->GetWorldRank() + 2,
                     MPI_COMM_WORLD); ///< Send rightside array to the rightside process
            MPI_Recv(&vdata_2[(A->GetLocalNx() - 1) * A->GetLocalNy()], A->GetLocalNy(), MPI_DOUBLE,
                     A->GetWorldRank() + A->GetPy(), 4 * A->GetWorldRank() + 3, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE); ///< Receive rightside array from the right
        }
        // If not part of leftmost column - Send and receive data from the left
        if (A->GetWorldRank() / A->GetPy() != 0) {
            MPI_Recv(&udata_2[0], A->GetLocalNy(), MPI_DOUBLE, A->GetWorldRank() - A->GetPy(), MPIconst1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE); ///< Receive leftside array from the leftside process
            MPI_Send(&udata_2[A->GetLocalNy()], A->GetLocalNy(), MPI_DOUBLE, A->GetWorldRank() - A->GetPy(),
                     MPIconst1 + 1, MPI_COMM_WORLD); ///< Send leftside array to the leftside process
            MPI_Recv(&vdata_2[0], A->GetLocalNy(), MPI_DOUBLE, A->GetWorldRank() - A->GetPy(), MPIconst1 + 2,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE); ///< Receive leftside array from the leftside process
            MPI_Send(&vdata_2[A->GetLocalNy()], A->GetLocalNy(), MPI_DOUBLE, A->GetWorldRank() - A->GetPy(),
                     MPIconst1 + 3, MPI_COMM_WORLD); ///< Send leftside array to the leftside process

        }
        //If not part of topmost row - Send and receive data from the top
        if (A->GetWorldRank() % A->GetPy() != 0) {
            // Prepare an array of values from the 2nd row from the top
            for (int i = 0; i < A->GetLocalNx(); ++i) {
                utoparraysend[i] = udata_2[i * A->GetLocalNy() + 1];
                vtoparraysend[i] = vdata_2[i * A->GetLocalNy() + 1];
            }
            MPI_Send(utoparraysend, A->GetLocalNx(), MPI_DOUBLE, A->GetWorldRank() - 1, MPIconst2 - 4,
                     MPI_COMM_WORLD); ///< Send top array to the process above
            MPI_Recv(utoparrayreceive, A->GetLocalNx(), MPI_DOUBLE, A->GetWorldRank() - 1, MPIconst2 - 3,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE); ///< Receive the top array from the process above
            MPI_Send(vtoparraysend, A->GetLocalNx(), MPI_DOUBLE, A->GetWorldRank() - 1, MPIconst2 - 2,
                     MPI_COMM_WORLD); ///< Send top array to the process above
            MPI_Recv(vtoparrayreceive, A->GetLocalNx(), MPI_DOUBLE, A->GetWorldRank() - 1, MPIconst2 - 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE); ///< Receive the top array from the process above
            // After receiving, parse the data into the top row
            for (int i = 0; i < A->GetLocalNx(); ++i) {
                udata_2[i * A->GetLocalNy()] = utoparrayreceive[i];
                vdata_2[i * A->GetLocalNy()] = vtoparrayreceive[i];
            }
        }
        // If not part of bottom row - Send and receive data from the bottom
        if (A->GetWorldRank() % A->GetPy() != A->GetPy() - 1) {
            // Prepare an array of values from the 2nd row from the bottom
            for (int i = 0; i < A->GetLocalNx(); ++i) {
                ubtmarraysend[i] = udata_2[i * A->GetLocalNy() + (A->GetLocalNy() - 2)];
                vbtmarraysend[i] = vdata_2[i * A->GetLocalNy() + (A->GetLocalNy() - 2)];
            }
            MPI_Recv(ubtmarrayreceive, A->GetLocalNx(), MPI_DOUBLE, A->GetWorldRank() + 1, MPIconst2, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            MPI_Send(ubtmarraysend, A->GetLocalNx(), MPI_DOUBLE, A->GetWorldRank() + 1, MPIconst2 + 1, MPI_COMM_WORLD);
            MPI_Recv(vbtmarrayreceive, A->GetLocalNx(), MPI_DOUBLE, A->GetWorldRank() + 1, MPIconst2 + 2,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(vbtmarraysend, A->GetLocalNx(), MPI_DOUBLE, A->GetWorldRank() + 1, MPIconst2 + 3, MPI_COMM_WORLD);
            // After receiving, parse the data into the bottom row
            for (int i = 0; i < A->GetLocalNx(); ++i) {
                udata_2[i * A->GetLocalNy() + (A->GetLocalNy() - 1)] = ubtmarrayreceive[i];
                vdata_2[i * A->GetLocalNy() + (A->GetLocalNy() - 1)] = vbtmarrayreceive[i];
            }
        }

        // Swap pointers so that spatial loops are always using new time data
        utempptr = udata_;
        udata_ = udata_2;
        udata_2 = utempptr;
        vtempptr = vdata_;
        vdata_ = vdata_2;
        vdata_2 = vtempptr;

    } ///< timestep close bracket

    // Frees the memory for the unused case and temporary arrays which are not required anymore
    delete[] udata_2;
    delete[] vdata_2;
    delete[] utoparraysend;
    delete[] utoparrayreceive;
    delete[] ubtmarraysend;
    delete[] ubtmarrayreceive;
    delete[] vtoparraysend;
    delete[] vtoparrayreceive;
    delete[] vbtmarraysend;
    delete[] vbtmarrayreceive;

    // * * * * * * * * * * * * * Assembling the velocity field matrix * * * * * * * * * * * * * * //

    // First, 'remove' the edge columns/rows so that there is no overlapping data.
    // This is done by simply shortening localNx and localNy and placing the
    // relevant data (without the overlapping edges) in the 'top left' of the array.

    // If not part of rightmost column - remove rightmost column
    if (A->GetWorldRank() / A->GetPy() != A->GetPx() - 1) {
        A->localNx -= 1;
    }
    // If not part of leftmost column - remove leftmost column
    if (A->GetWorldRank() / A->GetPy() != 0) {
        A->localNx -= 1;
        for (int i = 0; i < A->GetLocalNx(); ++i) {
            for (int j = 0; j < A->GetLocalNy(); ++j) {
                udata_[j + A->GetLocalNy() * i] = udata_[j + (i + 1) * A->GetLocalNy()];
                vdata_[j + A->GetLocalNy() * i] = vdata_[j + (i + 1) * A->GetLocalNy()];
            }
        }
    }
    //If not part of topmost row - remove topmost row
    if (A->GetWorldRank() % A->GetPy() != 0) {
        for (int i = 0; i < A->GetLocalNx(); ++i) {
            for (int j = 0; j < A->GetLocalNy() - 1; ++j) {
                udata_[j + (A->GetLocalNy() - 1) * i] = udata_[j + 1 + i * A->GetLocalNy()];
                vdata_[j + (A->GetLocalNy() - 1) * i] = vdata_[j + 1 + i * A->GetLocalNy()];
            }
        }
        A->localNy -= 1;
    }
    // If not part of bottom row - remove bottommost row
    if (A->GetWorldRank() % A->GetPy() != A->GetPy() - 1) {
        A->localNy -= 1;
        for (int i = 0; i < A->GetLocalNx(); ++i) {
            for (int j = 0; j < A->GetLocalNy(); ++j) {
                udata_[j + A->GetLocalNy() * i] = udata_[j + i * (A->GetLocalNy() + 1)];
                vdata_[j + A->GetLocalNy() * i] = vdata_[j + i * (A->GetLocalNy() + 1)];
            }
        }
    }

    // Create an array to receive the data from all processes
    double *ureceivebuffer = nullptr;
    double *vreceivebuffer = nullptr;

    // Create the receive counts, displacements, and local size vectors required for MPI_Gatherv and rearranging
    int myLen = A->GetLocalNx() * A->GetLocalNy();
    int *lengths = nullptr;
    int *displs = nullptr;
    int *localNys = nullptr;
    int *localNxs = nullptr;

    // Allocate dynamic memory to hold the information required for assembly in the rank 0 process
    if (A->GetWorldRank() == 0) {
        lengths = new int[A->GetWorldSize()]();
        displs = new int[A->GetWorldSize()]();
        localNys = new int[A->GetWorldSize()]();
        localNxs = new int[A->GetWorldSize()]();
        ureceivebuffer = new double[A->GetNx() * A->GetNy()];
        vreceivebuffer = new double[A->GetNx() * A->GetNy()];
    }

    // Gather the lengths and local sizes of all the local arrays onto the root process
    MPI_Gather(&myLen, 1, MPI_INT, lengths, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&A->localNx, 1, MPI_INT, localNxs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&A->localNy, 1, MPI_INT, localNys, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Calculate the displacement vector from the lengths
    if (A->GetWorldRank() == 0) {
        for (int i = 1; i < A->GetWorldSize(); ++i) {
            displs[i] = displs[i - 1] + lengths[i - 1];
        }
    }

    // Gather all the data into receivebuffer
    MPI_Gatherv(udata_, A->GetLocalNx() * A->GetLocalNy(), MPI_DOUBLE, ureceivebuffer, lengths, displs, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
    MPI_Gatherv(vdata_, A->GetLocalNx() * A->GetLocalNy(), MPI_DOUBLE, vreceivebuffer, lengths, displs, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);

    // Free the local data arrays used which are not required anymore
    delete[] udata_;
    delete[] vdata_;

    // Arrange the received data on rank 0
    // TODO: this doesn't seem to work 100% of the time for nPx != nPy(?)
    if (A->GetWorldRank() == 0) {
        ucombineddata_ = new double[A->GetNx() * A->GetNy()](); ///< Create the arrays to store arranged combined velocity data
        vcombineddata_ = new double[A->GetNx() * A->GetNy()]();
        int yskip = 0;
        int columndisp = 0;
        int totaldisp = 0;
        for (int p = 0; p < A->GetPx(); p++) { ///< skipping to the next process to the right
            for (int k = 0; k < localNxs[p * A->GetPy()]; k++) { ///< moving to the next column
                columndisp = 0; ///< reset the column displacement counter when on a new column
                for (int j = 0; j < A->GetPy(); j++) { ///< iterating over the processes below to assemble the entire column
                    for (int i = 0; i < localNys[j]; i++) { ///< retrieving a column from a local array
                        ucombineddata_[i + yskip] = ureceivebuffer[i + columndisp + k * localNys[j] + totaldisp];
                        vcombineddata_[i + yskip] = vreceivebuffer[i + columndisp + k * localNys[j] + totaldisp];
                    }
                    yskip += localNys[j]; // increase the yskip counter
                    columndisp += localNys[j] * localNxs[p * A->GetPy()]; ///< increase the column displacement counter
                }
            }
            // Generate the displacement for switching processes in the x-direction
            for (int l = 0; l < A->GetPy(); ++l) {
                totaldisp += localNys[l] * localNxs[p * A->GetPy()];
            }
        }
    }

    // Free the dynamic memory arrays used which are not required anymore for the arranging of data on rank 0
    if (A->GetWorldRank() == 0) {
        delete[] lengths;
        delete[] displs;
        delete[] localNxs;
        delete[] localNys;
        delete[] ureceivebuffer;
        delete[] vreceivebuffer;
    }

}


// Member function that prints the current velocity fields (in ucombineddata_ and vcombineddata_) to the file VelocityFields.txt
// Can only call this function after time integration
void Burgers::FilePrintVelFields() {

    // Only the rank 0 process with the combined data should print
    if (A->GetWorldRank() == 0) {

        // open file
        std::ofstream vMyFile("VelocityFields.txt");

        // check if filestream is OK
        if (vMyFile.good()) {
            vMyFile << "Current u Velocity Field: \n";   // think about flipping this to correct y axis

            for (int i = 0; i < A->GetNy(); ++i) { ///< i is the row index
                for (int j = 0; j < A->GetNx(); ++j) { ///< j is the column index
                    vMyFile << std::fixed << std::setprecision(4) << std::scientific
                            << ucombineddata_[A->GetNy() * j + i]
                            << ' ';
                }
                vMyFile << '\n';
            }

            vMyFile << "Current v Velocity Field: \n";

            for (int i = 0; i < A->GetNy(); ++i) { ///< i is the row index
                for (int j = 0; j < A->GetNx(); ++j) { ///< j is the column index
                    vMyFile << std::fixed << std::setprecision(4) << std::scientific
                            << vcombineddata_[A->GetNy() * j + i]
                            << ' ';
                }
                vMyFile << '\n';
            }

            std::cout << "Velocity Field printed to VelocityFields.txt. \n";
        } else std::cout << ("File could not be opened for writing.");

        vMyFile.close();
    }
}

// Member function that returns the energy of the velocity field - Can only call this function after time integration
double Burgers::EnergyOfVelField() {
        return 0.5 * (F77NAME(ddot)(A->GetNx() * A->GetNy(), ucombineddata_, 1, ucombineddata_, 1) +
                      F77NAME(ddot)(A->GetNx() * A->GetNy(), vcombineddata_, 1, vcombineddata_, 1)) * A->GetDx() *
               A->GetDy();

}
