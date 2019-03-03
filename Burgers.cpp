// Created by jwl16 on 18/02/19.

#include "Burgers.h"

// Constructor for Burgers that accepts an argument of class Model
Burgers::Burgers(Model &A) {

    xsize_ = A.Nx;
    ysize_ = A.Ny;

    // Use dynamic memory for data arrays - initialise values to 0 using ()
    udata_ = new double[xsize_ * ysize_]();
    vdata_ = new double[xsize_ * ysize_]();
    // Initialise two extra arrays to pass between during the
    udata_2 = new double[xsize_ * ysize_]();
    vdata_2 = new double[xsize_ * ysize_]();
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

    for (int i = 0; i < A.Nx; ++i) { ///< i is the column index
        for (int j = 0; j < A.Ny; ++j) { ///< j is the row index
            r = sqrt((A.x0 + i * A.dx) * (A.x0 + i * A.dx) + (A.y0 + j * A.dy) * (A.y0 + j * A.dy));
            if (r <= 1.0) {
                udata_[A.Ny * i + j] = 2.0 * pow(1.0 - r, 4) * (4.0 * r + 1);
                vdata_[A.Ny * i + j] = 2.0 * pow(1.0 - r, 4) * (4.0 * r + 1);
            }
        }
    }

}

void Burgers::DisplayuVelField(Model &A) {

    std::cout << "Current u Velocity Field: \n";

    for (int i = 0; i < A.Ny; ++i) { ///< i is the row index
        for (int j = 0; j < A.Nx; ++j) { ///< j is the column index
            std::cout << std::fixed << std::setprecision(3) << udata_[A.Ny * j + i] << ' ';
        }
        std::cout << '\n';
    }
}

void Burgers::DisplayvVelField(Model &A) {

    std::cout << "Current v Velocity Field: \n";

    for (int i = 0; i < A.Ny; ++i) { ///< i is the row index
        for (int j = 0; j < A.Nx; ++j) { ///< j is the column index
            std::cout << std::fixed << std::setprecision(3) << vdata_[A.Ny * j + i] << ' ';
        }
        std::cout << '\n';
    }
}

// Member function that integrates the velocity fields for the parameters specified in Model &A and with initial conditions from SetVelField
void Burgers::TimeIntegrateVelField(Model &A) {

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
        // Alternate between using data in udata_/vdata_ and udata_2/vdata_2
        switch (t % 2) {
            case 0 :
                for (int i = 1; i < A.Nx - 1; ++i) { ///< i is the column tracker
                    for (int j = 1; j < A.Ny - 1; ++j) { ///< j is the row tracker
                        udata_2[A.Ny * i + j] =
                                (const3 + bdt_dx * vdata_[A.Ny * i + j]) * udata_[A.Ny * i + (j - 1)] +
                                (const1 - bdt_dy * vdata_[A.Ny * i + j] +
                                 bdt_dx * (udata_[A.Ny * (i - 1) + j] - udata_[A.Ny * i + j])) *
                                udata_[A.Ny * i + j] + cdt_dysq * udata_[A.Ny * i + (j + 1)] +
                                const2 * udata_[A.Ny * (i - 1) + j] + cdt_dxsq * udata_[A.Ny * (i + 1) + j];

                        vdata_2[A.Ny * i + j] =
                                const3 * vdata_[A.Ny * i + (j - 1)] + (const1 - bdt_dx * udata_[A.Ny * i + j]
                                                                            + bdt_dy *
                                                                              (vdata_[A.Ny * i + (j - 1)] -
                                                                               vdata_[A.Ny * i + j])) *
                                                                           vdata_[A.Ny * i + j] +
                                cdt_dysq * vdata_[A.Ny * i + (j + 1)] +
                                (const2 + bdt_dx * udata_[A.Ny * i + j]) * vdata_[A.Ny * (i - 1) + j] +
                                cdt_dxsq * vdata_[A.Ny * (i + 1) + j];


                    }
                }

                break;
            case 1 :
                for (int i = 1; i < A.Nx - 1; ++i) { ///< i is the column tracker
                    for (int j = 1; j < A.Ny - 1; ++j) { ///< j is the row tracker
                        udata_[A.Ny * i + j] =
                                (const3 + bdt_dx * vdata_2[A.Ny * i + j]) * udata_2[A.Ny * i + (j - 1)] +
                                (const1 - bdt_dy * vdata_2[A.Ny * i + j]
                                 + bdt_dx * (udata_2[A.Ny * (i - 1) + j] - udata_2[A.Ny * i + j])) *
                                udata_2[A.Ny * i + j] +
                                cdt_dysq * udata_2[A.Ny * i + (j + 1)] +
                                const2 * udata_2[A.Ny * (i - 1) + j] +
                                cdt_dxsq * udata_2[A.Ny * (i + 1) + j];

                        vdata_[A.Ny * i + j] = const3 * vdata_2[A.Ny * i + (j - 1)] +
                                                    (const1 - bdt_dx * udata_2[A.Ny * i + j]
                                                     + bdt_dy * (vdata_2[A.Ny * i + (j - 1)] -
                                                                 vdata_2[A.Ny * i + j])) *
                                                    vdata_2[A.Ny * i + j] +
                                                    cdt_dysq * vdata_2[A.Ny * i + (j + 1)] +
                                                    (const2 + bdt_dx * udata_2[A.Ny * i + j]) *
                                                    vdata_2[A.Ny * (i - 1) + j] +
                                                    cdt_dxsq * vdata_2[A.Ny * (i + 1) + j];


                    }
                }
        

                break;

            default :
                std::cout << "Error in switching during time integration";
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








