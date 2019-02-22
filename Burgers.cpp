//
// Created by jwl16 on 18/02/19.
//

#include "Burgers.h"
//burger in chinese is called han bao. Macdonalds in chinese is called mai dang lao. xie xie ni. zai jian.

// Constructor for Burgers that accepts an argument of class Model
Burgers::Burgers(Model &A) {

    // size of matrix required would be gridpoints * (timesteps+1)
    xsize_ = A.Nx;
    ysize_ = A.Ny;

    // Use dynamic memory for data arrays - initialise values to 0 using ()
    udata_ = new double[xsize_ * ysize_]();
    vdata_ = new double[xsize_ * ysize_]();

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
void Burgers::SetVelField(Model& A) {

    double r; ///< Initialise r value

    for (int i = 0; i < A.Nx; ++i) { ///< i is the column index
        for (int j = 0; j < A.Ny; ++j) { ///< j is the row index
            r = sqrt((A.x0 + i * A.dx) * (A.x0 + i * A.dx) + (A.y0 + j * A.dy) * (A.y0 + j * A.dy));
            if (r <= 1.0) {
                udata_[A.Ny * i + j] = 2.0 * pow(1.0-r,4) * (4.0*r + 1);
                vdata_[A.Ny * i + j] = 2.0 * pow(1.0-r,4) * (4.0*r + 1);
            }
        }
    }

}

void Burgers::DisplayuVelField(Model &A) {

    std::cout << "Current u Velocity Field: \n";

    for (int i = 0; i < A.Ny; ++i) { ///< i is the row index
        for (int j = 0; j < A.Nx; ++j) { ///< j is the column index
            std::cout << std::fixed << std::setprecision(2) << udata_[A.Ny * j + i] << ' ';
        }
        std::cout << '\n';
    }
}

void Burgers::DisplayvVelField(Model &A) {

    std::cout << "Current v Velocity Field: \n";

    for (int i = 0; i < A.Ny; ++i) { ///< i is the row index
        for (int j = 0; j < A.Nx; ++j) { ///< j is the column index
            std::cout << std::fixed << std::setprecision(2) << vdata_[A.Ny * j + i] << ' ';
        }
        std::cout << '\n';
    }
}

void Burgers::TimeIntegrateVelField(Model &A) {

    // Initialise placeholder values required for the Forward Euler Explicit Method
    double u_n;

    for (int t = 0; t < A.Nt; ++t) {
        for (int i = 1; i < A.Nx - 1; ++i) {
            for (int j = 1; j < A.Ny - 1; ++j) {

                u_n = udata_[A.Ny * i + j]; ///< Store the current value for use in determining v

                udata_[A.Ny * i + j] += A.dt * (A.c * ((udata_[A.Ny * (i + 1) + j] - 2.0 * udata_[A.Ny * i + j] +
                                                        udata_[A.Ny * (i - 1) + j]) / (A.dx * A.dx) +
                                                       (udata_[A.Ny * i + (j + 1)] - 2.0 * udata_[A.Ny * i + j] +
                                                        udata_[A.Ny * i + (j - 1)]) / (A.dy * A.dy))
                                                - (A.ay + A.b * vdata_[A.Ny * i + j]) *
                                                  ((udata_[A.Ny * i + j] - udata_[A.Ny * i + (j - 1)]) / A.dy)
                                                - (A.ax + A.b * udata_[A.Ny * i + j]) *
                                                  ((udata_[A.Ny * i + j] - udata_[A.Ny * (i - 1) + j]) / A.dx));
                vdata_[A.Ny * i + j] += A.dt * (A.c * ((vdata_[A.Ny * (i + 1) + j] - 2.0 * vdata_[A.Ny * i + j] +
                                                        vdata_[A.Ny * (i - 1) + j]) / (A.dx * A.dx) +
                                                       (vdata_[A.Ny * i + (j + 1)] - 2.0 * vdata_[A.Ny * i + j] +
                                                        vdata_[A.Ny * i + (j - 1)]) / (A.dy * A.dy))
                                                - (A.ay + A.b * vdata_[A.Ny * i + j]) *
                                                  ((vdata_[A.Ny * i + j] - vdata_[A.Ny * i + (j - 1)]) / A.dy)
                                                - (A.ax + A.b * u_n) *
                                                  ((vdata_[A.Ny * i + j] - vdata_[A.Ny * (i - 1) + j]) / A.dx));

            }
        }
    }
}


void Burgers::PrintVelFields(Model &A) {
    //open file
    std::ofstream vMyFile("VelocityFields.txt");
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

        std::cout << "Velocity Files printed to VelocityFields.txt.";
    } else throw std::runtime_error("File could not be opened for writing.");

    vMyFile.close();
}






double Burgers::EnergyOfVelField(Model &A) {
    double* usquare;
    double* vsquare;
    usquare = new double[A.Nx*A.Ny]();
    vsquare = new double[A.Nx*A.Ny]();
    double usquaresum, vsquaresum;


    // Perform the square of the u and v terms
    F77NAME(dsbmv)('U', A.Nx*A.Ny , 0, 1.0, udata_ , 1, udata_ , 1, 0, usquare, 1);
    F77NAME(dsbmv)('U', A.Nx*A.Ny , 0, 1.0, vdata_ , 1, vdata_ , 1, 0, vsquare, 1); // check the unsigned long -> signed int conv

    // Sum the terms
    usquaresum=F77NAME(dasum)(A.Nx*A.Ny , usquare , 1);
    vsquaresum=F77NAME(dasum)(A.Nx*A.Ny , vsquare , 1);

    delete[] usquare;
    delete[] vsquare;
    return 0.5*(usquaresum + vsquaresum)*A.dx*A.dy;

}








