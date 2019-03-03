/**
* @class Burgers
* @brief Stores relevant data and functions for solving the Burgers equation
*/

#ifndef CLASS_BURGERS
#define CLASS_BURGERS

#include "Model.h"
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <climits>


#define F77NAME(x) x##_
extern "C" {

// double dot product
double F77NAME(ddot)(const int& n, const double* dx, const int& incx, const double* dy, const int& incy);

}


class Burgers {

public:
    // * * * * * * * * * * * * * * CONSTRUCTOR * * * * * * * * * * * * * * * * //

    Burgers(Model& A);

    // * * * * * * * * * * * * * * DESTRUCTOR  * * * * * * * * * * * * * * * * //

    ~Burgers();

    // * * * * * * * * * * * * * * MEMBER FUNCTIONS  * * * * * * * * * * * * * //

    // Calculates the initial velocity field from the input parameters for the initial condition:
    // At t=0; r=sqrt(x^2+y^2); u & v = 2 * (1 - r) ^4 * (4 * r + 1) for r <= 1 and = 0 for r > 1
    void SetVelField(Model& A);

    // Functions to display the u and v velocity field in the output, with the top left corner
    // being u(-L/2,-L/2), and the bottom right corner u(L/2,L/2)
    void DisplayuVelField(Model &A);
    void DisplayvVelField(Model &A);


    void TimeIntegrateVelField(Model &A);

    double EnergyOfVelField(Model &A);

    // Function to print the velocity fields to a file
    void FilePrintVelFields(Model &A);


    int xsize() const { return xsize_; } ///< Return x Size
    int ysize() const { return ysize_; } ///< Return y Size
    double* udata() const { return udata_; } ///< Returns Pointer to u data
    double* vdata() const { return vdata_; }; ///< Returns Pointer to v data
private:

    int xsize_; // long just in case of overflow from large size of matrix required
    int ysize_; // long just in case of overflow from large size of matrix required
    double* udata_;
    double* vdata_;
    double* udata_2;
    double* vdata_2;

};


#endif
