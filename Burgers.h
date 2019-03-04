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


// Define BLAS functions required
#define F77NAME(x) x##_
extern "C" {

// double dot product
double F77NAME(ddot)(const int& n, const double* dx, const int& incx, const double* dy, const int& incy);

}


class Burgers {
public:
    // * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

    Burgers(Model& B);

    // * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //

    ~Burgers();

    // * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //


    void SetVelField(); ///< Calculates the initial velocity field from the input parameters

    // Functions to display the u and v velocity field, with the top left corner
    // being u(-L/2,-L/2), and the bottom right corner u(L/2,L/2)
    void DisplayuVelField();
    void DisplayvVelField();

    void TimeIntegrateVelField(); ///< Function that time integrates the velocity fields according to Burgers'

    double EnergyOfVelField(); ///< Function that returns the total energy of the velocity fields

    void FilePrintVelFields(); ///< Function to print the velocity fields to a file "VelocityFields.txt"


    double* udata() const { return udata_; } ///< Returns Pointer to u data
    double* vdata() const { return vdata_; }; ///< Returns Pointer to v data

private:

    Model *A;
    double *udata_;
    double *vdata_;
    double *udata_2;
    double *vdata_2;
    double *ucombineddata_;
    double *vcombineddata_;
};


#endif
