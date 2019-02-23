#ifndef CLASS_BURGERS
#define CLASS_BURGERS

/**
* @class Burgers
* @brief Stores relevant data for solving the Burgers equation
*/
#include "Model.h"
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <climits>

#define F77NAME(x) x##_
extern "C" {
double F77NAME(dnrm2)(const int& n, const double* x, const int& incx);  // double Euclidean norm

// double absolute sum
double F77NAME(dasum)(const int& n, const double* x, const int& incx);

//banded matrix vector multiplication (y <- alpha*A*x + beta*y
void F77NAME(dsbmv)(const char& uplo, ///< specify whether upper or lower triangular
                   const int& n, ///< order of the matrix
                   const int& k, ///< number of super diagonals
                   const double& alpha, ///< scalar alpha
                   const double *A, ///< Symmetrical banded matrix A
                   const int& lda, ///< Leading dimension of A
                   const double *x, ///< Vector to multiply A by
                   const int& incx,
                   const double& beta,
                   double *y,
                   const int& incy);


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
    void PrintVelFields(Model &A);


    int xsize() const { return xsize_; } ///< Return x Size
    int ysize() const { return ysize_; } ///< Return y Size
    double* udata() const { return udata_; } ///< Returns Pointer to u data
    double* vdata() const { return vdata_; }; ///< Returns Pointer to v data
private:

    int xsize_; // long just in case of overflow from large size of matrix required
    int ysize_; // long just in case of overflow from large size of matrix required
    double* udata_;
    double* vdata_;


};


#endif
