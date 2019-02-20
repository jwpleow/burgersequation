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


    unsigned long xsize() const { return xsize_; } ///< Return x Size
    unsigned long ysize() const { return ysize_; } ///< Return y Size
    double* udata() const { return udata_; } ///< Returns Pointer to u data
    double* vdata() const { return vdata_; }; ///< Returns Pointer to v data
private:

    unsigned long xsize_; // long just in case of overflow from large size of matrix required
    unsigned long ysize_; // long just in case of overflow from large size of matrix required
    double* udata_;
    double* vdata_;

};


#endif
