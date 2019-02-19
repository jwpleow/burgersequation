#ifndef CLASS_BURGERS
#define CLASS_BURGERS

/**
* @class Burgers
* @brief Stores relevant data for solving the Burgers equation
*/
#include "Model.h"
#include <stdexcept>
#include <iostream>

class Burgers {
public:
    Burgers(Model& A);    ///< Constructor
    ~Burgers();     ///< Destructor

    void SetVelField();
    void TimeIntegrateVelField();


    long long int size() const; ///< Return Size
    double* udata() const; ///< Returns Pointer to u data
    double* vdata() const; ///< Returns Pointer to v data
private:

    long long int size_; // long just in case of overflow from large size of matrix required
    double* udata_;
    double* vdata_;

};


#endif
