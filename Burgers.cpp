//
// Created by jwl16 on 18/02/19.
//

#include "Burgers.h"


// Constructor for Burgers that accepts an argument of class Model
Burgers::Burgers(Model &A) {
size_ = A.Nx*A.Ny; // size of matrix required would be gridpoints * (timesteps+1)

}

// Destructor
Burgers::~Burgers() {};

/**
 * Set the velocity field
 * @param pVal number of something
 */
void Burgers::SetVelField() {

    if (size_ == 0)
        throw std::length_error("Cannot set velocity field for matrix of 0 size");

//    for (unsigned int i = 0; i < size_ * size_ * ; ++i)
//        data[i] = pVal;
//


}

void Burgers::TimeIntegrateVelField() {
}

// Function to return the size
long long int Burgers::size() const {
    return size_;
}

// Function to return the u data
double *Burgers::udata() const {
    return udata_;
}
// Function to return the v data
double *Burgers::vdata() const {
    return vdata_;
}


