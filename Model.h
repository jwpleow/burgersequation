#ifndef CLASS_MODEL
#define CLASS_MODEL

//#include "mpi.h"

/**
* @class Model
* @brief Stores the relevant parameters for the Burgers equation
*/

class Model {
    friend class Burgers;
public:

    // * * * * * * * * * * * * * * CONSTRUCTOR * * * * * * * * * * * * * * * * //
    Model(int argc, char *argv[]);
    // * * * * * * * * * * * * * * DESTRUCTOR  * * * * * * * * * * * * * * * * //
    ~Model();
    // * * * * * * * * * * * * * * MEMBER FUNCTIONS  * * * * * * * * * * * * * //
    void PrintParameters();

    bool IsValid(); ///< think about writing this and calling immediately after initialisation?

    // bool   IsVerbose() const { return verbose; } ///< tells user what it does (Extra information)
    // bool   IsHelp()    const { return help; } ///< print out available command line options

    double GetX0() const { return x0; }
    double GetY0() const { return y0; }
    double GetLx() const { return Lx; }
    double GetLy() const { return Ly; }
    double GetT() const { return T; }
    int GetNx() const { return Nx; }
    int GetNy() const { return Ny; }
    int GetNt() const { return Nt; }
    double GetDx() const { return dx; }
    double GetDy() const { return dy; }

    double GetDt() const { return dt; }

    double GetAx() const { return ax; }

    double GetAy() const { return ay; }

    double GetB() const { return b; }

    double GetC() const { return c; }

    // Add any other getters here...

private:
    void ParseParameters(int argc, char *argv[]);

    void ValidateParameters();

    //bool verbose;
    //bool help;

    // Numerics
    double x0; ///< lowest x-value in the domain
    double y0; ///< lowest y-value in the domain
    double Lx; ///< x domain size (-Lx/2 <= x <= Lx/2)
    double Ly; ///< y domain size (-Ly/2 <= x <= Ly/2)
    double T;  ///< final time
    unsigned long Nx;    ///< number of x direction discretisations / grid points
    unsigned long Ny;    ///< number of y direction discretisations / grid points
    unsigned long Nt;    ///< number of timesteps
    double dx; ///< x direction discretisation
    double dy; ///< y direction discretisation
    double dt; ///< time discretisation

    // Physics
    double ax; ///< coefficient a_x
    double ay; ///< coefficient a_y
    double b;  ///< coefficient b
    double c;  ///< coefficient c

    // Add any additional parameters here...
};


#endif
