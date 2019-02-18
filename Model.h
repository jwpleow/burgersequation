#ifndef CLASS_MODEL
#define CLASS_MODEL

#include "mpi.h"

class Model {
	public:
		Model(int argc, char* argv[]);
        ~Model();

        void PrintParameters();

        bool IsValid();

        // Getters
		bool   IsVerbose() const { return verbose; }
        bool   IsHelp()    const { return help; }
        double GetX0()     const { return x0; }
        double GetY0()     const { return y0; }
        double GetLx()     const { return Lx; }
        double GetLy()     const { return Ly; }
        double GetT()      const { return T; }
        int    GetNx()     const { return Nx; }
        int    GetNy()     const { return Ny; }
        int    GetNt()     const { return Nt; }
        double GetDx()     const { return dx; }
        double GetDy()     const { return dy; }
        double GetDt()     const { return dt; }
        double GetAx()     const { return ax; }
        double GetAy()     const { return ay; }
        double GetB()      const { return b; }
        double GetC()      const { return c; }

        // Add any other getters here...

	private:
        void ParseParameters(int argc, char* argv[]);
        void ValidateParameters();

	    bool verbose;
        bool help;

	    // Numerics
	    double x0;
	    double y0;
	    double Lx;
	    double Ly;
	    double T;
	    int    Nx;
	    int    Ny;
	    int    Nt;
	    double dx;
	    double dy;
	    double dt;

	    // Physics
	    double ax;
	    double ay;
	    double b;
	    double c;

        // Add any additional parameters here...
};


#endif
