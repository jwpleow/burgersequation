#ifndef CLASS_BURGERS
#define CLASS_BURGERS

/**
* @class Burgers
* @brief Stores relevant data for solving the Burgers equation
*/

class Burgers{

            public:
            Complex(double re = 0.0, double im = 0.0);    //constructor
            ~Complex();     //destructor
            double GetRe();         //function declarations
            double GetIm();
            void SetVelField(double pVal);
            double arg();
            double magn();
            friend std::ostream& operator<<(std::ostream& os, const Complex& pComp);
            Complex& operator=(const Complex&);
            bool operator==(const Complex&);
            bool operator==(double pVal);
            friend bool operator==(const double pLhs, const Complex& pRhs);
            bool operator!=(const Complex& pComp);
            Complex operator+(const Complex& pComp);
            Complex operator-(const Complex& pComp);
            Complex operator*(const Complex& pComp);
            private:
            double a;       // Data: a is the real part, b is the imaginary part

};















#endif
