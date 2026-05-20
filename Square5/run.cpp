#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include "Master/master.hpp"
using namespace std;

// Usage: ./run  N  J1*100  J2*100  hz*100  Tmin  Tmax  nT

int main(int argc, char* argv[])
{
    int    N    = atoi(argv[1]);
    double J1   = atof(argv[2]);
    double J2   = atof(argv[3]);
    double hz   = atof(argv[4]);
    double Tmin = atof(argv[5]);
    double Tmax = atof(argv[6]);
    int    nT   = atoi(argv[7]);

    // cout<<"start\n";
    input Data;
    Data.besis(N, J1, J2, 0.0, 0.0, hz);
    Data.mu_phi(true);   // diagonalize once, cache eigenspectrum

    // M is built or loaded once inside QFI, then
    // the full T sweep runs in memory — no repeated I/O
    // cout<<"to enter tQFI\n";
    Data.tQFI(M_PI, Tmin, Tmax, nT, true);

    return 0;
}
