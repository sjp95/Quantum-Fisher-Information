#ifndef VALUES_HPP_INCLUDED
#define VALUES_HPP_INCLUDED
#include <math.h>
#include "input.hpp"
#include <Eigen/Dense>

#include <complex>
using namespace Eigen;
using namespace std;

//=================== General input ====================//


void input::Values()
{

  Jx=MatrixXcd:: Zero(N,N);
  Jz=MatrixXcd:: Zero(N,N);

  // Trimer chain: 1D linear chain with J1-J1-J2 repeating pattern
  // Example for N=6: [0]--J1--[1]--J1--[2]--J2--[3]--J1--[4]--J1--[5]
  for (int i=0; i<N-1; i++)
  {
    int bond_index = i % 3;
    double J_coupling = (bond_index == 2) ? J2 : J1;

    Jx(i, i+1) = J_coupling;
    Jx(i+1, i) = J_coupling;
    Jz(i, i+1) = J_coupling;
    Jz(i+1, i) = J_coupling;
  }
  Jx(0,N-1) = J2;
  Jx(N-1,0) = J2;
  Jz(0,N-1) = J2;
  Jz(N-1,0) = J2;

  cout<< "==============================="<<endl;
  cout<< "Trimer Chain Couplings Generated"<<endl;
  cout<< "Chain pattern: J1--J1--J2 (repeating)"<<endl;
  cout<< "System size N: " << N << " (must be divisible by 3)"<<endl;
  cout<< "J1 (within-trimer): " << J1 << endl;
  cout<< "J2 (between-trimer): " << J2 << endl;
  cout<< "==============================="<<endl;

}
//======================================================//

#endif
