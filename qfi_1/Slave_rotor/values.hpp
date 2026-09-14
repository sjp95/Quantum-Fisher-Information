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

   J= MatrixXcd:: Zero(N,N);
   t= MatrixXcd:: Zero(N,N);
    for (int k=0;k<N;k++) //| === Basis loop (site basis)
    {
       J(k,(k+1)%N)=J1;
       J((k+1)%N,k)=J1;
       t(k,(k+1)%N)=t1;
       t((k+1)%N,k)=t1;
    }

}
//======================================================//

#endif
