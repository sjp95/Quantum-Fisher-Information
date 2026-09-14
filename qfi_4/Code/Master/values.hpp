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

   for(int i=0;i<N;i+=3)
   {
      if(i%3==0)
      {
            //==============// i
            //==============// i+1
            // t(i,i+1)=t1;
            // t(i+1,i)=t1;

            J(i,i+1)=J1;
            J(i+1,i)=J1;
            //==============// i+1
            //==============// i+2
            // t(i+1,i+2)=t1;
            // t(i+2,i+1)=t1;

            J(i+1,i+2)=J2;
            J(i+2,i+1)=J2;

            //==============// i+1
            //==============// i+3
            // if (i+3<N)
            // {
               // t(i+2,(i+3)%N)=t1;
               // t((i+3)%N,i+2)=t1;

               J(i+2,(i+3)%N)=J1;
               J((i+3)%N,i+2)=J1;
            // }
            
            //==============// i+2
            //==============// i+3+1
            // tt(i+2,(i+4)%N)=t2;
            // tt((i+4)%N,i+2)=t2;

            // JJ(i+2,(i+4)%N)=J2;
            // JJ((i+4)%N,i+2)=J2;
      }
      
   
   }

}
//======================================================//

#endif
