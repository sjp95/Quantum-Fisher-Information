#ifndef OIOJ_HPP_INCLUDED
#define OIOJ_HPP_INCLUDED
#include "../input.hpp"
#include "o.hpp"
#include <Eigen/Dense>
#include <omp.h>
using namespace std;
using namespace Eigen;

//===================================================================================================//
                                    //=== Inside the cluster Hoping ===//
//===================================================================================================//

//================================================================//
                            //== Up ==//
//================================================================//

void input::b_i_j( int i)
{
    int n = i;
    for(int k1=0; k1<N;k1++)
    for(int k2=k1+1; k2<N;k2++)
    {
        int id1 =0, id2=0;
        long long bi=0;
        for ( int j =0; j<N; j++)
        {
            int reminder = i%Ls;
            i=i/Ls;
            //======================//
            if(j==k1)
            {
                
                if (reminder>0)
                {
                    reminder-=1;
                    id1++;
                }
            }
            if(j==k2)
            {
                
                if (reminder<Ls-1)
                {
                    reminder+=1;
                    id2++;
                }

            }
            //======================//
            bi+=reminder*pow(Ls,j);
        }
        double sss=double(id1*id2);
        int q = int (bi);
        if(q!=n)
        {
            if(sss>0)
            B(k1,k2)+=conj(evs(n,0)*conj(evs(q,0))*sss);
        }
        i=n;
    }
}
    void input::Bij()
    {
        //===============================================//
        for(int k1=0; k1<N;k1++)
        {
            for(int k2=k1+1; k2<N;k2++)
            {
              B(k1,k2)=complex<double>(0,0);
            }
        }
        //===============================================//
        #pragma omp parallel for
        for(int i=0;i<le;i++)
        {
            b_i_j(i);
        }
        //===============================================//
        for(int k1=0; k1<N;k1++)
        for(int k2=k1+1; k2<N;k2++)
        {
            B(k2,k1)=conj(B(k1,k2));
        }
        //===============================================//
    }



//===========================================================================================//
//===========================================================================================//
#endif