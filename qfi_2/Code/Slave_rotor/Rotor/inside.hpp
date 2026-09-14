#ifndef INSIDE_HPP_INCLUDED
#define INSIDE_HPP_INCLUDED
#include "../input.hpp"
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

//===================================================================================================//
                                    //=== Inside the cluster Hoping ===//
//===================================================================================================//


void input::insidehoping( int i)
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
            {
                if(N!=1)
                {
                    H(q,n)=J(k1,k2)*sss*0.5;
                    H(n,q)=conj(H(q,n));
    
                    // H(q,n)=t(k1,k2)*X(k1,k2)*sss;
                    // H(n,q)=conj(H(q,n));

                }
                else
                {
                    H(q,n)=J(k1,k2)*sss*0.5;
                    H(n,q)=conj(H(q,n));

                }
            }
        }
        
        i=n;
    }
}
//===========================================================================================//
//===========================================================================================//
#endif
