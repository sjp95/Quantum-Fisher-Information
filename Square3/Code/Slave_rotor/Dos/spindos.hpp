#ifndef SPINDOS_HPP_INCLUDED
#define SPINDOS_HPP_INCLUDED
//===============================================================//
//===============================================================//
#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include <math.h>
#include <fstream>
#include <vector>
#include "../input.hpp"
#include "../matworks.hpp"
#include "../Rotor/Hspin.hpp"

using namespace std;
using namespace Eigen;

typedef std::pair<Eigen::MatrixXcd, Eigen::VectorXd> Spectrum;

//============================================================================================//
                             //=== <m|O|n> calculation ===//
//============================================================================================//

char title_s_d[100];
//================================================================//
                            //==Up==//
//================================================================//
void input::rotor( int i ,int fm, int fn)
{
    int n = i;
    for(int k1=0; k1<N; k1++)
    {
        long long bi1=0;
        int id=0;
        for ( int j =0; j<N; j++)
        {
            int reminder = i%Ls;
            i=i/Ls;
            //======================//
            if(j==k1)
            {
                 if (reminder<Ls-1)
                {
                    reminder+=1;
                    id++;
                }
            }
            //======================//
            bi1+=reminder*pow(Ls,j);
            //=======================//
        }
        //======================================//
        double sss=double(id);
        int q =int(bi1);// Decimal(bi1);
        if(q!=n)
        if(sss>0)
            DO_s(k1)+=conj(evss(q,fm))*evss(n,fn);
        //======================================//
        i=n;
    }
}
//===================================================================//
           //................ ********** .................//
//===================================================================//
 void input::Spin_dos()
 {
    evss = MatrixXcd :: Zero(le,le);
    DO_s = VectorXcd :: Zero(N);
    //------------------------------------------------------------//
    //------------------------------------------------------------//
     int tut=int(U*1000);                         //| ==== File | F
    //  int tut0=int(A*1000);
     sprintf(title_s_d,"../Data/Spin_%d.dat",tut); //| ==== File | F
     std::ofstream file;                          //| ==== Name | I
     file.open(title_s_d);                        //| ==== File | L
    //------------------------------------------------------------//
    //------------------------------------------------------------//
     Hspin();
     Spectrum E0 = Eigenspectrum(H);
     evss=E0.first;
    //  cout<<evss<<endl;
     //........................................//
     double Z=0; //|--------- Partition function
     for(int i=0; i<E0.second.size(); i++)
     {
         Z+=exp(-(E0.second(i)-E0.second(0))/T); //|------- Partition function calculation
     }
     //........................................//
     //=========================================================================//
     for(int fm=0;fm<le;fm++)
     for(int fn=0;fn<le;fn++)
     {
        DO_s = VectorXcd :: Zero(N);
        double pp1= exp((E0.second(0)-E0.second(fm))/T)/Z, pp2= exp((E0.second(0)-E0.second(fn))/T)/Z ;
        if(pp1 >pow(10,-6)|| pp2 >pow(10,-6))
         {
        //     /*--------------------------*/
            for(int i=0;i<le;i++)
            {
                rotor(i,fm,fn);   //|---- <m|O|n> for up spin
            }
            /*--------------------------*/
            double aaa=0;
            VectorXcd DO0_s = VectorXcd :: Zero(N);
            for (int i = 0; i < DO_s.size(); i++)
            {
                aaa+=pow(abs(DO_s(i)),2);
                DO0_s(i%N)+=0.5*DO_s(i);
            }
            aaa=aaa/DO_s.size();
            if(aaa>pow(10,-4))
            file<< DO0_s.transpose()<<"  "<< pp1 <<"  "<< pp2<<"  "<< E0.second(fm) <<"  "<< E0.second(fn)<<endl;
            // cout<<"spindos"<< DO0_s.transpose()<<"  "<< pp1 <<"  "<< pp2<<"  "<< E0.second(fm) <<"  "<< E0.second(fn)<<endl;

        }
     }
    file.close();
 }

//========================================================================================//
//========================================================================================//


#endif

//===============================================================//
//===============================================================//
