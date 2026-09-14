#ifndef SLAVE_ROTOR_HPP_INCLUDED
#define SLAVE_ROTOR_HPP_INCLUDED
#include <math.h>
#include <complex>
#include "input.hpp"
#include "matworks.hpp"
#include "Rotor/Hspin.hpp"
#include "Rotor/QFI.hpp"
// #include "Rotor/spinconverge.hpp"
// #include "Rotor/oioj.hpp"

//#include "Dos/dos.hpp"
#include "values.hpp"
using namespace std;
using namespace Eigen;
char title3[100];
char title4[100];
char title5[100];

void input::besis(int NN1,double J0,double hx0, double hy0, double hz0)
{
    N=NN1;
    J1=J0;
    hx= hx0;
    hy= hy0;
    hz= hz0;

    Ls=2*lz+1;
    ls = 0;                 // Besis count start
    le = pow(Ls,N);         // Besis count End
    //==============================//
    H= MatrixXcd :: Zero(le,le);
    J= MatrixXcd :: Zero(N,N);
    t= MatrixXcd::Zero(N,N);
    evs = MatrixXcd :: Zero(le,le);
    es = VectorXd :: Zero(le);
    //==============================//
    Values();
    // for (int k=0;k<N;k++) //| === Basis loop (site basis)
    // {
    //     J(k,(k+1)%N)=J1;
    //     J((k+1)%N,k)=J1;
    //     t(k,(k+1)%N)=t1;
    //     t((k+1)%N,k)=t1;
    // }
}



//==========================================================================================//
void input::output()
{
//     cout<<"\nU: "<<U<<"\nPhi: "<< O.transpose()<<"\nRotor occupation: "<< Stotal<<endl;
}

//==========================================================================================//    

void input:: mu_phi()
{
    
    Hspin();  
    pair<MatrixXcd, VectorXd> e = Eigenspectrum(H);
    es=e.second;
    evs=e.first;

    //cout<< H <<endl;
    QFI(M_PI);
    
}
//=======================================================//
#endif
