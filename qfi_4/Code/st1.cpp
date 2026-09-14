#include <iostream>
#include <cmath>
#include <math.h> 
#include <complex>
#include <fstream>
//#include <omp.h>
//#include "mkl_lapacke.h"
#include "Master/master.hpp"
using namespace std;

char title[100];

int main(int argc,char* argv[])
{
    // std::ofstream file;               //| ==== File | L
    // file.open("U-phi.dat");

    int N=atoi(argv[1]);
    double T=atof(argv[2]);
    double J2=atof(argv[3]);

    input Data;
    Data.T=T;
    Data.J2=J2;
    Data.besis(N,0.27,0.0,0.0,0.0);
   // Data.Values();
    Data.mu_phi();
    
    // for (double i = 0; i <= 100; i+=1)
    // {
    //     Data.U_phi(i/10.0);
    //     file<< i/10.0<<"    "<<Data.O.transpose().real()<<endl;
    // }

    //int i=atof(argv[1]);
    // int u=i-8;
   
    // sprintf(title,"U-phi.dat"); //| ==== Name | I
    // std::ofstream file(title, std::fstream::out | std::fstream::app); 
    // {
    //     Data.mu_phi();
    //     //file<< u<<"    "<<Data.O.transpose().real()<<endl;
    //     file<<"\n"<<endl;
    // }


    // sprintf(title,"U-phi.dat"); //| ==== Name | I
    // std::ofstream file(title, std::fstream::out | std::fstream::app); 
    // {
    //     Data.U_phi(i/10.0);
    //     file<< i/10.0<<"    "<<Data.O.transpose().real()<<endl;
    //     file<<"\n"<<endl;
    // }



    //file.close();
}
