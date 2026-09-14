#ifndef FERMIDOS_HPP_INCLUDED
#define FERMIDOS_HPP_INCLUDED
#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include <math.h>
#include <fstream>
#include <vector>
#include "../input.hpp"
#include "../matworks.hpp"
#include "../Fermi/Hfermi.hpp"


using namespace std;
using namespace Eigen;
char title_f_d[100];

typedef std::pair<Eigen::MatrixXcd, Eigen::VectorXd> Spectrum;
//===============================================================//
//===============================================================//
void input::Fermi_dos()
{
    //------------------------------------------------------------//
    //------------------------------------------------------------//
     int tut=int(U*1000);                         //| ==== File | F
    //  int tut0=int(A*1000);
     sprintf(title_f_d,"../Data/Fermi_%d.dat",tut); //| ==== File | F
     std::ofstream file;                           //| ==== Name | I
     file.open(title_f_d);                         //| ==== File | L
    //------------------------------------------------------------//
    //------------------------------------------------------------//
   
    //=======================================================================================//
    //=======================================================================================//
    Hfermi();
    //===============================================//
    pair<MatrixXcd, VectorXd> efo = Eigenspectrum(Hf); //=== Eigen value and eigen vector calculation
    MatrixXcd evf = efo.first;   //| === Stortion of Eigen vector in globel array
    VectorXd ef =efo.second;    //| === Stortion of Eigen value in globel array
    //===============================================//
    VectorXd RR = VectorXd :: Zero(N);
    //---------------------------------------------------------------------------------------------//
    for(int f=0;f<ef.size();f++)
    { 
        double aaa=0;
        for (int k=0;k<Nf*Nf;k++)
        {
            int x=int(k%Nf);  //---x vatue
            int y=k/Nf;      //---y value
            //=======================================================================//
            //=======================================================================//
            int p=(y+1)%Nf;  //-- nearest neighbour along y
            int q=(x+p*Nf)%(Nf*Nf); //---- y axis hoping point
            //---------------------------------------------------//
            int r=(x+1)%Nf; //-- nearest neighbour along x
            int s=(r+y*Nf)%(Nf*Nf); //---- x axis hoping point 
            //=========================================================//




            //=========================================================//
            if(N==16)
            {
                int y1=y%4;
                int x1=x%4;

                int a=(x1+y1*4)%N;

                int p1=(y1+1)%4;
                int q1=(x1+p1*4)%N;

                int r1=(x1+1)%4;
                int s1=(r1+y1*4)%N;

                RR(a)+= pow(abs(evf(k,f)),2);
            }
            //=======================================================================//






            //=========================================================//
            if(N==8)
            {
                int y1=y%2;
                int x1=x%4;

                int a=(x1+y1*4)%N;

                int p1=(y1+1)%2;
                int q1=(x1+p1*4)%N;

                int r1=(x1+1)%4;
                int s1=(r1+y1*4)%N;

                RR(a)+= pow(abs(evf(k,f)),2);
            }
            //=======================================================================//



            //=========================================================//
            if(N==4)
            {
                int y1=y%2;
                int x1=x%2;

                int a=(x1+y1*2)%N;

                int p1=(y1+1)%2;
                int q1=(x1+p1*2)%N;

                int r1=(x1+1)%2;
                int s1=(r1+y1*2)%N;

                RR(a)+= pow(abs(evf(k,f)),2);
            }
            //=======================================================================//







            //=======================================================================//        
            if(N==2)
            {
                if((x+y)%2==0) //+- hoping
                {
                    RR(0)+= pow(abs(evf(k,f)),2);
                }
                    else
                {
                    RR(1)+= pow(abs(evf(k,f)),2);
                }
                //=======================================================================//
                //=======================================================================//        
                
            }
            //==========================================================================//
            if(N==1)
            {
                RR(0)+= pow(abs(evf(k,f)),2);
            }
            //==============================================================================//
            
        }
        RR=RR*N/(Nf*Nf);
        for (int i = 0; i < RR.size(); i++)
        {
            aaa+=RR(i);
        }
        if(aaa>0.000001)
        file<< RR.transpose() << "     "<< ef(f)<< endl;
        //cout<<"fermidos" <<RR.transpose() << "     "<< ef(f)<< endl;
    }
    //===============================================//
    
    //================================================================================/
    
    file.close();
}
//===============================================================//
//===============================================================//
#endif










 // double kx=0.0,ky=0.0;
    // double ns=0.0;
    // //======================================================================================//
    // vector<double>ef;
    // kx=-M_PI; //---------------------------------------- starting value asinement of zone
    // ky=-M_PI/sqrt(3); //-------------------------------- starting value of zone
    
    // double sx=(M_PI-2*M_PI/Nf); //------------------------ ending value of zone
    // double sy=(M_PI/sqrt(3) -2*M_PI/(Nf*sqrt(3))); //----- ending value of zone
    // double erx=M_PI/Nf/2 , ery=M_PI/(Nf*sqrt(3))/2; //----- error factor to remove error  

    // kx=-M_PI; 
    //     do
    //     {
    //         ky=-M_PI/sqrt(3);
    //         do
    //             {






    //                 ky+=2*M_PI/(sqrt(3)*Nf);
    //             }while (ky<=sy+ery);
    //         kx+=2*M_PI/Nf;
    //     } while (kx<=sx+erx);