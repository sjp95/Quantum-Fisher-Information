#ifndef DOS_HPP_INCLUDED
#define DOS_HPP_INCLUDED
#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include <math.h>
#include <fstream>
#include <vector>
#include "../input.hpp"
#include "spindos.hpp"
#include "fermidos.hpp"


using namespace std;
using namespace Eigen;
char title_dos[100];

// double delta(double x0, double x)
// {
//     double L=0.15;
//     double et=4.0/(Nf*Nf);
//     double D= (L/(2*M_PI))/(et+(x-x0)*(x-x0)); //-------- Delta function
//     return D;
// }

typedef std::pair<Eigen::MatrixXcd, Eigen::VectorXd> Spectrum;
//===============================================================//
//===============================================================//
void input::DoS()
{
    //------------------------------------------------------------//
    //------------------------------------------------------------//
     int tut=int(U*1000);                         //| ==== File | F
    //  int tut0=int(A*1000);
     sprintf(title_dos,"../Data/DOS_%d.dat",tut); //| ==== File | F
     std::ofstream file;                         //| ==== Name | I
     file.open(title_dos);                       //| ==== File | L
    
    //------------------------------------------------------------//
    
    
    //------------------------------------------------------------//
    //------------------------------------------------------------//
    Spin_dos();  //|-------- Calculation of spin part
    Fermi_dos(); //|-------- Calculation of fermi part
    //.....................................................//
    double nnnppp=double(Nf*Nf);
    double w = -U*0.5-3; //-------- initial frequency;
    do
    {
        ifstream files;           //|------ FOR ROTOR FILE
        files.open(title_s_d);   //|------ ROTOR FILE READ
        VectorXcd mon = VectorXcd :: Zero(N); 
        VectorXd Fr = VectorXd :: Zero(N);

        VectorXcd ds = VectorXcd :: Zero(N);
        double PZm=0, PZn=0, ef=0,em=0,en=0; 






        if(N==16)
        {
            while(files>>mon(0)>>mon(1)>>mon(2)>>mon(3)>>mon(4)>>mon(5)>>mon(6)>>mon(7)>>mon(8)>>mon(9)>>mon(10)>>mon(11)>>mon(12)>>mon(13)>>mon(14)>>mon(15)>>PZm>>PZn>>em>>en)
            {
                ifstream filef;           //|------ FOR SPINON FILE
                filef.open(title_f_d);   //|------ SPINON FILE READ
                ef=0;
                Fr = VectorXd :: Zero(N);
                while(filef>>Fr(0)>>Fr(1)>>Fr(2)>>Fr(3)>>Fr(4)>>Fr(5)>>Fr(6)>>Fr(7)>>Fr(8)>>Fr(9)>>Fr(10)>>Fr(11)>>Fr(12)>>Fr(13)>>Fr(14)>>Fr(15)>>ef)
                {
                    ds(0)+=Fr(0)*pow(abs(mon(0)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(1)+=Fr(1)*pow(abs(mon(1)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(2)+=Fr(2)*pow(abs(mon(2)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(3)+=Fr(3)*pow(abs(mon(3)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(4)+=Fr(4)*pow(abs(mon(4)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(5)+=Fr(5)*pow(abs(mon(5)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(6)+=Fr(6)*pow(abs(mon(6)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(7)+=Fr(7)*pow(abs(mon(7)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    


                    ds(8)+=Fr(8)*pow(abs(mon(8)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(9)+=Fr(9)*pow(abs(mon(9)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(10)+=Fr(10)*pow(abs(mon(10)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(11)+=Fr(11)*pow(abs(mon(11)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(12)+=Fr(12)*pow(abs(mon(12)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(13)+=Fr(13)*pow(abs(mon(13)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(14)+=Fr(14)*pow(abs(mon(14)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(15)+=Fr(15)*pow(abs(mon(15)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);


                }
            }
            double ddd=0;
            for(int p=0;p<ds.size();p++)
            {
             ddd+=ds(p).real();
             ds(p)=ds(p)/nnnppp;
            } 
            ddd=ddd/(Nf*Nf);
            file<< w << "   " << ddd << "   " << ds.transpose().real() <<endl;  

        }














        if(N==8)
        {
            while(files>>mon(0)>>mon(1)>>mon(2)>>mon(3)>>mon(4)>>mon(5)>>mon(6)>>mon(7)>>PZm>>PZn>>em>>en)
            {
                ifstream filef;           //|------ FOR SPINON FILE
                filef.open(title_f_d);   //|------ SPINON FILE READ
                ef=0;
                Fr = VectorXd :: Zero(N);
                while(filef>>Fr(0)>>Fr(1)>>Fr(2)>>Fr(3)>>Fr(4)>>Fr(5)>>Fr(6)>>Fr(7)>>ef)
                {
                    ds(0)+=Fr(0)*pow(abs(mon(0)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(1)+=Fr(1)*pow(abs(mon(1)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(2)+=Fr(2)*pow(abs(mon(2)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(3)+=Fr(3)*pow(abs(mon(3)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(4)+=Fr(4)*pow(abs(mon(4)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(5)+=Fr(5)*pow(abs(mon(5)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(6)+=Fr(6)*pow(abs(mon(6)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(7)+=Fr(7)*pow(abs(mon(7)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    
                }
            }
            double ddd=0;
            for(int p=0;p<ds.size();p++)
            {
             ddd+=ds(p).real();
             ds(p)=ds(p)/nnnppp;
            } 
            ddd=ddd/(Nf*Nf);
            file<< w << "   " << ddd << "   " << ds.transpose().real() <<endl;  

        }




        if(N==4)
        {
            while(files>>mon(0)>>mon(1)>>mon(2)>>mon(3)>>PZm>>PZn>>em>>en)
            {
                ifstream filef;           //|------ FOR SPINON FILE
                filef.open(title_f_d);   //|------ SPINON FILE READ
                ef=0;
                Fr = VectorXd :: Zero(N);
                while(filef>>Fr(0)>>Fr(1)>>Fr(2)>>Fr(3)>>ef)
                {
                    ds(0)+=Fr(0)*pow(abs(mon(0)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(1)+=Fr(1)*pow(abs(mon(1)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(2)+=Fr(2)*pow(abs(mon(2)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(3)+=Fr(3)*pow(abs(mon(3)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    
                }
            }
            double ddd=0;
            for(int p=0;p<ds.size();p++)
            {
             ddd+=ds(p).real();
             ds(p)=ds(p)/nnnppp;
            } 
            ddd=ddd/(Nf*Nf);
            file<< w << "   " << ddd << "   " << ds.transpose().real() <<endl;  

        }





        if(N==2)
        {
            while(files>>mon(0)>>mon(1)>>PZm>>PZn>>em>>en)
            {
                ifstream filef;           //|------ FOR SPINON FILE
                filef.open(title_f_d);   //|------ SPINON FILE READ
                ef=0;
                Fr = VectorXd :: Zero(N);
                while(filef>>Fr(0)>>Fr(1)>>ef)
                {
                    ds(0)+=Fr(0)*pow(abs(mon(0)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                    ds(1)+=Fr(1)*pow(abs(mon(1)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);                    
                }
            }
            double ddd=0;
            for(int p=0;p<ds.size();p++)
            {
             ddd+=ds(p).real();
             ds(p)=ds(p)/nnnppp;
            } 
            ddd=ddd/(Nf*Nf);
            file<< w << "   " << ddd << "   " << ds.transpose().real() <<endl;  

        }
        if(N==1)
        {
            while(files>>mon(0)>>PZm>>PZn>>em>>en)
            {
                ifstream filef;           //|------ FOR SPINON FILE
                filef.open(title_f_d);   //|------ SPINON FILE READ
                ef=0;
                Fr = VectorXd :: Zero(N);
                while(filef>>Fr(0)>>ef)
                {
                    ds(0)+=Fr(0)*pow(abs(mon(0)),2)*(Fermi(ef,muf)*PZn+(1-Fermi(ef,muf))*PZm)*delta(w,ef-muf+em-en);
                }
            }
            double ddd=0;
            for(int p=0;p<ds.size();p++)
            {
             ddd+=ds(p).real();
             ds(p)=ds(p)/nnnppp;
            } 
            ddd=ddd/(Nf*Nf);
            file<< w << "   " << ddd << "   " << ds.transpose().real() <<endl;  

        }

        w+=0.01;//.0/(2.0*Nf*Nf);
    } while (w<U*0.5+3);
    
    
    file.close(); //------------------ file closing
    remove(title_s_d); //------------ delete rotor file
    remove(title_f_d); //----------- delete spinon file
}
//===============================================================//
//===============================================================//
#endif
