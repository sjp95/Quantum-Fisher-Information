#ifndef SPINCOVGERGE_HPP_INCLUDED
#define SPINCOVGERGE_HPP_INCLUDED
#include "../input.hpp"
#include "Hspin.hpp"
#include "o.hpp"
#include "ocupation.hpp"
#include <Eigen/Dense>
#include "../matworks.hpp"
#define roun_d(x,d) ((floor(((x)*pow(10,d))+.5))/pow(10,d)) //======= round of  number
using namespace std;
using namespace Eigen;


char title1[100];
char title2[100];


void input::O_converge()
{
    int count = 0;
    VectorXcd Os0 = VectorXcd :: Zero(N);
    double ppr=0.0;
    int t1=int(U);
    
      do
        {

            Os0=O;
            Hspin();  
            evs=Eigenvectors(H);  
            OO();
            ppr=0.0;
            for(int i=0;i<Os0.size();i++)
            {
                ppr+=abs(O(i)-Os0(i));
            }
            ppr=double(ppr/Os0.size());
            count++;
            if(count>1000)
            {
                cout<< mus <<"    "<<ppr<< " O_convo" <<endl;
                break;
            }
            sprintf(title1,"CURRENT_STATUS1_%d.dat",t1); //| ==== Name | I
            std::ofstream file1(title1, std::fstream::out | std::fstream::app); 
            file1<< count <<"    "<<abs(ppr) << endl;
    
        }
        while(abs(ppr)>0.001);
}

//=================================================================//
void input::Spin_converge()
{
    //================Bisection start=================
    //=================================================
    double ml=-U-10; //----------- minimum energy
    double mu=U+10;//--- maximum energy
    mus=(ml+mu)*0.5;
    int t2=int(U);

    B= MatrixXcd :: Identity(N,N);
    B= B*(1.0/(sqrt(2)+U));
    O_converge();
    ocupation();


    int count=0;
    do
    {
        mus=(ml+mu)*0.5;
        B= MatrixXcd :: Identity(N,N);
        B= B*(1.0/(sqrt(2)+U));
        O_converge();
        ocupation();
        Stotal=Stotal/double(N);
        if(Stotal>Spin_fill)     //|
            mu=mus;              //|
        if (Stotal<Spin_fill)    //|===== Bisection
            ml=mus;              //|
        if(count>5000)           //|
        {
            cout<< "Bisection falied- "<<U<< " _ " << Stotal <<endl;
            break;
        }
        count++;
        sprintf(title2,"CURRENT_STATUS2_%d.dat",t2); //| ==== Name | I
        std::ofstream file2(title2, std::fstream::out | std::fstream::app); 
        file2<< count <<"  "<< Stotal<<"  "<< Spin_fill<<"    "<<abs(Stotal-Spin_fill) << endl;

        //Stotal=roun_d(Stotal,3);
        //cout<< count << "   "<< O.transpose()<< " _ " << Stotal << endl;
    }while(abs(Stotal-Spin_fill)>0.0001);
    //================Bisection end=================//
    //==============================================//
     //Q= MatrixXcd :: Identity(2*N,2*N);
     //O_converge();
}

//========================================================================================//
//========================================================================================//

#endif
