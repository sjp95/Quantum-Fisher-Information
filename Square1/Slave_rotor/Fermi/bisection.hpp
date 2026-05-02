#ifndef BISECTION_HPP_INCLUDED
#define BISECTION_HPP_INCLUDED
#include "../input.hpp"
#include "Hfermi.hpp"
#include "jijj.hpp"
#include <Eigen/Dense>
//#define roun_d(x,d) ((floor(((x)*pow(10,d))+.5))/pow(10,d)) //======= round of  number
using namespace std;
using namespace Eigen;

void input::Bisection()
{
    //================Bisection start=================
    //=================================================
    double ml=ef(0); //----------- minimum energy
    double mu=ef(ef.size()-1);//--- maximum energy
    double mm=(mu+ml)/2+0.0;//---------average energy
    double fm=0;
    for (int i=0;i<ef.size();i++)
        fm+=1/(exp((ef(i)-mm)/T)+1); //-------- feeling at avg energy
        fm=fm/double(Nf*Nf);
    int count=0;
    while(abs(fm-filling)>0.000001)
    {
        mm=(mu+ml)/2;
        fm=0;
        for (int i=0;i<ef.size();i++)
            fm+=1/(exp((ef(i)-mm)/T)+1); //----- feeling count
            fm=fm/double(Nf*Nf);
        if(fm>filling)     //|
            mu=mm;         //|
        if (fm<filling)    //|===== Bisection
            ml=mm;         //|
        if(count>10000)    //|
        {
            cout<< "Bisection falied: " <<endl;
            break;
        }
        count++;
    }
    //================Bisection end=================//
    //==============================================//
    muf=mm;
}

//========================================================================================//
//========================================================================================//

#endif
