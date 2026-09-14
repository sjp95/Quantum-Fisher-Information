#ifndef JIJJ_HPP_INCLUDED
#define JIJJ_HPP_INCLUDED
#include "../input.hpp"
#include "Hfermi.hpp"
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;


//==============================================================//
//==============================================================//
void input::xij(int i)
{
    MatrixXcd Jz =MatrixXcd :: Zero(N,N);
    MatrixXcd Jz1 =MatrixXcd :: Zero(N,N);
    MatrixXcd Jz2 =MatrixXcd :: Zero(N,N);

    double bond =0;
    double bond1 =0;
    double bond2 =0;


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








      if(N==16)
      {

          int y1=y%4;
          int x1=x%4;

          int a=(x1+y1*4)%N;

          int p1=(y1+1)%4;
          int q1=(x1+p1*4)%N;

          int r1=(x1+1)%4;
          int s1=(r1+y1*4)%N;
          //=========================================================//

          if (k!=q)
          { 


              if(y1==3)
              {
                  Jz2(a,q1)+=conj(evf(q,i))*evf(k,i)+conj(evf(q+Nf*Nf,i))*evf(k+Nf*Nf,i);
                  Jz2(q1,a)+=conj(evf(k,i))*evf(q,i)+conj(evf(k+Nf*Nf,i))*evf(q+Nf*Nf,i);
                  if(x1==0)
                  bond1+=1;

                  // cout<<a<<"  "<<q1<<"  "<<k<<"  "<<q<<endl;
  
              }
              else
              {
                  Jz1(a,q1)+=conj(evf(q,i))*evf(k,i)+conj(evf(q+Nf*Nf,i))*evf(k+Nf*Nf,i);
                  Jz1(q1,a)+=conj(evf(k,i))*evf(q,i)+conj(evf(k+Nf*Nf,i))*evf(q+Nf*Nf,i);
                  // cout<<a<<"  "<<q1<<"  "<<k<<"  "<<q<<endl;
                  // cout<<k<<"  "<<q<<"   "<<conj(evf(q,i))*evf(k,i)<<"   "<<conj(evf(q+Nf*Nf,i))*evf(k+Nf*Nf,i)<<endl;

              }

          }
          //=======================================================================//
          //=======================================================================//        
          if(k!=s)
          {  
              if(x1==3)
              {
                  Jz2(a,s1)+=conj(evf(s,i))*evf(k,i)+conj(evf(s+Nf*Nf,i))*evf(k+Nf*Nf,i);
                  Jz2(s1,a)+=conj(evf(k,i))*evf(s,i)+conj(evf(k+Nf*Nf,i))*evf(s+Nf*Nf,i);
                  if(y1==0)
                  bond2+=1;
                  // cout<<a<<"  "<<s1<<"  "<<k<<"  "<<s<<endl;
                  // cout<<k<<"  "<<s<<"   "<<conj(evf(s,i))*evf(k,i)<<"   "<<conj(evf(s+Nf*Nf,i))*evf(k+Nf*Nf,i)<<endl;

              }
              else
              {

                  Jz1(a,s1)+=conj(evf(s,i))*evf(k,i)+conj(evf(s+Nf*Nf,i))*evf(k+Nf*Nf,i);
                  Jz1(s1,a)+=conj(evf(k,i))*evf(s,i)+conj(evf(k+Nf*Nf,i))*evf(s+Nf*Nf,i);
                  // cout<<a<<"  "<<s1<<"  "<<k<<"  "<<s<<endl;
                  // cout<<i<<"   "<<k<<"  "<<s<<"   "<<conj(evf(s,i))*evf(k,i)<<"   "<<conj(evf(s+Nf*Nf,i))*evf(k+Nf*Nf,i)<<endl;

              }
          }
      

      }


       //=======================================================================//
          //=======================================================================//        














      if(N==8)
      {

          int y1=y%2;
          int x1=x%4;

          int a=(x1+y1*4)%N;

          int p1=(y1+1)%2;
          int q1=(x1+p1*4)%N;

          int r1=(x1+1)%4;
          int s1=(r1+y1*4)%N;
          //=========================================================//

          if (k!=q)
          { 


              if(y1==1)
              {
                  Jz2(a,q1)+=conj(evf(q,i))*evf(k,i)+conj(evf(q+Nf*Nf,i))*evf(k+Nf*Nf,i);
                  Jz2(q1,a)+=conj(evf(k,i))*evf(q,i)+conj(evf(k+Nf*Nf,i))*evf(q+Nf*Nf,i);
                  if(x1==0)
                  bond1+=1;

                  // cout<<a<<"  "<<q1<<"  "<<k<<"  "<<q<<endl;



  
              }
              else
              {
                  Jz1(a,q1)+=conj(evf(q,i))*evf(k,i)+conj(evf(q+Nf*Nf,i))*evf(k+Nf*Nf,i);
                  Jz1(q1,a)+=conj(evf(k,i))*evf(q,i)+conj(evf(k+Nf*Nf,i))*evf(q+Nf*Nf,i);
                  // cout<<a<<"  "<<q1<<"  "<<k<<"  "<<q<<endl;
                  // cout<<k<<"  "<<q<<"   "<<conj(evf(q,i))*evf(k,i)<<"   "<<conj(evf(q+Nf*Nf,i))*evf(k+Nf*Nf,i)<<endl;

              }

          }
          //=======================================================================//
          //=======================================================================//        
          if(k!=s)
          {  
              if(x1==3)
              {
                  Jz2(a,s1)+=conj(evf(s,i))*evf(k,i)+conj(evf(s+Nf*Nf,i))*evf(k+Nf*Nf,i);
                  Jz2(s1,a)+=conj(evf(k,i))*evf(s,i)+conj(evf(k+Nf*Nf,i))*evf(s+Nf*Nf,i);
                  if(y1==0)
                  bond2+=1;
                  // cout<<a<<"  "<<s1<<"  "<<k<<"  "<<s<<endl;
                  // cout<<k<<"  "<<s<<"   "<<conj(evf(s,i))*evf(k,i)<<"   "<<conj(evf(s+Nf*Nf,i))*evf(k+Nf*Nf,i)<<endl;

              }
              else
              {

                  Jz1(a,s1)+=conj(evf(s,i))*evf(k,i)+conj(evf(s+Nf*Nf,i))*evf(k+Nf*Nf,i);
                  Jz1(s1,a)+=conj(evf(k,i))*evf(s,i)+conj(evf(k+Nf*Nf,i))*evf(s+Nf*Nf,i);
                  // cout<<a<<"  "<<s1<<"  "<<k<<"  "<<s<<endl;
                  // cout<<i<<"   "<<k<<"  "<<s<<"   "<<conj(evf(s,i))*evf(k,i)<<"   "<<conj(evf(s+Nf*Nf,i))*evf(k+Nf*Nf,i)<<endl;

              }
          }
      

      }


       //=======================================================================//
          //=======================================================================//        







      if(N==4)
      {

          int y1=y%2;
          int x1=x%2;

          int a=(x1+y1*2)%N;

          int p1=(y1+1)%2;
          int q1=(x1+p1*2)%N;

          int r1=(x1+1)%2;
          int s1=(r1+y1*2)%N;
          //=========================================================//

          if (k!=q)
          { 


              if(y1==1)
              {
                  Jz2(a,q1)+=conj(evf(q,i))*evf(k,i)+conj(evf(q+Nf*Nf,i))*evf(k+Nf*Nf,i);
                  Jz2(q1,a)+=conj(evf(k,i))*evf(q,i)+conj(evf(k+Nf*Nf,i))*evf(q+Nf*Nf,i);
                  if(x1==0)
                  bond1+=1;

                  // cout<<a<<"  "<<q1<<"  "<<k<<"  "<<q<<endl;
                  // cout<<k<<"  "<<q<<"   "<<conj(evf(q,i))*evf(k,i)<<"   "<<conj(evf(q+Nf*Nf,i))*evf(k+Nf*Nf,i)<<endl;
                  // cout<<k<<"  "<<q<<"   "<<Jz2(a,q1)<<"   "<<Jz2(q1,a)<<endl;



  
              }
              else
              {
                  Jz1(a,q1)+=conj(evf(q,i))*evf(k,i)+conj(evf(q+Nf*Nf,i))*evf(k+Nf*Nf,i);
                  Jz1(q1,a)+=conj(evf(k,i))*evf(q,i)+conj(evf(k+Nf*Nf,i))*evf(q+Nf*Nf,i);
                  // cout<<a<<"  "<<q1<<"  "<<k<<"  "<<q<<endl;
                  // cout<<k<<"  "<<q<<"   "<<conj(evf(q,i))*evf(k,i)<<"   "<<conj(evf(q+Nf*Nf,i))*evf(k+Nf*Nf,i)<<endl;

              }

          }
          //=======================================================================//
          //=======================================================================//        
          if(k!=s)
          {  
              if(x1==1)
              {
                  Jz2(a,s1)+=conj(evf(s,i))*evf(k,i)+conj(evf(s+Nf*Nf,i))*evf(k+Nf*Nf,i);
                  Jz2(s1,a)+=conj(evf(k,i))*evf(s,i)+conj(evf(k+Nf*Nf,i))*evf(s+Nf*Nf,i);
                  if(y1==0)
                  bond2+=1;
                  // cout<<a<<"  "<<s1<<"  "<<k<<"  "<<s<<endl;
                  // cout<<k<<"  "<<s<<"   "<<conj(evf(s,i))*evf(k,i)<<"   "<<conj(evf(s+Nf*Nf,i))*evf(k+Nf*Nf,i)<<endl;

              }
              else
              {

                  Jz1(a,s1)+=conj(evf(s,i))*evf(k,i)+conj(evf(s+Nf*Nf,i))*evf(k+Nf*Nf,i);
                  Jz1(s1,a)+=conj(evf(k,i))*evf(s,i)+conj(evf(k+Nf*Nf,i))*evf(s+Nf*Nf,i);
                  // cout<<a<<"  "<<s1<<"  "<<k<<"  "<<s<<endl;
                  // cout<<i<<"   "<<k<<"  "<<s<<"   "<<conj(evf(s,i))*evf(k,i)<<"   "<<conj(evf(s+Nf*Nf,i))*evf(k+Nf*Nf,i)<<endl;

              }
          }
      

      }


       //=======================================================================//
          //=======================================================================//        



      if(N==2)
      {
          if (k!=q)
          { 
            if((x+y)%2==0) //+- hoping
            {
              Jz(1,0)+=conj(evf(q,i))*evf(k,i)+conj(evf(q+Nf*Nf,i))*evf(k+Nf*Nf,i);
              bond+=1;
            }
              else
            {
              Jz(1,0)+=conj(evf(k,i))*evf(q,i)+conj(evf(k+Nf*Nf,i))*evf(q+Nf*Nf,i);
              bond+=1;
            }
          }
          //=======================================================================//
          //=======================================================================//        
          if(k!=s)
          {  
            if((x+y)%2==0)//+- hoping
            {            
              Jz(1,0)+=conj(evf(s,i))*evf(k,i)+conj(evf(s+Nf*Nf,i))*evf(k+Nf*Nf,i);
              bond+=1;
            }
            else
            {
              Jz(1,0)+=conj(evf(k,i))*evf(s,i)+conj(evf(k+Nf*Nf,i))*evf(s+Nf*Nf,i);
              bond+=1;
            }
          }
      }
      //==========================================================================//
      if(N==1)
      {
          if (k!=q)
          {           
            Jz(0,0)+=conj(evf(q,i))*evf(k,i)+conj(evf(q+Nf*Nf,i))*evf(k+Nf*Nf,i);
            bond+=1;
          }
          //=======================================================================//
          //=======================================================================//
          if(k!=s)
          {           
            Jz(0,0)+=conj(evf(s,i))*evf(k,i)+conj(evf(s+Nf*Nf,i))*evf(k+Nf*Nf,i);
            bond+=1;
          }
      }
      //==============================================================================//
      
    }
    if(N==2)
    {
     Jz(0,1)=conj(Jz(1,0));
    }
    X+= Jz*Fermi(ef(i),muf)/bond;


    Jz1=Jz1/bond1;
    Jz2=Jz2/bond1;




    // Jz1=4*Jz1/(Nf*Nf);
    // Jz2=4*Jz2/(Nf*Nf);




    X2+= Jz2*Fermi(ef(i),muf);
    X1+= Jz1*Fermi(ef(i),muf);







}


void input::Xij()
{
    X=MatrixXcd :: Zero(N,N);
    X1=MatrixXcd :: Zero(N,N);
    X2=MatrixXcd :: Zero(N,N);
    for (int i=0; i< Hf.col(0).size();i++)
    {
      xij(i);
    }
    // cout<<X2.real()<<endl;
}
//==============================================================//
//==============================================================//
#endif




//=======================================================================//
//=======================================================================//
    //   int u=(r+p*Nf)%(Nf*Nf); //-- next nearest neighbour along 
    //   if(k!=u)
    //   { 
    //     if((x+y)%2==0)
    //     { 
    //       Hf(u,k)=complex<double> t(0,1)*Q(0,2);
    //       Hf(k,u)=complex<double> conj(t(0,1)*Q(0,2));
    //       Hf(u+Nf*Nf,k+Nf*Nf)=complex<double> t(0,1)*Q(0+N,2+N);  
    //       Hf(k+Nf*Nf,u+Nf*Nf)=complex<double> conj(t(0,1)*Q(0+N,2+N));
    //     }
    //     else
    //     {
    //       Hf(u,k)=complex<double> t(0,1)*Q(0,2);
    //       Hf(k,u)=complex<double> conj(t(0,1)*Q(0,2));
    //       Hf(u+Nf*Nf,k+Nf*Nf)=complex<double> t(0,1)*Q(0+N,2+N);  
    //       Hf(k+Nf*Nf,u+Nf*Nf)=complex<double> conj(t(0,1)*Q(0+N,2+N));
    //     } 
    //   }

    //   if(q!=s)// q=y+1, s= x+1
    //   { 
    //     if((x+y)%2!=0)
    //     { 
    //       Hf(q,s)=complex<double> t(0,1)*Q(1,3);
    //       Hf(s,q,)=complex<double> conj(t(0,1)*Q(1,3));
    //       Hf(q+Nf*Nf,s+Nf*Nf)=complex<double> t(0,1)*Q(1+N,3+N);  
    //       Hf(s+Nf*Nf,q+Nf*Nf)=complex<double> conj(t(0,1)*Q(1+N,3+N));
    //     }
    //     else
    //     {
    //        Hf(q,s)=complex<double> t(0,1)*Q(1,3);
    //       Hf(s,q,)=complex<double> conj(t(0,1)*Q(1,3));
    //       Hf(q+Nf*Nf,s+Nf*Nf)=complex<double> t(0,1)*Q(1+N,3+N);  
    //       Hf(s+Nf*Nf,q+Nf*Nf)=complex<double> conj(t(0,1)*Q(1+N,3+N));
    //     } 
    //   }
      
    //   // if((x+y)%2==0)
    //   // {
    //   //   Hf(k,k)=cd(dlt)*0.5;
    //   //   Hf(k+N*N,k+N*N)=cd(dlt)*0.5;
    //   // }
    //   // else
    //   // {
    //   //   Hf(k,k)=-cd(dlt)*0.5;
    //   //   Hf(k+N*N,k+N*N)=-cd(dlt)*0.5;
    //   // }
    // //   Hf(k,k)=-h*0.5;
    // //   Hf(k+N*N,k+N*N)=h*0.5;
    // }
