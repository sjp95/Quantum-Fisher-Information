//||==================================================================================================||//
//||              ****** Hedder file to compute onbody tightbinding ouxilary ******                   ||//
//||                            ******* Fermionic Hamiltonian *******                                 ||//                                                           ||*/
//||==================================================================================================||//

//===========================//|| H
#ifndef HFERMI_HPP_INCLUDED  //|| E
#define HFERMI_HPP_INCLUDED  //|| A
#include "../input.hpp"         //|| D
#include <Eigen/Dense>       //|| E
using namespace std;         //|| R
using namespace Eigen;       //|| F I
//===========================//|| L E


//==============================================================//
//==============================================================//
void input::Hfermi() //| === Function to reat hamiltonian
{
    Hf = MatrixXcd :: Zero(2*Nf*Nf,2*Nf*Nf); //| === Initial Zero 
    // cout<<B<<endl;

    for (int k=0;k<Nf*Nf;k++) //| === Basis loop (site basis)
    {
      //=================================================================================================//
      int x=int(k%Nf);  //| --- x vatue of site k
      int y=k/Nf;       //| --- y value of site k
      //=============================================================================//
      //=============================================================================//
      int p=(y+1)%Nf;  //-- nearest neighbour along y
      int q=(x+p*Nf)%(Nf*Nf); //---- y axis hoping point or new site  H|k> = t*Q |q>
      //----------------------------------------------------------------------------//
      int r=(x+1)%Nf; //-- nearest neighbour along x
      int s=(r+y*Nf)%(Nf*Nf); //---- x axis hoping point or new site  H|k> = t*Q |s>
      //=============================================================================//











      if(N==16)
      {

          int y1=y%4;
          int x1=x%4;

          int a=(x1+y1*4)%N;

          int p1=(y1+1)%4;
          int q1=(x1+p1*4)%N;

          int r1=(x1+1)%4;
          int s1=(r1+y1*4)%N;

          if (k!=q)
          { 
            // if(y%2==0)
            // {
              if(y1==3)
              {
                  Hf(k,q)+=complex<double> (t2(a,q1)*B(5,6));
                  Hf(q,k)+=complex<double> (conj(t2(a,q1)*B(5,6)));


                  Hf(k+Nf*Nf,q+Nf*Nf)+=complex<double> (t2(a,q1)*B(5,6));
                  Hf(q+Nf*Nf,k+Nf*Nf)+=complex<double> (conj(t2(a,q1)*B(5,6)));

                  // cout<<k<<"  "<<q<<"  "<<a<<"  "<<q1<<endl;

              }
              else
              {
                  Hf(k,q)+=complex<double> (t1(a,q1)*B(a,q1));
                  Hf(q,k)+=complex<double> (conj(t1(a,q1)*B(a,q1)));


                  Hf(k+Nf*Nf,q+Nf*Nf)+=complex<double> (t1(a,q1)*B(a,q1));
                  Hf(q+Nf*Nf,k+Nf*Nf)+=complex<double> (conj(t1(a,q1)*B(a,q1)));

                  // cout<<k<<"  "<<q<<"  "<<a<<"  "<<q1<<endl;

              }
            // }
          

          }
          //=======================================================================//
          //=======================================================================//
          if(k!=s)
          {  

                if(x1==3)
                {
                  Hf(k,s)+=complex<double> (t2(a,s1)*B(5,6));
                  Hf(s,k)+=complex<double> (conj(t2(a,s1)*B(5,6)));


                  Hf(k+Nf*Nf,s+Nf*Nf)+=complex<double> (t2(a,s1)*B(5,6));
                  Hf(s+Nf*Nf,k+Nf*Nf)+=complex<double> (conj(t2(a,s1)*B(5,6)));

                  // cout<<k<<"  "<<s<<"  "<<a<<"  "<<s1<<endl;

                }
                else
                {

                  Hf(k,s)+=complex<double> (t1(a,s1)*B(a,s1));
                  Hf(s,k)+=complex<double> (conj(t1(a,s1)*B(a,s1)));  


                  Hf(k+Nf*Nf,s+Nf*Nf)+=complex<double> (t1(a,s1)*B(a,s1));
                  Hf(s+Nf*Nf,k+Nf*Nf)+=complex<double> (conj(t1(a,s1)*B(a,s1)));  

                  // cout<<k<<"  "<<s<<"  "<<a<<"  "<<s1<<endl;

                }


          }
      }


























      if(N==8)
      {

          int y1=y%2;
          int x1=x%4;

          int a=(x1+y1*4)%N;

          int p1=(y1+1)%2;
          int q1=(x1+p1*4)%N;

          int r1=(x1+1)%4;
          int s1=(r1+y1*4)%N;

          if (k!=q)
          { 
            // if(y%2==0)
            // {
              if(y1==1)
              {
                  Hf(k,q)+=complex<double> (t2(a,q1)*B(5,6));
                  Hf(q,k)+=complex<double> (conj(t2(a,q1)*B(5,6)));


                  Hf(k+Nf*Nf,q+Nf*Nf)+=complex<double> (t2(a,q1)*B(5,6));
                  Hf(q+Nf*Nf,k+Nf*Nf)+=complex<double> (conj(t2(a,q1)*B(5,6)));

                  // cout<<k<<"  "<<q<<"  "<<a<<"  "<<q1<<endl;

              }
              else
              {
                  Hf(k,q)+=complex<double> (t1(a,q1)*B(a,q1));
                  Hf(q,k)+=complex<double> (conj(t1(a,q1)*B(a,q1)));


                  Hf(k+Nf*Nf,q+Nf*Nf)+=complex<double> (t1(a,q1)*B(a,q1));
                  Hf(q+Nf*Nf,k+Nf*Nf)+=complex<double> (conj(t1(a,q1)*B(a,q1)));

                  // cout<<k<<"  "<<q<<"  "<<a<<"  "<<q1<<endl;

              }
            // }
          

          }
          //=======================================================================//
          //=======================================================================//
          if(k!=s)
          {  

                if(x1==3)
                {
                  Hf(k,s)+=complex<double> (t2(a,s1)*B(5,6));
                  Hf(s,k)+=complex<double> (conj(t2(a,s1)*B(5,6)));


                  Hf(k+Nf*Nf,s+Nf*Nf)+=complex<double> (t2(a,s1)*B(5,6));
                  Hf(s+Nf*Nf,k+Nf*Nf)+=complex<double> (conj(t2(a,s1)*B(5,6)));

                  // cout<<k<<"  "<<s<<"  "<<a<<"  "<<s1<<endl;

                }
                else
                {

                  Hf(k,s)+=complex<double> (t1(a,s1)*B(a,s1));
                  Hf(s,k)+=complex<double> (conj(t1(a,s1)*B(a,s1)));  


                  Hf(k+Nf*Nf,s+Nf*Nf)+=complex<double> (t1(a,s1)*B(a,s1));
                  Hf(s+Nf*Nf,k+Nf*Nf)+=complex<double> (conj(t1(a,s1)*B(a,s1)));  

                  // cout<<k<<"  "<<s<<"  "<<a<<"  "<<s1<<endl;

                }


          }
      }



      //=====================================================================================================================//


      if(N==4)
      {


          int y1=y%2;
          int x1=x%2;

          int a=(x1+y1*2)%N;

          int p1=(y1+1)%2;
          int q1=(x1+p1*2)%N;

          int r1=(x1+1)%2;
          int s1=(r1+y1*2)%N;



          if (k!=q)
          { 
            // if(y%2==0)
            // {
              if(y1==1)
              {
                  Hf(k,q)+=complex<double> (t2(a,q1)*B(a,q1));
                  Hf(q,k)+=complex<double> (conj(t2(a,q1)*B(a,q1)));


                  Hf(k+Nf*Nf,q+Nf*Nf)+=complex<double> (t2(a,q1)*B(a,q1));
                  Hf(q+Nf*Nf,k+Nf*Nf)+=complex<double> (conj(t2(a,q1)*B(a,q1)));

                  // cout<<k<<"  "<<q<<"  "<<a<<"  "<<q1<<endl;

              }
              else
              {
                  Hf(k,q)+=complex<double> (t1(a,q1)*B(a,q1));
                  Hf(q,k)+=complex<double> (conj(t1(a,q1)*B(a,q1)));


                  Hf(k+Nf*Nf,q+Nf*Nf)+=complex<double> (t1(a,q1)*B(a,q1));
                  Hf(q+Nf*Nf,k+Nf*Nf)+=complex<double> (conj(t1(a,q1)*B(a,q1)));

                  // cout<<k<<"  "<<q<<"  "<<a<<"  "<<q1<<endl;

              }
            // }
          

          }
          //=======================================================================//
          //=======================================================================//
          if(k!=s)
          {  

                if(x1==1)
                {
                  Hf(k,s)+=complex<double> (t2(a,s1)*B(a,s1));
                  Hf(s,k)+=complex<double> (conj(t2(a,s1)*B(a,s1)));


                  Hf(k+Nf*Nf,s+Nf*Nf)+=complex<double> (t2(a,s1)*B(a,s1));
                  Hf(s+Nf*Nf,k+Nf*Nf)+=complex<double> (conj(t2(a,s1)*B(a,s1)));

                  // cout<<k<<"  "<<s<<"  "<<a<<"  "<<s1<<endl;

                }
                else
                {

                  Hf(k,s)+=complex<double> (t1(a,s1)*B(a,s1));
                  Hf(s,k)+=complex<double> (conj(t1(a,s1)*B(a,s1)));  


                  Hf(k+Nf*Nf,s+Nf*Nf)+=complex<double> (t1(a,s1)*B(a,s1));
                  Hf(s+Nf*Nf,k+Nf*Nf)+=complex<double> (conj(t1(a,s1)*B(a,s1)));  

                  // cout<<k<<"  "<<s<<"  "<<a<<"  "<<s1<<endl;

                }


          }
      }



      //=====================================================================================================================//



      if(N==2)
      {
        if (k!=q)
        { 
          if((x+y)%2==0) //+- hoping
          {
            Hf(q,k)=complex<double> (conj(t(0,1)*B(0,1)));                  //| == <q|H|k> = t*Q up spin
            Hf(k,q)=complex<double> (t(0,1)*B(0,1));                        //| == <k|H|q> = conjucate(t*Q) up spin
            Hf(q+Nf*Nf,k+Nf*Nf)=complex<double> (conj(t(0,1)*B(0,1)));  //| == <q|H|k> = t*Q down spin
            Hf(k+Nf*Nf,q+Nf*Nf)=complex<double> (t(0,1)*B(0,1));        //| == <q|H|k> = conjucate(t*Q) dowm spin
          }
            else
          {
            Hf(q,k)=complex<double> (t(0,1)*B(0,1));                       //| == <q|H|k> = t*Q up spin
            Hf(k,q)=complex<double> (conj(t(0,1)*B(0,1)));                 //| == <k|H|q> = conjucate(t*Q) up spin
            Hf(q+Nf*Nf,k+Nf*Nf)=complex<double> (t(0,1)*B(0,1));       //| == <q|H|k> = t*Q down spin
            Hf(k+Nf*Nf,q+Nf*Nf)=complex<double> (conj(t(0,1)*B(0,1))); //| == <q|H|k> = conjucate(t*Q) dowm spin
          }
        }
        //=======================================================================//
        //=======================================================================//
        if(k!=s)
        {  
          if((x+y)%2==0)//+- hoping
          {            
            Hf(s,k)=complex<double> (conj(t(0,1)*B(0,1)));              //| == <s|H|k> = t*Q up spin
            Hf(k,s)=complex<double> (t(0,1)*B(0,1));                    //| == <k|H|s> = conjucate (t*Q) up spin
            Hf(s+Nf*Nf,k+Nf*Nf)=complex<double> (conj(t(0,1)*B(0,1)));  //| == <s|H|k> = t*Q down spin 
            Hf(k+Nf*Nf,s+Nf*Nf)=complex<double> (t(0,1)*B(0,1));        //| == <k|H|s> = conjugate (t*Q) down spin 
          }
          else
          {
            Hf(s,k)=complex<double> (t(0,1)*B(0,1));                  //| == <s|H|k> = t*Q up spin
            Hf(k,s)=complex<double> (conj(t(0,1)*B(0,1)));            //| == <k|H|s> = conjucate (t*Q) up spin
            Hf(s+Nf*Nf,k+Nf*Nf)=complex<double> (t(0,1)*B(0,1));      //| == <s|H|k> = t*Q down spin 
            Hf(k+Nf*Nf,s+Nf*Nf)=complex<double> (conj(t(0,1)*B(0,1)));//| == <k|H|s> = conjugate (t*Q) down spin  
          }
        }
      }
      //=====================================================================================================================//
      if(N==1)
      {
        if (k!=q)
        {           
          Hf(q,k)=complex<double> (conj(t(0,0)*B(0,0)));              //| == <q|H|k> = t*Q up spin
          Hf(k,q)=complex<double> (t(0,0)*B(0,0));                    //| == <k|H|q> = conjucate(t*Q) up spin
          Hf(q+Nf*Nf,k+Nf*Nf)=complex<double> (conj(t(0,0)*B(0,0)));  //| == <q|H|k> = t*Q down spin
          Hf(k+Nf*Nf,q+Nf*Nf)=complex<double> (t(0,0)*B(0,0));        //| == <q|H|k> = conjucate(t*Q) dowm spin
        }
        //=======================================================================//
        //=======================================================================//
        if(k!=s)
        {                    
          Hf(s,k)=complex<double> (conj(t(0,0)*B(0,0)));              //| == <s|H|k> = t*Q up spin
          Hf(k,s)=complex<double> (t(0,0)*B(0,0));                    //| == <k|H|s> = conjucate (t*Q) up spin
          Hf(s+Nf*Nf,k+Nf*Nf)=complex<double> (conj(t(0,0)*B(0,0)));  //| == <s|H|k> = t*Q down spin 
          Hf(k+Nf*Nf,s+Nf*Nf)=complex<double> (t(0,0)*B(0,0));        //| == <k|H|s> = conjugate (t*Q) down spin           
        }
      }
      //=====================================================================================================================//



    }
}
//==============================================================//
//==============================================================//
#endif


