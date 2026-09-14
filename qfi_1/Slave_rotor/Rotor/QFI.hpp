#ifndef QFI_HPP_INCLUDED
#define QFI_HPP_INCLUDED
#include "../input.hpp"
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <complex>
#include <omp.h>
using namespace std;
using namespace Eigen;


//================================================//
            //=== QFI ===//
//================================================//


void input::QFI(double q)
{ 
    std::vector<complex<double>> SZZ(le, std::complex<double>(0.0, 0.0));
    complex <double> II=complex <double>(0,1.0);

    omp_set_num_threads(4);
    #pragma omp parallel for
    for ( int i =0; i<le; i++)
    {
        int n=i;
        for ( int j =0; j<N; j++)
        {
            double reminder = n%Ls;
            n=n/Ls;
            double nt = double(reminder-1.0/2.0);
            SZZ[i]+= nt*exp(II*q*double(j));
        }
       // cout<< i<<" "<<SZZ[i].real() <<endl;
    }

    double QFI_sum_re = 0.0;
    double QFI_sum_im = 0.0;

    #pragma omp parallel for collapse(2) reduction(+:QFI_sum_re, QFI_sum_im)
    for (int l1 = 0; l1 < le; l1++)
    {
        for (int l2 = 0; l2 < le; l2++)
        {
            complex<double> element_QFI(0.0,0.0);
            for (int i = 0; i < le; i++)
            {
                element_QFI += SZZ[i] * conj(evs(i,l1)) * evs(i,l2);
            }

            complex<double> contrib =
                2.0 * pow(exp(-(es(l1)-es(0))/T) - exp(-(es(l2)-es(0))/T), 2.0)
                / (exp(-(es(l1)-es(0))/T) + exp(-(es(l2)-es(0))/T))
                * pow(abs(element_QFI), 2.0);

            QFI_sum_re += contrib.real();
            QFI_sum_im += contrib.imag();
        }
    }

    complex<double> QFI_sum(QFI_sum_re, QFI_sum_im);
    cout << "QFI_sum: " << QFI_sum/double(N) << endl;
    


    
}

#endif
