#ifndef SLAVE_ROTOR_HPP_INCLUDED
#define SLAVE_ROTOR_HPP_INCLUDED
#include <math.h>
#include <complex>
#include "input.hpp"
#include "matworks.hpp"
#include "Operators/Hspin.hpp"
#include <Eigen/Dense>
#include <highfive/H5File.hpp>
#include <highfive/eigen.hpp>
#include "creat_directory.hpp"
#include "Operators/sx.hpp"
#include "Operators/sz.hpp"
#include "Operators/QFI.hpp"
#include "values.hpp"
#include <cuda_runtime.h>
#include <cusolverDn.h>

// #include "Operators/GPU.hpp"

using namespace std;
using namespace Eigen;
char title3[100];
char title4[100];
char title5[100];

void input::besis(int NN1,double Jjx,double Jjz,double hx0, double hy0, double hz0)
{
    N=NN1;
    J1=Jjx;
    J2=Jjz;
    hx= hx0;
    hy= hy0;
    hz= hz0;

    Ls=2*lz+1;
    ls = 0;                 // Besis count start
    le = pow(Ls,N);         // Besis count End
    //==============================//
    H= MatrixXcd :: Zero(le,le);
    Jx= MatrixXcd :: Zero(N,N);
    Jz= MatrixXcd :: Zero(N,N); 
    //t= MatrixXcd::Zero(N,N);
    evs = MatrixXcd :: Zero(le,le);
    es = VectorXd :: Zero(le);
    //==============================//
     // Precompute powers
    Lspow.resize(N);
    Lspow[0] = 1;

    for(int i=1;i<N;i++)
        Lspow[i] = Lspow[i-1] * Ls;
    cout<< "==============================="<<endl;
    cout<< "Data Alocation Done"<<endl;
    cout<< "==============================="<<endl;
    
    Values();
}



//==========================================================================================//
void input::output()
{
//     cout<<"\nU: "<<U<<"\nPhi: "<< O.transpose()<<"\nRotor occupation: "<< Stotal<<endl;
}

//==========================================================================================//    
void input::basis_print()
{
    cout<< "==============================="<<endl;
    cout<< "Basis Generated"<<endl;
    cout<< "==============================="<<endl;
    for (int i=0; i<le; i++)
    {
        int n=i;
        VectorXi nnt = VectorXi :: Zero(N);
        for ( int j =0; j<N; j++)
        {
            double reminder = n%Ls;
            n=n/Ls;
            nnt(j)=int(reminder);
        }
        cout<< i <<"  |"<< nnt.transpose()<<">"<<endl;
        
    }
    cout<< "==============================="<<endl;
}

//==========================================================================================//

void input::ground_state(bool clean)
{
    createDirectory("../Data/Eigen");

    std::string outfile =
    "../Data/Eigen/GroundState_" +
    std::to_string(N) + "_" +
    fmt(J1) + "_" +
    fmt(J2) + "_" +
    fmt(hz) + ".h5";

    bool loaded_from_file = false;

    cout << "===============================\n";
    cout << "Ground State Mode\n";
    cout << "hz: " << hz << "  Jx: " << J1 << "  Jz: " << J2 << "\n";
    cout << "Output file: " << outfile << "\n";
    cout << "===============================\n";

    //----------------------------------------------------------
    // Try loading cached ground state
    //----------------------------------------------------------
    if (std::filesystem::exists(outfile))
    {
        try
        {
            HighFive::File file(outfile, HighFive::File::ReadOnly);
            file.getDataSet("gs_energy").read(gs_energy);
            file.getDataSet("gs_vector").read(gs_vector);
            MatrixXcd().swap(H);
            loaded_from_file = true;
            cout << "Ground state loaded from HDF5.\n";
        }
        catch (const std::exception& e)
        {
            cerr << "HDF5 Read Error: " << e.what() << " — recomputing.\n";
        }
    }

    //----------------------------------------------------------
    // Compute via Eigen SelfAdjointEigenSolver
    //----------------------------------------------------------
    if (!loaded_from_file)
    {
        Hspin();
        cout << "H generated. H(0,0): " << H(0,0) << "\n";
        cout << "Diagonalizing ground state on GPU (cuSOLVER)...\n";

        MatrixXd H_real = H.real();
        MatrixXcd().swap(H);

        int n = H_real.rows();

        //----------------------------------------------------------
        // Host storage (column-major, cuSOLVER expects this)
        //----------------------------------------------------------
        std::vector<double> h_A(n * n);
        // Eigen is column-major by default so this is a direct copy
        memcpy(h_A.data(), H_real.data(), sizeof(double) * n * n);

        std::vector<double> h_W(n);        // all eigenvalues (we read only [0])
        std::vector<double> h_V(n * n);    // eigenvectors (we read only col 0)

        //----------------------------------------------------------
        // GPU allocations
        //----------------------------------------------------------
        double* d_A   = nullptr;
        double* d_W   = nullptr;
        int*    d_info = nullptr;

        cudaMalloc(&d_A,    sizeof(double) * n * n);
        cudaMalloc(&d_W,    sizeof(double) * n);
        cudaMalloc(&d_info, sizeof(int));

        cudaMemcpy(d_A, h_A.data(), sizeof(double) * n * n, cudaMemcpyHostToDevice);

        //----------------------------------------------------------
        // cuSOLVER setup
        //----------------------------------------------------------
        cusolverDnHandle_t handle;
        cusolverDnCreate(&handle);

        // Ask cuSOLVER how much workspace it needs
        int lwork = 0;
        cusolverDnDsyevd_bufferSize(
            handle,
            CUSOLVER_EIG_MODE_VECTOR,   // compute eigenvectors
            CUBLAS_FILL_MODE_LOWER,     // use lower triangle of H
            n, d_A, n, d_W, &lwork);

        double* d_work = nullptr;
        cudaMalloc(&d_work, sizeof(double) * lwork);

        //----------------------------------------------------------
        // Diagonalize — full spectrum on GPU, cheapest cuSOLVER path
        // cusolverDnDsyevdx (subset) is available but requires
        // cuSOLVER >= 11.0 and more setup; full diag is simpler and
        // GPU-side cost is still O(n^3) either way for dense matrices.
        //----------------------------------------------------------
        cusolverDnDsyevd(
            handle,
            CUSOLVER_EIG_MODE_VECTOR,
            CUBLAS_FILL_MODE_LOWER,
            n, d_A, n, d_W,
            d_work, lwork, d_info);

        // Check for convergence
        int info_host = 0;
        cudaMemcpy(&info_host, d_info, sizeof(int), cudaMemcpyDeviceToHost);
        if (info_host != 0)
            cerr << "cuSOLVER dsyevd failed, info = " << info_host << "\n";

        // Copy back only what we need: eigenvalue[0] and eigenvector col 0
        cudaMemcpy(h_W.data(), d_W, sizeof(double) * n,     cudaMemcpyDeviceToHost);
        cudaMemcpy(h_V.data(), d_A, sizeof(double) * n * n, cudaMemcpyDeviceToHost);
        // Note: after dsyevd, d_A is overwritten with eigenvectors

        //----------------------------------------------------------
        // Extract ground state
        //----------------------------------------------------------
        gs_energy = h_W[0];
        gs_vector.resize(n);
        for (int i = 0; i < n; i++)
            gs_vector(i) = complex<double>(h_V[i], 0.0);  // col 0 = first n elements

        cout << "Ground state energy: " << gs_energy << "\n";
        cout << "===============================\n";

        //----------------------------------------------------------
        // Cleanup
        //----------------------------------------------------------
        cudaFree(d_A);
        cudaFree(d_W);
        cudaFree(d_work);
        cudaFree(d_info);
        cusolverDnDestroy(handle);

        //----------------------------------------------------------
        // Cache to HDF5
        //----------------------------------------------------------
        if (!clean)
        {
            try
            {
                HighFive::File file(outfile, HighFive::File::Overwrite);
                file.createDataSet("gs_energy", gs_energy);
                file.createDataSet("gs_vector", gs_vector);
                cout << "Ground state saved to HDF5.\n";
            }
            catch (const std::exception& e)
            {
                cerr << "HDF5 Write Error: " << e.what() << "\n";
            }
        }
    }
}


//==========================================================================================//

void input::mu_phi(bool clean)
{
    bool tut = createDirectory("../Data/Eigen");

    std::string outfile =
    "../Data/Eigen/EigenSpectrum_" +
    std::to_string(N) + "_" +
    // std::to_string(static_cast<int>(std::round(J1 * 1e6))) + "_" +
    // std::to_string(static_cast<int>(std::round(J2 * 1e6))) + "_" +
    // std::to_string(static_cast<int>(std::round(hz * 1e6))) + ".h5";
    fmt(J1) + "_" +
    fmt(J2) + "_" +
    fmt(hz) + ".h5";
    bool loaded_from_file = false;

    //basis_print();
    //Hspin();

    //cout << "===============================\n";
    cout << "===============================\n";
    cout << "hz: " << hz << endl;
    cout << "Jx: " << J1 << endl;
    cout << "Jz: " << J2 << endl;
    cout << "Output file: " << outfile << endl;
    cout << "===============================\n";
    // cout << H << endl;
    //cout << "===============================\n";
    //==========================================================
    // Try loading existing eigenspectrum
    //==========================================================
    if (std::filesystem::exists(outfile))
    {
        try
        {
            cout << "===============================\n";
            cout << "Existing HDF5 file found.\n";
            cout << "Loading eigenspectrum...\n";
            cout << "===============================\n";

            HighFive::File file(outfile, HighFive::File::ReadOnly);

            file.getDataSet("eigenvalues").read(es);
            file.getDataSet("eigenvectors").read(evs);

            // Free Hamiltonian memory
            MatrixXcd().swap(H);

            loaded_from_file = true;

            cout << "HDF5 read successful\n";
        }
        catch (const std::exception& e)
        {
            cerr << "HDF5 Read Error: " << e.what() << endl;
            cerr << "Recomputing eigenspectrum...\n";
        }
    }

    //==========================================================
    // Compute if file absent or read failed
    //==========================================================
    if (!loaded_from_file)
    {
        Hspin();

        cout << "===============================\n";
        cout << "H generated\n";
        cout << "H(0,0): " << H(0,0) << endl;
        cout << "===============================\n";

        // basis_print();
        // cout << H << endl;

        cout << "===============================\n";
        cout << "Diagonalizing Hamiltonian...\n";
        cout << "===============================\n";

        // pair<MatrixXcd, VectorXd> e = Eigenspectrum(H);
        //
        // es = e.second;
        // evs = e.first;

        gpu_diagonalise(H, es, evs);

        if (!clean)
        {try
        {
            HighFive::File file(outfile, HighFive::File::Overwrite);

            file.createDataSet("eigenvalues", es);
            file.createDataSet("eigenvectors", evs);

            cout << "HDF5 write successful\n";
        }
        catch (const std::exception& e)
        {
            cerr << "HDF5 Write Error: " << e.what() << endl;
        }}

        // Free Hamiltonian after diagonalization
        MatrixXcd().swap(H);
    }

    //==========================================================
    // Observables
    //==========================================================
    // Sx();
    // Sz();
    // QFI(M_PI);


}
//=======================================================//
#endif
