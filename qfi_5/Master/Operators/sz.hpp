#ifndef Sz_HPP_INCLUDED
#define Sz_HPP_INCLUDED

#include "../input.hpp"
#include <Eigen/Dense>
// #include <highfive/H5File.hpp>
// #include <highfive/eigen.hpp>
// #include <vector>
// #include <cmath>
// #include <complex>
#include <omp.h>
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <filesystem>
// #include "../creat_directory.hpp"

// #include <cublas_v2.h>
#include <cuda_runtime.h>
// #include <cusolverDn.h>

#include "GPU.hpp"
#include "GPU.cuh"

// #ifndef Sz_HPP_INCLUDED
// #define Sz_HPP_INCLUDED
//
// #include "../input.hpp"
// #include <Eigen/Dense>
// #include <omp.h>
// #include "GPU.cuh"

using namespace std;
using namespace Eigen;

void input::gSz_GPU()
{
    createDirectory("../Data/Magnetization/Sz");

    std::string outfile =
    "../Data/Magnetization/Sz/Sz_" +
    std::to_string(N) + "_" +
    fmt(J1) + "_" +
    fmt(J2) + "_" +
    fmt(hz) + ".dat";

    std::ofstream file(outfile);
    file << std::scientific << std::setprecision(17);

    // Precompute squared magnitudes of ground state eigenvector on CPU
    VectorXd weights = evs.col(0).cwiseAbs2();

    // GPU acceleration
    double *d_weights, *d_Sz;

    cudaMalloc(&d_weights, le * sizeof(double));
    cudaMalloc(&d_Sz, sizeof(double));

    // Upload precomputed weights
    cudaMemcpy(d_weights, weights.data(), le * sizeof(double), cudaMemcpyHostToDevice);

    // Launch GPU kernel
    launch_sz_gpu(d_weights, d_Sz, le, N, Ls);

    // Download result
    double h_Sz = 0.0;
    cudaMemcpy(&h_Sz, d_Sz, sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_weights);
    cudaFree(d_Sz);

    file << hz << " "
         << J1 << " "
         << J2 << " "
         << h_Sz << endl;
}

// void input::gSz_CPU()
// {
//     createDirectory("../Data/Magnetization/Sz");
//
//     std::string outfile =
//     "../Data/Magnetization/Sz/Sz_" +
//     std::to_string(N) + "_" +
//     fmt(J1) + "_" +
//     fmt(J2) + "_" +
//     fmt(hz) + ".dat";
//
//     std::ofstream file(outfile);
//     file << std::scientific << std::setprecision(17);
//
//     double total = 0.0;
//
//     #pragma omp parallel for reduction(+:total) schedule(static)
//     for (int n = ls; n < le; n++)
//     {
//         int state = n;
//         double weight = std::norm(evs(n,0));
//         double local = 0.0;
//
//         for (int i = 0; i < N; i++)
//         {
//             int reminder = state % Ls;
//             state /= Ls;
//             local += (reminder - 0.5);
//         }
//
//         total += local * weight;
//     }
//
//     file << hz << " "
//          << J1 << " "
//          << J2 << " "
//          << total/double(N) << endl;
// }

void input::gSz_CPU()
{
    createDirectory("../Data/Magnetization/Sz");

    std::string outfile =
    "../Data/Magnetization/Sz/Sz_" +
    std::to_string(N) + "_" +
    fmt(J1) + "_" +
    fmt(J2) + "_" +
    fmt(hz) + ".dat";

    std::ofstream file(outfile);
    file << std::scientific << std::setprecision(17);

    double total = 0.0;

    #pragma omp parallel for reduction(+:total) schedule(static)
    for (int n = ls; n < le; n++)
    {
        int state = n;
        double weight = std::norm(gs_vector(n));
        double local = 0.0;

        for (int i = 0; i < N; i++)
        {
            int reminder = state % Ls;
            state /= Ls;
            local += (reminder - 0.5);
        }

        total += local * weight;
    }

    file << hz << " "
         << J1 << " "
         << J2 << " "
         << total/double(N) << endl;
}

#endif
