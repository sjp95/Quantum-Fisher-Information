#ifndef QFI_HPP_INCLUDED
#define QFI_HPP_INCLUDED

#include "../input.hpp"
#include <Eigen/Dense>
#include <highfive/H5File.hpp>
#include <highfive/eigen.hpp>
#include <vector>
#include <cmath>
#include <complex>
#include <omp.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include "../creat_directory.hpp"

#include "GPU.hpp"
#include "GPU.cuh"

using namespace std;
using namespace Eigen;

//================================================//
// Quantum Fisher Information for spin chain
//
// Generator:
// O(qx) = sum_x Sx(x) exp(-i qx x)
//
// Site mapping:
// x = site  (linear chain)
//
// Matrix M(l1,l2) = <l1|O(qx)|l2> is temperature
// independent and cached to HDF5. QFI sum is then
// cheap to re-run for any T.
//================================================//

void input::QFI(double qx)
{
    createDirectory("../Data/QFI");

    //------------------------------------------------------------
    // Matrix file  — T-independent, no qx in name (almost always pi)
    //------------------------------------------------------------
    std::string matfile =
        "../Data/QFI/QFI_Matrix_" +
        std::to_string(N) + "_" +
        fmt(J1) + "_" +
        fmt(J2) + "_" +
        fmt(hz) + ".h5";

    //------------------------------------------------------------
    // Output file — T-independent name, one line per T appended
    //------------------------------------------------------------
    std::string outfile =
        "../Data/QFI/QFI_" +
        std::to_string(N) + "_" +
        fmt(J1) + "_" +
        fmt(J2) + "_" +
        fmt(hz) + ".dat";

    std::ofstream file(outfile, std::ios::app);
    file << std::scientific << std::setprecision(17);

    const complex<double> II(0.0, 1.0);
    const double beta = 1.0 / T;

    //============================================================
    // PHASE 1 — Build or load M(l1,l2) = <l1|O(qx)|l2>
    //============================================================

    MatrixXcd M(le, le);
    bool matrix_loaded = false;

    //------------------------------------------------------------
    // Try loading from HDF5
    //------------------------------------------------------------
    if (std::filesystem::exists(matfile))
    {
        try
        {
            cout << "===============================\n";
            cout << "QFI matrix file found.\n";
            cout << "Loading M from HDF5...\n";
            cout << "===============================\n";

            HighFive::File h5(matfile, HighFive::File::ReadOnly);
            h5.getDataSet("M").read(M);

            matrix_loaded = true;
            cout << "QFI matrix loaded successfully.\n";
        }
        catch (const std::exception& e)
        {
            cerr << "HDF5 Read Error: " << e.what() << endl;
            cerr << "Rebuilding QFI matrix...\n";
        }
    }

    //------------------------------------------------------------
    // Build M if not loaded
    //------------------------------------------------------------
    if (!matrix_loaded)
    {
        cout << "===============================\n";
        cout << "Building QFI operator matrix...\n";
        cout << "===============================\n";

        //--------------------------------------------------------
        // Step 1 — Build P(q, l2) in computational basis:
        //
        // Apply O(qx) to each eigenvector column:
        // P(q, l2) = sum_{m, x: flip_x(m)=q} 0.5 * exp(-i qx x) * evs(m, l2)
        //
        // For each basis state m and each site x:
        //   q = flip_x(m)  using Lspow[x]
        //   P(q, l2) += 0.5 * phase(x) * evs(m, l2)
        //
        // Parallelise over l2 — each column is independent.
        //--------------------------------------------------------

        MatrixXcd P = MatrixXcd::Zero(le, le);

        std::vector<std::vector<int>> spins(le, std::vector<int>(N));
        for (int m = 0; m < le; m++) {
            int temp = m;
            for (int j = 0; j < N; j++) {
                spins[m][j] = temp % Ls;
                temp       /= Ls;
            }
        }

        // Precompute Fourier phases per site
        std::vector<complex<double>> phase(N);
        for (int x = 0; x < N; x++)
            phase[x] = 0.5 * exp(-II * qx * double(x));

        #pragma omp parallel for schedule(dynamic)
        for (int l2 = 0; l2 < le; l2++) {
            for (int m = 0; m < le; m++) {
                for (int x = 0; x < N; x++) {
                    int q = (spins[m][x] == 0) ? m + Lspow[x] : m - Lspow[x];
                    P(q, l2) += phase[x] * evs(m, l2);
                }
            }
        }

        //--------------------------------------------------------
        // Step 2 — Project into eigenbasis
        // M(l1,l2) = <l1|O|l2> = evs† · P
        //--------------------------------------------------------
        cout << "===============================\n";
        cout << "Projecting into eigenbasis using gpu...\n";
        cout << "===============================\n";

        gpu_zgemm_adjoint(le, evs, P, M);

        //--------------------------------------------------------
        // Step 3 — Write M to HDF5
        //--------------------------------------------------------
        try
        {
            HighFive::File h5(matfile, HighFive::File::Overwrite);
            h5.createDataSet("M", M);
            cout << "QFI matrix written to HDF5 successfully.\n";
        }
        catch (const std::exception& e)
        {
            cerr << "HDF5 Write Error: " << e.what() << endl;
        }
    }

    //============================================================
    // PHASE 2 — Thermal QFI sum using M
    //
    // F_Q = (4 / N*Z) * sum_{l1 < l2}
    //           (p1-p2)^2 / (p1+p2) * |M(l1,l2)|^2
    //============================================================

    //------------------------------------------------------------
    // Partition function
    //------------------------------------------------------------
    double Z = 0.0;

    #pragma omp parallel for reduction(+:Z)
    for (int state = 0; state < le; state++)
        Z += exp(-(es(state) - es(0)) * beta);

    //------------------------------------------------------------
    // Double sum over eigenpairs
    //------------------------------------------------------------
    double QFI_sum = 0.0;

    #pragma omp parallel for reduction(+:QFI_sum) schedule(dynamic)
    for (int l1 = 0; l1 < le; l1++)
    {
        for (int l2 = l1 + 1; l2 < le; l2++)
        {
            double p1 = exp(-(es(l1) - es(0)) * beta);
            double p2 = exp(-(es(l2) - es(0)) * beta);

            double numerator   = (p1 - p2) * (p1 - p2);
            double denominator = p1 + p2;

            if (denominator < 1e-14) continue;
            if (numerator / denominator < 1e-12) continue;

            QFI_sum +=
                4.0 *
                norm(M(l1, l2)) *
                numerator / denominator;
        }
    }

    double qfi_density = QFI_sum / (double(N) * Z);

    cout << "===============================\n";
    cout << "QFI = " << qfi_density << endl;
    cout << "===============================\n";

    file << T   << " "
         << hz  << " "
         << J1  << " "
         << J2  << " "
         << qfi_density << endl;

    file.close();
}

void input::tQFI(double qx, double Tmin, double Tmax, int nT, bool clean)
{
    createDirectory("../Data/QFI");

    //------------------------------------------------------------
    // Matrix file — T-independent
    //------------------------------------------------------------
    // std::string matfile =
    //     "../Data/QFI/QFI_Matrix_" +
    //     std::to_string(N) + "_" +
    //     std::to_string(static_cast<int>(std::round(J1 * 100.0))) + "_" +
    //     std::to_string(static_cast<int>(std::round(J2 * 100.0))) + "_" +
    //     std::to_string(static_cast<int>(std::round(hz * 100.0))) + ".h5";
    std::string matfile =
        "../Data/QFI/QFI_Matrix" +
        std::to_string(N) + "_" +
        fmt(J1) + "_" +
        fmt(J2) + "_" +
        fmt(hz) + ".h5";

    //------------------------------------------------------------
    // Output file — one line per T, truncated fresh each run
    //------------------------------------------------------------
    // std::string outfile =
    //     "../Data/QFI/QFI_" +
    //     std::to_string(N) + "_" +
    //     std::to_string(static_cast<int>(std::round(J1 * 100.0))) + "_" +
    //     std::to_string(static_cast<int>(std::round(J2 * 100.0))) + "_" +
    //     std::to_string(static_cast<int>(std::round(hz * 100.0))) + ".dat";
    std::string outfile =
        "../Data/QFI/QFI_" +
        std::to_string(N) + "_" +
        fmt(J1) + "_" +
        fmt(J2) + "_" +
        fmt(hz) + ".dat";

    // Truncate once at the start — clean slate for this run
    { std::ofstream trunc(outfile, std::ios::trunc); }

    const complex<double> II(0.0, 1.0);

    //============================================================
    // PHASE 1 — Build or load M(l1,l2) = <l1|O(qx)|l2>
    //           Happens ONCE regardless of how many T values
    //============================================================

    MatrixXcd M(le, le);
    bool matrix_loaded = false;

    //------------------------------------------------------------
    // Try loading from HDF5
    //------------------------------------------------------------
    if (std::filesystem::exists(matfile))
    {
        try
        {
            cout << "===============================\n";
            cout << "QFI matrix file found.\n";
            cout << "Loading M from HDF5...\n";
            cout << "===============================\n";

            HighFive::File h5(matfile, HighFive::File::ReadOnly);
            h5.getDataSet("M").read(M);

            matrix_loaded = true;
            cout << "QFI matrix loaded successfully.\n";
        }
        catch (const std::exception& e)
        {
            cerr << "HDF5 Read Error: " << e.what() << endl;
            cerr << "Rebuilding QFI matrix...\n";
        }
    }

    //------------------------------------------------------------
    // Build M if not loaded
    //------------------------------------------------------------
    if (!matrix_loaded)
    {
        cout << "===============================\n";
        cout << "Building QFI operator matrix in tQFI...\n";
        cout << "===============================\n";

        MatrixXcd P = MatrixXcd::Zero(le, le);

        std::vector<std::vector<int>> spins(le, std::vector<int>(N));
        for (int m = 0; m < le; m++) {
            int temp = m;
            for (int j = 0; j < N; j++) {
                spins[m][j] = temp % Ls;
                temp       /= Ls;
            }
        }

        // Precompute Fourier phases per site
        std::vector<complex<double>> phase(N);
        for (int x = 0; x < N; x++)
            phase[x] = 0.5 * exp(-II * qx * double(x));

        #pragma omp parallel for schedule(dynamic)
        for (int l2 = 0; l2 < le; l2++) {
            for (int m = 0; m < le; m++) {
                for (int x = 0; x < N; x++) {
                    int q = (spins[m][x] == 0) ? m + Lspow[x] : m - Lspow[x];
                    P(q, l2) += phase[x] * evs(m, l2);
                }
            }
        }

        cout << "===============================\n";
        cout << "Projecting into eigenbasis using gpu...\n";
        cout << "===============================\n";

        // omp_set_num_threads(1);
        // M.noalias() = evs.adjoint() * P;
        // omp_set_num_threads(omp_get_max_threads());

        gpu_zgemm_adjoint(le, evs, P, M);

        try
        {
            HighFive::File h5(matfile, HighFive::File::Overwrite);
            h5.createDataSet("M", M);
            cout << "QFI matrix written to HDF5 successfully.\n";
        }
        catch (const std::exception& e)
        {
            cerr << "HDF5 Write Error: " << e.what() << endl;
        }
    }

    //============================================================
    // PHASE 2 — Precompute |M(l1,l2)|^2 once, reuse for all T
    //           Avoids repeated norm() calls inside the T loop
    //============================================================

    MatrixXd M2 = M.cwiseAbs2();

    // M is no longer needed — free memory before T loop
    MatrixXcd().swap(M);

    //============================================================
    // PHASE 3 — T sweep entirely in memory, no further disk I/O
    //
    // F_Q = (4 / N*Z) * sum_{l1 < l2}
    //           (p1-p2)^2 / (p1+p2) * |M(l1,l2)|^2
    //============================================================

    std::ofstream file(outfile, std::ios::app);
    file << std::scientific << std::setprecision(17);

    // Upload once (dM2 may already be on GPU from multiply step)
    double *des, *dM2_gpu;
    cudaMalloc(&des,      le * sizeof(double));
    cudaMalloc(&dM2_gpu,  le * le * sizeof(double));
    cudaMemcpy(des,      es.data(),  le * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dM2_gpu,  M2.data(),  le * le * sizeof(double), cudaMemcpyHostToDevice);

    // Output buffers
    double *d_Z, *d_QFI;
    cudaMalloc(&d_Z,   nT * sizeof(double));
    cudaMalloc(&d_QFI, nT * sizeof(double));

    // Launch — one block per T, 256 threads per block
    // qfi_tloop_kernel<<<nT, 256>>>(dM2_gpu, des, d_Z, d_QFI, le, N, nT, Tmin, Tmax);
    launch_qfi_tloop(dM2_gpu, des, d_Z, d_QFI, le, N, nT, Tmin, Tmax);

    // Single download of all results
    std::vector<double> h_Z(nT), h_QFI(nT);
    cudaMemcpy(h_Z.data(),   d_Z,   nT * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_QFI.data(), d_QFI, nT * sizeof(double), cudaMemcpyDeviceToHost);

    // check starts --->
    int check_i = 0;
    double t     = Tmin * pow(Tmax / Tmin, double(check_i) / double(nT - 1));
    double beta  = 1.0 / t;

    double Z_cpu = 0.0;
    for (int s = 0; s < le; s++)
        Z_cpu += exp(-(es(s) - es(0)) * beta);

    double QFI_cpu = 0.0;
    for (int l1 = 0; l1 < le; l1++) {
        for (int l2 = l1+1; l2 < le; l2++) {
            double p1 = exp(-(es(l1) - es(0)) * beta);
            double p2 = exp(-(es(l2) - es(0)) * beta);
            double denom = p1 + p2;
            if (denom < 1e-14) continue;
            double contrib = (p1-p2)*(p1-p2)/denom;
            if (contrib < 1e-12) continue;
            QFI_cpu += 4.0 * M2(l1, l2) * contrib;
        }
    }
    QFI_cpu /= double(N) * Z_cpu;

    std::cout << "T=" << t << "\n"
            << "  GPU Z   = " << h_Z[check_i]   << "  CPU Z   = " << Z_cpu   << "\n"
            << "  GPU QFI = " << h_QFI[check_i] << "  CPU QFI = " << QFI_cpu << "\n";

    // Print first 10 pairs
    int count = 0;
    for (int l1 = 0; l1 < le && count < 10; l1++) {
        for (int l2 = l1+1; l2 < le && count < 10; l2++) {
            double p1 = exp(-(es(l1) - es(0)) * beta);
            double p2 = exp(-(es(l2) - es(0)) * beta);
            double denom = p1 + p2;
            if (denom < 1e-14) continue;
            double contrib = (p1-p2)*(p1-p2)/denom;
            printf("CPU pair (%d,%d): p1=%.6e p2=%.6e contrib=%.6e M2=%.6e term=%.6e\n",
                l1, l2, p1, p2, contrib, M2(l1,l2), 4.0*M2(l1,l2)*contrib);
            count++;
        }
    }
    // check ends <-----

    // Write output — now just file I/O
    for (int i = 0; i < nT; i++)
    {
        double t = Tmin * pow(Tmax / Tmin, double(i) / double(nT - 1));
        file << t         << " "
            << hz        << " "
            << J1        << " "
            << J2        << " "
            << h_QFI[i]  << "\n";
    }

    cudaFree(des); cudaFree(dM2_gpu);
    cudaFree(d_Z); cudaFree(d_QFI);

    file.close();

    if (clean)
    {try {
        if (std::filesystem::remove(matfile)) {
            std::cout << "File " << matfile << " deleted successfully\n";
        } else {
            std::cout << "File " << matfile << " does not exist\n";
        }
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error: " << e.what() << '\n';
    }}
}

#endif
