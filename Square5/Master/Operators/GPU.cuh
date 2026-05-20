#pragma once

// Plain C++ declaration — safe to include from .cpp/.hpp files
void launch_qfi_tloop(
    const double* dM2,
    const double* des,
    double* d_Z,
    double* d_QFI,
    int le, int N, int nT,
    double Tmin, double Tmax);
