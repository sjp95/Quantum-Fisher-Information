#include "GPU.cuh"

__global__ void qfi_tloop_kernel(
    const double* __restrict__ dM2,   // le×le, row-major
    const double* __restrict__ des,   // eigenvalues, length le
    double* __restrict__ out_Z,       // one value per T
    double* __restrict__ out_QFI,     // one value per T
    int    le,
    int    N,
    int    nT,
    double Tmin,
    double Tmax)
{
    // Each block handles one T value
    int ti = blockIdx.x;
    if (ti >= nT) return;

    double t    = Tmin * pow(Tmax / Tmin, (double)ti / (double)(nT - 1));
    double beta = 1.0 / t;
    double es0  = des[0];

    // ---- Partition function: reduce over states ----
    __shared__ double sZ[256];
    double local_Z = 0.0;
    for (int s = threadIdx.x; s < le; s += blockDim.x)
        local_Z += exp(-(des[s] - es0) * beta);

    sZ[threadIdx.x] = local_Z;
    __syncthreads();
    for (int s = blockDim.x/2; s > 0; s >>= 1) {
        if (threadIdx.x < s) sZ[threadIdx.x] += sZ[threadIdx.x + s];
        __syncthreads();
    }
    double Z = sZ[0];

    // ---- QFI triangular sum ----
    __shared__ double sQ[256];
    double local_Q = 0.0;

    long long npairs = (long long)le * (le - 1) / 2;
    for (long long k = threadIdx.x; k < npairs; k += blockDim.x)
    {
        // Recover l1, l2 from flat index
        long long l1 = (long long)((2.0*le - 1.0 -
                        sqrt((2.0*le - 1.0)*(2.0*le - 1.0) - 8.0*k)) / 2.0);
        long long l2 = k - l1*(le - 1) + l1*(l1 - 1)/2 + l1 + 1;

        double p1    = exp(-(des[l1] - es0) * beta);
        double p2    = exp(-(des[l2] - es0) * beta);
        double denom = p1 + p2;
        if (denom < 1e-14) continue;
        double contrib = (p1 - p2) * (p1 - p2) / denom;
        if (contrib < 1e-12) continue;
        // local_Q += 4.0 * dM2[l1 * le + l2] * contrib;
        local_Q += 4.0 * dM2[l2 * le + l1] * contrib;
    }

    sQ[threadIdx.x] = local_Q;
    __syncthreads();
    for (int s = blockDim.x/2; s > 0; s >>= 1) {
        if (threadIdx.x < s) sQ[threadIdx.x] += sQ[threadIdx.x + s];
        __syncthreads();
    }

    if (threadIdx.x == 0) {
        out_Z[ti]   = Z;
        out_QFI[ti] = sQ[0] / ((double)N * Z);
    }
}


void launch_qfi_tloop(
    const double* dM2,
    const double* des,
    double* d_Z,
    double* d_QFI,
    int le, int N, int nT,
    double Tmin, double Tmax)
{
    qfi_tloop_kernel<<<nT, 256>>>(dM2, des, d_Z, d_QFI, le, N, nT, Tmin, Tmax);
    cudaDeviceSynchronize();
}
