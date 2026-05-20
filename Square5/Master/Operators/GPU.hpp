#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>

static void gpu_zgemm_adjoint(
    int n,
    const MatrixXcd& A,   // evs  (n×n)
    const MatrixXcd& B,   // P    (n×n)
    MatrixXcd&       C)   // M    (n×n), output
{
    cublasHandle_t handle;
    cublasCreate(&handle);

    // Allocate device memory
    cuDoubleComplex *dA, *dB, *dC;
    size_t bytes = n * n * sizeof(cuDoubleComplex);
    cudaMalloc(&dA, bytes);
    cudaMalloc(&dB, bytes);
    cudaMalloc(&dC, bytes);

    // Upload (Eigen is column-major, cuBLAS expects column-major — matches)
    cublasSetMatrix(n, n, sizeof(cuDoubleComplex),
                    A.data(), n, dA, n);
    cublasSetMatrix(n, n, sizeof(cuDoubleComplex),
                    B.data(), n, dB, n);

    cuDoubleComplex alpha = {1.0, 0.0};
    cuDoubleComplex beta  = {0.0, 0.0};

    // M = A^H * B  →  CUBLAS_OP_C = conjugate transpose
    cublasZgemm(handle,
                CUBLAS_OP_C, CUBLAS_OP_N,
                n, n, n,
                &alpha,
                dA, n,   // A^H
                dB, n,   // B
                &beta,
                dC, n);  // C = output

    // Download result
    C.resize(n, n);
    cublasGetMatrix(n, n, sizeof(cuDoubleComplex),
                    dC, n, C.data(), n);

    cudaFree(dA); cudaFree(dB); cudaFree(dC);
    cublasDestroy(handle);
}

void gpu_diagonalise(
    MatrixXcd&  H,    // input Hamiltonian (Hermitian), overwritten
    VectorXd&   evals,
    MatrixXcd&  evecs)
{
    cusolverDnHandle_t handle;
    cusolverDnCreate(&handle);

    int n = H.rows();
    size_t bytes = n * n * sizeof(cuDoubleComplex);

    cuDoubleComplex* dA;
    double*          dW;   // eigenvalues
    cudaMalloc(&dA, bytes);
    cudaMalloc(&dW, n * sizeof(double));

    cublasSetMatrix(n, n, sizeof(cuDoubleComplex),
                    H.data(), n, dA, n);

    // Query workspace size
    int lwork = 0;
    cusolverDnZheevd_bufferSize(
        handle,
        CUSOLVER_EIG_MODE_VECTOR,  // compute eigenvectors
        CUBLAS_FILL_MODE_LOWER,    // use lower triangle
        n, dA, n, dW, &lwork);

    cuDoubleComplex* dWork;
    int*             dInfo;
    cudaMalloc(&dWork, lwork * sizeof(cuDoubleComplex));
    cudaMalloc(&dInfo, sizeof(int));

    // Solve
    cusolverDnZheevd(
        handle,
        CUSOLVER_EIG_MODE_VECTOR,
        CUBLAS_FILL_MODE_LOWER,
        n, dA, n, dW,
        dWork, lwork, dInfo);

    // Check convergence
    int info = 0;
    cudaMemcpy(&info, dInfo, sizeof(int), cudaMemcpyDeviceToHost);
    if (info != 0)
        throw std::runtime_error("cusolverDnZheevd failed, info=" + std::to_string(info));

    // Download
    evals.resize(n);
    evecs.resize(n, n);
    cudaMemcpy(evals.data(), dW, n * sizeof(double), cudaMemcpyDeviceToHost);
    cublasGetMatrix(n, n, sizeof(cuDoubleComplex),
                    dA, n, evecs.data(), n);

    cudaFree(dA); cudaFree(dW); cudaFree(dWork); cudaFree(dInfo);
    cusolverDnDestroy(handle);
}
