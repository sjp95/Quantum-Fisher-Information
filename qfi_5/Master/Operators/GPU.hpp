#ifndef GPU_HPP
#define GPU_HPP

void gpu_zgemm_adjoint(
    int n,
    const MatrixXcd& A,   // evs  (n×n)
    const MatrixXcd& B,   // P    (n×n)
    MatrixXcd&       C);   // M    (n×n), output

void gpu_diagonalise(
    MatrixXcd&  H,    // input Hamiltonian (Hermitian), overwritten
    VectorXd&   evals,
    MatrixXcd&  evecs);

#endif
