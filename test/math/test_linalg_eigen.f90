program test_linalg_eigen
   !! `eigen` must work with and without its optional eigenvectors argument,
   !! and return the known spectrum of a symmetric matrix.
   use yaeos__constants, only: pr
   use yaeos__math_linalg, only: eigen
   use testing_aux, only: assert, test_title
   implicit none
   real(pr) :: A(2, 2)
   real(pr) :: values(2)
   real(pr) :: vectors(2, 2)
   real(pr) :: residual(2)

   write(*, *) test_title("Symmetric eigenproblem")

   ! [[2, 1], [1, 2]] has eigenvalues 1 and 3.
   A = reshape([2.0_pr, 1.0_pr, 1.0_pr, 2.0_pr], [2, 2])
   call eigen(A, values, vectors)
   call assert(maxval(abs(values - [1.0_pr, 3.0_pr])) < 1e-12_pr, "eigenvalues in ascending order")
   residual = matmul(reshape([2.0_pr, 1.0_pr, 1.0_pr, 2.0_pr], [2, 2]), vectors(:, 2)) - 3.0_pr * vectors(:, 2)
   call assert(maxval(abs(residual)) < 1e-12_pr, "A v = lambda v for the largest eigenvalue")

   A = reshape([2.0_pr, 1.0_pr, 1.0_pr, 2.0_pr], [2, 2])
   call eigen(A, values)
   call assert(maxval(abs(values - [1.0_pr, 3.0_pr])) < 1e-12_pr, "eigenvalues without the optional eigenvectors")
end program test_linalg_eigen
