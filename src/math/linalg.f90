module yaeos__math_linalg
   !! Wrapper module around LAPACK's `dgesv`
   use yaeos_constants, only: pr
   implicit none

contains
   function solve_system(a, b) result(x)
      !! Solve a linear sytem AX = b
      real(pr), intent(in) :: b(:)
      real(pr), intent(in) :: a(size(b), size(b))

      real(pr) :: x(size(b))

      real(8) :: a_lapack(size(b), size(b)), b_lapack(size(b))
      integer :: n, nrhs, lda, ipiv(size(b)), ldb, info

      interface
         subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            integer :: n
            integer :: nrhs
            real(8) :: a(n, n)
            integer :: lda
            integer :: ipiv(n)
            real(8) :: b(n)
            integer :: ldb
            integer :: info
         end subroutine
      end interface

      n = size(a, dim=1)
      nrhs = 1
      lda = n
      ldb = n

      a_lapack = a
      b_lapack = b
      call dgesv(n, nrhs, a_lapack, lda, ipiv, b_lapack, ldb, info)

      if (info > 0) error stop 1

      x = b_lapack
   end function solve_system
end module