module yaeos__math_linalg
   !! Wrapper module around LAPACK's `dgesv`
   use yaeos__constants, only: pr
   implicit none

contains
   function solve_system(a, b) result(x)
      !! Solve a linear sytem AX = b
      real(pr), intent(in) :: b(:)
      real(pr), intent(in) :: a(size(b), size(b))
      integer, parameter :: dp = selected_real_kind(15)

      real(pr) :: x(size(b))

      real(dp) :: a_lapack(size(b), size(b)), b_lapack(size(b))
      integer :: n, nrhs, lda, ipiv(size(b)), ldb, info

      interface
         subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import dp
            integer :: n
            integer :: nrhs
            real(dp) :: a(n, n)
            integer :: lda
            integer :: ipiv(n)
            real(dp) :: b(n)
            integer :: ldb
            integer :: info
         end subroutine dgesv
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

   subroutine cubic_roots(a, b, c, d, r1, r2, r3, flag)
      real(pr), parameter :: pi=atan(1.0_pr)*4
      real(pr), intent(in) :: a, b, c, d
      real(pr), intent(out) :: r1, r2, r3
      integer, intent(out) :: flag

      real(pr) :: p, q, u, v
      real(pr) :: disc

      real(pr) :: z1, z2, z3, theta

      p = c/a - b**2/(3*a**2)
      q = d/a - b*c/(3*a**2) + 2*b**3/(27*a**3)

      disc = q**2 + 4*p**3 / 27

      z1 = 0
      z2 = 0
      z3 = 0
      r1 = 0
      r2 = 0
      r3 = 0

      if (abs(disc) < 1e-15) then
         flag = 0
         z1 = 3*q/p
         z2 = -3*q/(2*p)
         z3 = z2
      elseif (disc < 0) then
         flag = -1
         theta = acos(0.5_pr * 3 * q / p * sqrt(-3/p))
         z1 = 2 * sqrt(-p/3) * cos(theta/3)
         z2 = 2 * sqrt(-p/3) * cos((theta + 2*pi)/3)
         z3 = 2 * sqrt(-p/3) * cos((theta + 4*pi)/3)
      elseif (disc > 0) then
         flag = 1
         u = ((-q + sqrt(disc))/2)
         v = ((-q - sqrt(disc))/2)

         u = sign(abs(u)**(1.0_pr/3.0_pr), u)
         v = sign(abs(v)**(1.0_pr/3.0_pr), v)
         z1 = u + v
      endif

      r1 = z1 - b/(3*a)
      r2 = z2 - b/(3*a)
      r3 = z3 - b/(3*a)
   end subroutine

end module yaeos__math_linalg
