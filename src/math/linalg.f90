module yaeos__math_linalg
   !! Wrapper module around LAPACK's `dgesv`
   use yaeos__constants, only: pr
   implicit none

contains
   function solve_system(a, b) result(x)
      !! Solve a linear sytem AX = b
      use iso_fortran_env, only: error_unit
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

      ! if (info > 0) write(error_unit, *) "WARN: dgesv error"

      x = b_lapack
   end function solve_system

   subroutine cubic_roots(parameters, real_roots, complex_roots, flag)
      use yaeos__auxiliar, only: sort
      real(pr), parameter :: pi=atan(1.0_pr) * 4.0_pr
      real(pr), intent(in) :: parameters(4)
      real(pr), intent(out) :: real_roots(3)
      complex(pr), intent(out) :: complex_roots(3)
      integer, intent(out) :: flag
      !! flag that identifies which case the solution is
      !! - `0`: 3 real rotos, one of them repeated (use real_roots(1) and real_roots(2))
      !! - `1`: 1 real root, 2 complex roots.
      !!   Use real_roots(1) and complex_roots(1) and complex_roots(2)
      !! - `-1`: 3 real roots, all different

      real(pr) :: p, q, u, v, nan
      real(pr) :: disc, theta

      nan = 0
      nan = nan/nan

      associate(&
         a => parameters(1), b => parameters(2), &
         c => parameters(3), d => parameters(4)&
         )

         p = c/a - b**2/(3*a**2)
         q = d/a - b*c/(3*a**2) + 2*b**3/(27*a**3)

         disc = q**2 + 4 * p**3 / 27
         real_roots = nan
         complex_roots = nan

         if (abs(disc) < 1e-15) then
            flag = 0
            real_roots(1) = 3*q/p
            real_roots(2) = -3*q/(2*p)
            real_roots(3) = real_roots(2)
         elseif (disc < 0) then
            flag = -1
            theta = acos(0.5_pr * 3 * q / p * sqrt(-3/p))
            real_roots(1) = 2 * sqrt(-p/3) * cos(theta/3)
            real_roots(2) = 2 * sqrt(-p/3) * cos((theta + 2*pi)/3)
            real_roots(3) = 2 * sqrt(-p/3) * cos((theta + 4*pi)/3)
            call sort(real_roots)
         elseif (disc > 0) then
            flag = 1
            u = ((-q + sqrt(disc))/2)
            v = ((-q - sqrt(disc))/2)

            u = sign(abs(u)**(1.0_pr/3.0_pr), u)
            v = sign(abs(v)**(1.0_pr/3.0_pr), v)
            real_roots(1) = u + v
         endif
         real_roots = real_roots - b/(3*a)
      end associate
   end subroutine cubic_roots

   subroutine cubic_roots_rosendo(parameters, real_roots, complex_roots, flag)
      use yaeos__auxiliar, only: sort
      real(pr), parameter :: pi=atan(1.0_pr) * 4.0_pr
      real(pr), intent(in) :: parameters(4)
      real(pr), intent(out) :: real_roots(3)
      complex(pr), intent(out) :: complex_roots(3)
      integer, intent(out) :: flag

      real(16) :: d1, d2, d3, Q, R, A, B, theta, alp, bet, gam
      integer :: i

      d1 = parameters(2) / parameters(1)
      d2 = parameters(3) / parameters(1)
      d3 = parameters(4) / parameters(1)

      Q = (d1**2 - 3*d2) / 9.0_16
      R = (2*d1**3 - 9*d1*d2 + 27*d3) / 54.0_16

      if (R**2 <= Q**3) then
         theta = acos(R / sqrt(Q**3))

         real_roots(1) = -2 * sqrt(Q) * cos(theta / 3.0_16) - d1 / 3.0_16
         real_roots(2) = -2 * sqrt(Q) * cos((theta + 2 * pi) / 3.0_16) - d1 / 3.0_16
         real_roots(3) = -2 * sqrt(Q) * cos((theta - 2 * pi) / 3.0_16) - d1 / 3.0_16

         ! Correction??
         ! do i=1,100
         !    real_roots(1) = -d1 - (real_roots(2) + real_roots(3))
         !    real_roots(2) = (d2 - real_roots(1) * real_roots(3)) / (real_roots(1) + real_roots(3))
         !    real_roots(3) = -d3 / (real_roots(1) * real_roots(2))
         ! end do

         call sort(real_roots)
         flag = -1
      else
         A = - sign((abs(R) + sqrt(R**2 - Q**3))**(1.0_16/3.0_16), R)

         if (abs(A) < 1e-6) then
            A = 0.0_16
            B = 0.0_16
         else
            B = Q / A
         end if

         real_roots = (A + B) - d1 / 3.0_16
         flag = 1
      end if
   end subroutine

   subroutine eigen(A, eigenvalues, eigenvectors)
      !! # eigen
      !!
      !! # Description
      !! Calculate the eigenvalues and eigenvectors of a real symmetric matrix
      !! `A` using LAPACK's `dsyev`. The eigenvectors are stored in the columns
      !! of `eigenvectors`. The eigenvalues are stored in `eigenvalues`.

      ! use stdlib_linalg, only: eigh, linalg_state_type
      ! type(linalg_state_type) :: stat
      real(pr), intent(in out) :: A(:, :)
      real(pr), intent(out) :: eigenvalues(:)
      real(pr), optional, intent(out) :: eigenvectors(:, :)
      
      interface
         subroutine dsyev(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO) 	
            character(1) :: JOBZ
            character(1) :: UPLO
            integer :: N
            double precision :: A(lda, *)
            integer :: lda
            double precision :: W(*)
            double precision :: WORK(*)
            integer :: LWORK
            integer :: INFO
         end subroutine
      end interface
      integer :: istat
      real(pr) :: work(size(A, 1)*10)
      real(pr) :: Ain(size(A, 1), size(A, 1))
      integer :: n

      n = size(A, 1)
      Ain = A

      call dsyev(jobz='V', uplo='U', n=n, A=Ain, lda=n, w=eigenvalues, &
         work=work, lwork=size(work, 1), info=istat &
      )
      eigenvectors = Ain

      ! call eigh(&
      !    A=A, lambda=eigenvalues, vectors=eigenvectors, err=stat  &
      ! )
      
      ! if (.not. stat%ok()) then
      !    write(*, *) stat%print_msg()
      ! end if
   end subroutine

end module yaeos__math_linalg
