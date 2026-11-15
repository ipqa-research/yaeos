module yaeos__math_linalg
   !! # Linear algebra helpers
   !! Dense linear solves, eigenproblems, and cubic-polynomial root utilities.
   !!
   !! # Description
   !! The module wraps LAPACK for linear systems and symmetric eigenproblems,
   !! and provides cubic-root solvers used by EOS implementations. `solve_system`
   !! solves \(Ax=b\); `eigen` computes eigenpairs of a real symmetric matrix.
   !!
   !! # Examples
   !! ```fortran
   !! x = solve_system(A, b)
   !! call eigen(A, values, vectors)
   !! ```
   use yaeos__constants, only: pr
   implicit none

contains
   function solve_system(a, b) result(x)
      !! Solve the dense linear system \(A x=b\) with LAPACK `dgesv`.
      !! The routine expects a square nonsingular matrix; LAPACK's status is
      !! currently not returned to the caller.
      use iso_fortran_env, only: error_unit
      real(pr), intent(in) :: b(:) !! Right-hand-side vector.
      real(pr), intent(in) :: a(size(b), size(b)) !! Square coefficient matrix.
      integer, parameter :: dp = selected_real_kind(15)

      real(pr) :: x(size(b))

      real(dp) :: a_lapack(size(b), size(b))
      real(dp) :: b_lapack(size(b))
      integer :: n
      integer :: nrhs
      integer :: lda
      integer :: ipiv(size(b))
      integer :: ldb
      integer :: info !! Solver status code.

      interface
         subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            !! # dgesv
            !!
            !! # Description
            !! LAPACK interface for dense linear solves using LU factorization.
            import dp
            integer :: n !! Order of the square coefficient matrix.
            integer :: nrhs !! Number of right-hand sides.
            real(dp) :: a(n, n) !! Coefficient matrix; overwritten by LU factors.
            integer :: lda !! Leading dimension of `a`.
            integer :: ipiv(n) !! Pivot indices from the factorization.
            real(dp) :: b(n) !! Right-hand side; overwritten by the solution.
            integer :: ldb !! Leading dimension of `b`.
            integer :: info !! LAPACK status code.
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
      !! # cubic_roots
      !!
      !! # Description
      !! Solve a real cubic polynomial with coefficients \((a,b,c,d)\).
      !!
      !! Uses the discriminant to distinguish three real roots, repeated real
      !! roots, and one real root with a complex-conjugate pair. Unused slots in
      !! the output arrays are set to NaN.
      use yaeos__auxiliar, only: sort
      real(pr), parameter :: pi=atan(1.0_pr) * 4.0_pr
      real(pr), intent(in) :: parameters(4) !! Polynomial coefficients \((a,b,c,d)\) in \(ax^3+bx^2+cx+d=0\).
      real(pr), intent(out) :: real_roots(3) !! Real roots; unused entries are NaN.
      complex(pr), intent(out) :: complex_roots(3) !! Reserved complex-root output; this solver currently leaves NaN sentinels.
      integer, intent(out) :: flag !! Root case: `0` repeated real roots, `1` one real plus conjugate pair, `-1` three distinct real roots.

      real(pr) :: p
      real(pr) :: q
      real(pr) :: u
      real(pr) :: v
      real(pr) :: nan
      real(pr) :: disc
      real(pr) :: theta

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
      !! # cubic_roots_rosendo
      !!
      !! # Description
      !! Solve a cubic polynomial using the trigonometric/Cardano formulation.
      !! This alternate implementation evaluates intermediate values in
      !! extended precision and returns real roots plus a root-case flag.
      use yaeos__auxiliar, only: sort
      real(pr), parameter :: pi=atan(1.0_pr) * 4.0_pr
      real(pr), intent(in) :: parameters(4) !! Coefficients \((a,b,c,d)\) of \(ax^3+bx^2+cx+d=0\).
      real(pr), intent(out) :: real_roots(3) !! Real roots; unused entries are not populated by this routine.
      complex(pr), intent(out) :: complex_roots(3) !! Present for interface compatibility; not assigned here.
      integer, intent(out) :: flag !! `-1` for three real roots, `1` for one real root.

      real(16) :: d1
      real(16) :: d2
      real(16) :: d3
      real(16) :: Q
      real(16) :: R !! Cardano intermediate \(R=(2d_1^3-9d_1d_2+27d_3)/54\) (not the gas constant).
      real(16) :: A
      real(16) :: B
      real(16) :: theta
      real(16) :: alp
      real(16) :: bet
      real(16) :: gam
      integer :: i !! Loop/index variable.

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
      !!
      !! # Examples
      !! ```fortran
      !! call eigen(A, lambda, vectors)
      !! ```
      !!
      !! # Notes
      !! The current implementation requests eigenvectors from LAPACK and
      !! assigns them to `eigenvectors`; provide this optional argument.

      ! use stdlib_linalg, only: eigh, linalg_state_type
      ! type(linalg_state_type) :: stat
      real(pr), intent(in out) :: A(:, :) !! Symmetric matrix; overwritten by LAPACK work data.
      real(pr), intent(out) :: eigenvalues(:) !! Eigenvalues in ascending order.
      real(pr), optional, intent(out) :: eigenvectors(:, :) !! Eigenvectors stored as columns.
      
      interface
         subroutine dsyev(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)
            !! LAPACK interface for a real symmetric eigenproblem.
            character(1) :: JOBZ !! Whether eigenvectors are computed.
            character(1) :: UPLO !! Triangle of the input matrix that is referenced.
            integer :: N !! Matrix order.
            double precision :: A(lda, *) !! Matrix/work array.
            integer :: lda !! Leading dimension of `A`.
            double precision :: W(*) !! Eigenvalues.
            double precision :: WORK(*) !! Workspace array.
            integer :: LWORK !! Workspace length.
            integer :: INFO !! LAPACK status code.
         end subroutine
      end interface
      integer :: istat
      real(pr) :: work(size(A, 1)*10) !! Numerical workspace.
      real(pr) :: Ain(size(A, 1), size(A, 1))
      integer :: n

      n = size(A, 1)
      Ain = A

      call dsyev(jobz='V', uplo='U', n=n, A=Ain, lda=n, w=eigenvalues, &
         work=work, lwork=size(work, 1), info=istat &
      )
      if (present(eigenvectors)) eigenvectors = Ain

      ! call eigh(&
      !    A=A, lambda=eigenvalues, vectors=eigenvectors, err=stat  &
      ! )
      
      ! if (.not. stat%ok()) then
      !    write(*, *) stat%print_msg()
      ! end if
   end subroutine

end module yaeos__math_linalg
