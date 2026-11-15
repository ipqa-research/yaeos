program test_math
   use yaeos__constants, only: pr
   use auxiliar_functions, only: allclose
   use testing_aux, only: assert, test_title
   implicit none

   write(*, *) test_title("Math")
   call test_dxk_dn()
   call test_dx_to_dn()
   call test_sq_error()
   call test_cardano_method()
   call test_rosendo_method()
   call test_newton_method()
   write(*, *) ""

contains

   subroutine test_dxk_dn()
      use yaeos__math, only: derivative_dxk_dni, derivative_d2xk_dnidnj

      integer, parameter :: nc = 3

      real(pr) :: n(nc), n_tot
      real(pr) :: dxk_dn(nc, nc), dxk_dni
      real(pr) :: d2xk_dn(nc, nc, nc), d2xk_dninj

      integer i, j, k

      n = [10.0_pr, 2.543_pr, 3.123_pr]
      n_tot = sum(n)

      dxk_dn = derivative_dxk_dni(n)
      d2xk_dn = derivative_d2xk_dnidnj(n)

      ! First derivative
      do k = 1, nc
         do i = 1, nc
            if (k == i) then
               dxk_dni = (n_tot - n(i)) / n_tot**2
            else
               dxk_dni = -n(k) / n_tot**2
            end if

            call assert(abs(dxk_dn(k, i) - dxk_dni) < 1e-10_pr, &
               "Test dxk_dn d2xk_dn2: dxk_dni vs analytical")
         end do
      end do

      ! Second derivative
      do k=1, nc
         do i=1, nc
            do j=1, nc
               if (i==k .and. j==k) then
                  d2xk_dninj = -2 * (n_tot - n(i)) / n_tot**3
               else if (i==k) then
                  d2xk_dninj = (2 * n(i) - n_tot) / n_tot**3
               else if (j==k) then
                  d2xk_dninj = (2 * n(j) - n_tot) / n_tot**3
               else
                  d2xk_dninj = 2 * n(k) / n_tot**3
               end if

               call assert(abs(d2xk_dn(k, i, j) - d2xk_dninj) < 1e-10_pr, &
                  "Test dxk_dn d2xk_dn2: d2xk_dninj vs analytical")
            end do
         end do
      end do
   end subroutine test_dxk_dn

   subroutine test_dx_to_dn()
      use yaeos__math, only: dx_to_dn

      real(pr) :: dns(3)

      dns = dx_to_dn(&
         [0.2_pr, 0.7_pr, 0.1_pr], &
         [1053.0323_pr, -237.758335_pr, -3665.7490_pr]&
         )

      call assert(allclose(dns, [1375.4315_pr, 84.6409421_pr, -3343.3497_pr], 1e-4_pr), &
         "Test dx_to_dn: dns vs reference")
   end subroutine test_dx_to_dn

   subroutine test_sq_error()
      use yaeos__math, only: sq_error


      real(pr) :: sim(5), exps(5), errors_sq(5)

      sim = [0.5_pr, 1.5_pr, 0.4_pr, 1.85_pr, 1.66_pr]
      exps = [0.45_pr, 1.51_pr, 0.42_pr, 1.88_pr, 1.68_pr]

      errors_sq = sq_error(exps, sim)

      call assert(allclose(errors_sq, (sim - exps)**2, 1e10_pr), &
         "Test sq_error: sq_error vs (sim-exps)**2")

   end subroutine test_sq_error

   subroutine test_cardano_method()
      use yaeos__math_linalg, only: pr, cubic_roots
      real(pr) :: p(4)
      real(pr) :: rr(3)
      complex(pr) :: cr(3)
      real(pr) :: num
      integer :: flag

      p = [1, -10, 35, -50]
      call cubic_roots(p, rr, cr, flag)
      call assert(flag == 1, &
         "Test cardano: flag == 1")
      call assert(abs(rr(1) - 5.0_pr) < 1e-10_pr, &
         "Test cardano: roots vs expected")

      p = [0.1, -2.6, 1., 1.]
      call cubic_roots(p, rr, cr, flag)
      call assert(flag == -1, &
         "Test cardano: flag == -1")
      call assert(maxval(abs(rr - [-0.454216, 0.8601986, 25.594016])) < 1e-5_pr, &
         "Test cardano: roots vs expected")

      ! =======================================================================
      ! 3 Real roots
      ! x1 = 1, x2 = 3, x3 = 4
      ! -----------------------------------------------------------------------
      p = [1._pr, -8._pr, 19._pr, -12._pr]

      call cubic_roots(p, rr, cr, flag)

      call assert(maxval(rr - [1.0_pr, 3.0_pr, 4.0_pr]) < 1e-7, &
         "Test cardano: roots vs expected")
      call assert(flag == -1, &
         "Test cardano: flag == -1")

      ! =======================================================================
      ! 1 real root different, 2 equal real roots
      ! x1 = 1, x2 = x3 = 4
      ! -----------------------------------------------------------------------
      p = [1._pr, -9._pr, 24._pr, -16._pr]

      call cubic_roots(p, rr, cr, flag)

      call assert(maxval(rr - [1.0_pr, 4.0_pr, 4.0_pr]) < 1e-7, &
         "Test cardano: roots vs expected")
      call assert(flag == 0, &
         "Test cardano: flag == 0")

      ! =======================================================================
      ! 1 real root, 2 complex
      ! x1 = 1, x2 = i, x3 = -i
      ! -----------------------------------------------------------------------
      p = [1._pr, -1._pr, 1._pr, -1._pr]

      call cubic_roots(p, rr, cr, flag)

      call assert((rr(1) - 1.0_pr) < 1e-7, &
         "Test cardano: roots vs expected")
      call assert(flag == 1, &
         "Test cardano: flag == 1")
   end subroutine test_cardano_method

   subroutine test_rosendo_method()
      use yaeos__math_linalg, only: pr, cubic_roots_rosendo


      real(pr) :: p(4)
      real(pr) :: a, b
      real(pr) :: roots(3)
      complex(pr) :: c_roots(3)
      integer :: flag

      ! =======================================================================
      ! 3 Real roots
      ! x1 = 1, x2 = 3, x3 = 4
      ! -----------------------------------------------------------------------
      p = [1._pr, -8._pr, 19._pr, -12._pr]

      call cubic_roots_rosendo(p, roots, c_roots, flag)

      call assert(maxval(roots - [1.0_pr, 3.0_pr, 4.0_pr]) < 1e-7, &
         "Test rosendo: roots vs expected")
      call assert(flag == -1, &
         "Test rosendo: flag == -1")

      ! =======================================================================
      ! 1 real root different, 2 equal real roots
      ! x1 = 1, x2 = x3 = 4
      ! -----------------------------------------------------------------------
      p = [1._pr, -9._pr, 24._pr, -16._pr]

      call cubic_roots_rosendo(p, roots, c_roots, flag)

      call assert(maxval(roots - [1.0_pr, 4.0_pr, 4.0_pr]) < 1e-7, &
         "Test rosendo: roots vs expected")
      call assert(flag == -1, &
         "Test rosendo: flag == -1")

      ! =======================================================================
      ! 1 real root, 2 complex
      ! x1 = 1, x2 = i, x3 = -i
      ! -----------------------------------------------------------------------
      p = [1._pr, -1._pr, 1._pr, -1._pr]

      call cubic_roots_rosendo(p, roots, c_roots, flag)

      call assert(maxval(roots - [1.0_pr, 1.0_pr, 1.0_pr]) < 1e-7, &
         "Test rosendo: roots vs expected")
      call assert(flag == 1, &
         "Test rosendo: flag == 1")

      ! =======================================================================
      ! 1 real root and two small complex roots
      ! x1 = 1, x2 = a+bi, x3 = a-bi
      ! -----------------------------------------------------------------------
      a = 1._pr
      b = 1e-6_pr

      p = [1._pr, -(2*a + 1), (a**2 + b**2 + 2 * a), -(a**2 + b**2)]

      call cubic_roots_rosendo(p, roots, c_roots, flag)

      call assert(maxval(roots - [1.0_pr, 1.0_pr, 1.0_pr]) < 1e-7, &
         "Test rosendo: roots vs expected")
      call assert(flag == 1, &
         "Test rosendo: flag == 1")
   end subroutine test_rosendo_method

   subroutine test_newton_method()
      use yaeos__math, only: pr, newton
      real(pr) :: x
      real(pr) :: tol=1e-5
      integer :: max_iters = 100
      logical :: failed

      x = 0.5
      call newton(newton_foo, x, tol, max_iters, failed)
      call assert(.not. failed, &
         "Test newton: converged")
      call assert(abs(x - sqrt(2._pr)) < tol, &
         "Test newton: root equals sqrt(2)")
   end subroutine test_newton_method

   subroutine newton_foo(xx, f, df)
      !! Objective function for `test_newton_method`: f(x) = x**2 - 2
      real(pr), intent(in) :: xx
      real(pr), intent(out) :: f
      real(pr), intent(out) :: df
      f = xx**2 - 2
      df = 2*xx
   end subroutine newton_foo
end program test_math
