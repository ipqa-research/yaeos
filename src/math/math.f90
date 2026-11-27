module yaeos__math
   !! # Mathematical methods for `yaeos`
   !!
   !! # Description
   !! This module provides all the relevant mathematical functions used in this
   !! library. Most important ones are:
   !!
   !! - newton: Newton solving method
   !! - solve_system: Solving linear system Ax = b
   !! - continuation: Continuation method for line tracing
   !!
   !! # Examples
   !!
   !! ## Squared error calculation
   !! ```fortran
   !!  use yaeos__constants, only: pr
   !!  use yaeos__math, only: sq_error
   !!  print *, sq_error(2.5_pr, 3.0_pr)
   !! ```
   !!
   !! ```fortran
   !!  use yaeos__constants, only: pr
   !!  use yaeos__math, only: sq_error
   !!  real(pr) :: x(2) = [2.5_pr, 5.0_pr], y(2) = [3.0_pr, 4.5_pr]
   !!  ! It also works with arrays
   !!  print *, sq_error(x, y)
   !! ```

   use yaeos__math_continuation, only: continuation
   use yaeos__math_linalg, only: solve_system, cubic_roots
   use yaeos__constants, only: pr

   implicit none

   type :: Point
      !! # Point
      !! Represents the single point of a line segment.
      !!
      !! # Description
      !! Representation of a point in the 2D space with its coordinates. It is
      !! used to represent either the intersection of two lines segments or the
      !! self-intersection of a single line. In the first case the `i` and `j`
      !! attributes represent the indices of the line segments that intersected.
      !! In the second case both represent the index of the line segment that
      !! self-intersected.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! ```
      real(pr) :: x !! X coordinate
      real(pr) :: y !! Y coordinate
      integer :: i !! Index of the first line segment
      integer :: j !! Index of the second line segment
   end type Point

   abstract interface
      subroutine f_1d(x, f, df)
         !! Callback interface for a scalar residual and its derivative.
         import pr
         real(pr), intent(in) :: x !! Current scalar iterate.
         real(pr), intent(out) :: f !! Residual value \(f(x)\).
         real(pr), intent(out) :: df !! Derivative \(df/dx\).
      end subroutine f_1d
   end interface

   interface newton
      module procedure :: newton_1d
   end interface newton

contains
   elemental real(pr) function sq_error(exp, pred)
      !! # Squared error between two values.
      !!
      !! # Description
      !! Returns the squared relative error \(((exp-pred)/exp)^2\). `exp`
      !! must be nonzero; the result is dimensionless and elementally
      !! applicable to scalar or array values.
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  error = sq_error(true_value, model_value)
      !! ```
      use yaeos__constants, only: pr
      real(pr), intent(in) :: exp !! Reference or expected value; must be nonzero.
      real(pr), intent(in) :: pred !! Predicted value being compared with `exp`.
      sq_error = ((exp - pred)/exp)**2
   end function sq_error

   function central_numdiff(f, x, eps) result(df)
      !! Approximate the Jacobian of a vector function by centered differences.
      !! Column `i` is evaluated using a perturbation \(dx_i=x_i\,eps\); thus
      !! this routine requires nonzero entries in `x` for finite differences.
      interface
         function f(x)
            !! Vector-valued function differentiated by `central_numdiff`.
            import pr
            real(pr), intent(in) :: x(:) !! Function input vector.
            real(pr) :: f(size(x)) !! Function output vector.
         end function f
      end interface
      real(pr), intent(in) :: x(:) !! Evaluation point.
      real(pr), intent(in) :: eps !! Relative finite-difference step multiplier.
      real(pr) :: df(size(x), size(x)) !! Approximate Jacobian matrix.

      integer :: i !! Loop/index variable.
      real(pr) :: f1(size(x))
      real(pr) :: f2(size(x))
      real(pr) :: dx(size(x))

      do i=1,size(x)
         dx = 0
         dx(i) = x(i) * eps
         df(:, i) = (f(x + dx) - f(x - dx))/(2*dx(i))
      end do
   end function central_numdiff

   function dx_to_dn(x, dx) result(dn)
      !! # dx_to_dn
      !!
      !! # Description
      !! Convert the mole fraction derivatives of a quantity (calculated
      !! so they do not sum to 1) to mole number derivatives (where the mole
      !! fractions do sum to one). Requires the derivatives and the mole fractions
      !! of the mixture.
      !! From [Chemicals (Python)](https://chemicals.readthedocs.io/chemicals.utils.html?highlight=dxs_to_dns#chemicals.utils.dxs_to_dns)
      use yaeos__constants, only: pr
      real(pr), intent(in) :: x(:) !! Mole-fraction vector.
      real(pr), intent(in) :: dx(:) !! Derivative of the property with respect to mole fractions.
      real(pr) :: dn(size(x))

      real(pr) :: sum_xdx

      dn = 0

      sum_xdx = sum(x * dx)

      dn = dx - sum_xdx
   end function dx_to_dn

   function derivative_dxk_dni(n) result(dxk_dni)
      !! # derivative_dxk_dni
      !!
      !! # Description
      !! Calculate the mole fraction first derivatives respect to mole numbers
      !!
      real(pr), intent(in) :: n(:) !! Component mole amounts.
      real(pr) :: dxk_dni(size(n), size(n)) !! \(\partial x_k/\partial n_i\), indexed `(k,i)`.

      real(pr) :: n_tot !! Total number of moles.
      integer :: nc !! Number of components.
      integer :: k
      integer :: i !! Loop/index variable.

      n_tot = sum(n)
      nc = size(n)

      dxk_dni = 0.0_pr
      do concurrent(k=1:nc, i=1:nc)
         if (k == i) then
            dxk_dni(k,i) = (n_tot - n(i)) / n_tot**2
         else
            dxk_dni(k,i) = -n(k) / n_tot**2
         end if
      end do
   end function derivative_dxk_dni

   function derivative_d2xk_dnidnj(n) result(d2xk_dnidnj)
      !! # derivative_d2xk_dnidnj
      !!
      !! # Description
      !! Calculate the mole fraction second derivatives respect to mole numbers
      !!
      real(pr), intent(in) :: n(:) !! Component mole amounts.
      real(pr) :: d2xk_dnidnj(size(n), size(n), size(n)) !! \(\partial^2x_k/\partial n_i\partial n_j\), indexed `(k,i,j)`.

      real(pr) :: n_tot !! Total number of moles.
      integer :: nc !! Number of components.
      integer :: k !! Loop index `k`.
      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.

      n_tot = sum(n)
      nc = size(n)

      d2xk_dnidnj = 0.0_pr
      do concurrent (k=1:nc, i=1:nc, j=1:nc)
         if (i==k .and. j==k) then
            d2xk_dnidnj(k,i,j) = -2 * (n_tot - n(i)) / n_tot**3
         else if (i==k) then
            d2xk_dnidnj(k,i,j) = (2 * n(i) - n_tot) / n_tot**3
         else if (j==k) then
            d2xk_dnidnj(k,i,j) = (2 * n(j) - n_tot) / n_tot**3
         else
            d2xk_dnidnj(k,i,j) = 2 * n(k) / n_tot**3
         end if
      end do
   end function derivative_d2xk_dnidnj

   subroutine newton_1d(f, x, tol, max_iters, failed)
      !! # newton_1d
      !!
      !! # Description
      !! Find a root of a scalar function by damped Newton iteration.
      !! Each step is halved until its magnitude is at most half the current
      !! iterate's magnitude. This is a step-size safeguard, not a residual
      !! line search.
      procedure(f_1d) :: f !! Function and derivative callback.
      real(pr), intent(in out) :: x !! Initial guess, replaced by the root estimate.
      real(pr), intent(in) :: tol !! Absolute tolerance for residual and step size.
      integer, intent(in) :: max_iters !! Maximum Newton iterations.
      logical, intent(out) :: failed !! True when the iteration limit is exhausted.

      integer :: i !! Loop/index variable.
      real(pr) :: fval
      real(pr) :: df
      real(pr) :: step


      fval = 10
      step = 10

      do i=1, max_iters
         if (abs(fval) < tol .or. abs(step) < tol)  exit
         call f(x, fval, df)


         step = fval/df
         do while (abs(step) > 0.5 * abs(x))
            step = step/2
         end do

         x = x - step

         failed = i >= max_iters
      end do

   end subroutine newton_1d

   elemental function interpol(x1, x2, y1, y2, x_obj) result(y)
      !! # interpol
      !!
      !! # Description
      !! Linear interpolation.
      !!
      !! Calculates the linear interpolation between two points at a desired
      !! x value with the equation:
      !! \[
      !!    y = \frac{y_2 - y_1}{x_2 - x_1} \cdot (x_{obj})  - x_1 + y_1
      !! \]
      !!
      !! Since this function is defined as `elemental` it will also interpolate
      !! a set of vectors.
      !!
      !! Examples of usage:
      !!
      !! ```fortran
      !! x1 = 2
      !! x2 = 5
      !! y1 = 2
      !! y2 = 9
      !! y = interpol(x1, x2, y1, y2, 2.3)
      !! ```
      !!
      !! ```fortran
      !! x1 = 2
      !! x2 = 5
      !! y1 = [2, 6]
      !! y2 = [9, 15]
      !! y = interpol(x1, x2, y1, y2, 2.3)
      !! ```
      real(pr), intent(in) :: x1 !! First point x value
      real(pr), intent(in) :: x2 !! Second point x value
      real(pr), intent(in) :: y1 !! First point y value
      real(pr), intent(in) :: y2 !! Second point y value
      real(pr), intent(in) :: x_obj !! Desired x value to interpolate
      real(pr) :: y !! y value at `x_obj`
      y = (y2 - y1)/(x2 - x1)*(x_obj - x1) + y1
   end function interpol

   function intersect_one_line(lx, ly) result(intersections)
      !! # intersect_one_line
      !! Find the intersections of a single line with itself.
      !!
      !! # Description
      !! This function finds the self-intersections in a line. This is
      !! determined by checking all possible pairs of lines segments.
      !! The iteration starts from the first segment of the line and compares
      !! it with all subsequent segments to find intersections. Then it goes
      !! to the next segment and repeats the process. The intersections are
      !! stored in an array of `point` type, which contains the coordinates
      !! of the intersection points.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! ```
      !!
      !! # References
      !!
      real(pr), intent(in) :: lx(:) !! X coordinates of the polyline vertices.
      real(pr), intent(in) :: ly(:) !! Y coordinates of the polyline vertices.
      type(point), allocatable :: intersections(:) !! Detected non-adjacent segment intersections.

      real(pr) :: s
      real(pr) :: t
      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.

      real(pr) :: x
      real(pr) :: y
      real(pr) :: xold
      real(pr) :: yold

      xold = 9999
      yold = 9999

      allocate (intersections(0))
      line1: do i = 2, size(lx) - 1
         line2: do j = i + 2, size(lx)
            associate ( &
               x1 => lx(i - 1), x2 => lx(i), &
               x3 => lx(j), x4 => lx(j - 1), &
               y1 => ly(i - 1), y2 => ly(i), &
               y3 => ly(j), y4 => ly(j - 1))

               call intersects(x1, x2, x3, x4, y1, y2, y3, y4, s, t)
               if (0 <= s .and. s <= 1 .and. 0 <= t .and. t <= 1) then
                  x = s*(x2 - x1) + x1
                  y = s*(y2 - y1) + y1
                  if (abs(x - xold) > 1 .or. abs(y - yold) > 1) then
                     xold = x
                     yold = y
                     ! Use earliest point for the "other" line
                     intersections = [intersections, point(x, y, i, j - 1)]
                  end if
               end if
            end associate
         end do line2
      end do line1
      if (size(intersections) > 3) then
         deallocate(intersections)
         allocate(intersections(0))
      end if
   end function intersect_one_line

   subroutine intersects(x1, x2, x3, x4, y1, y2, y3, y4, s, t)
      !! # intersects
      !! Calculate the intersection between two line segments.
      !!
      !! # Description
      !! This subroutine calculates the intersection point of two line segments
      !! defined by their endpoints (x1, y1) to (x2, y2) and (x3, y3) to (x4, y4).
      !! If the segments intersect, the parameters s and t will contain the
      !! normalized distances along each segment to the intersection point.
      !! The intersection point can be calculated as:
      !! \[
      !!    x = s \cdot (x2 - x1) + x1
      !!    y = s \cdot (y2 - y1) + y1
      !! \]
      !! If the segments do not intersect, s and t will be outside
      !! the range [0, 1].
      !!
      !! # Examples
      !!
      !! ```fortran
      !! ```
      !!
      !! # References
      !!
      real(pr), intent(in) :: x1 !! X coordinate of the first endpoint of segment 1.
      real(pr), intent(in) :: x2 !! X coordinate of the second endpoint of segment 1.
      real(pr), intent(in) :: x3 !! X coordinate of the first endpoint of segment 2.
      real(pr), intent(in) :: x4 !! X coordinate of the second endpoint of segment 2.
      real(pr), intent(in) :: y1 !! Y coordinate of the first endpoint of segment 1.
      real(pr), intent(in) :: y2 !! Y coordinate of the second endpoint of segment 1.
      real(pr), intent(in) :: y3 !! Y coordinate of the first endpoint of segment 2.
      real(pr), intent(in) :: y4 !! Y coordinate of the second endpoint of segment 2.
      real(pr), intent(out) :: s !! Normalized position along segment 1; `s` and `t` both in [0,1] iff the segments intersect.
      real(pr), intent(out) :: t !! Normalized position along segment 2.

      real(pr) :: A(2, 2)
      real(pr) :: b(2)
      real(pr) :: tmp !! Temporary intermediate value.

      A(1, :) = [x2 - x1, x3 - x4]
      A(2, :) = [y2 - y1, y3 - y4]
      b = [x3 - x1, y3 - y1]

      b = solve_system(a, b)
      s = b(1)
      t = b(2)
   end subroutine intersects

   subroutine levenberg_marquardt(&
      fun, tol, X, F, info &
      )
      !! Solve a nonlinear least-squares system with MINPACK's `lmdif1`.
      !! The routine updates `X` and returns residuals and MINPACK status.
      use minpack_module, only: lmdif1, lmdif
      use yaeos__auxiliar, only: optval
      use yaeos__math_linalg, only: solve_system

      interface
         subroutine fun(m, n, x, fvec, iflag)
            !! # fun
            !!
            !! # Description
            !! MINPACK residual callback for rectangular least-squares systems.
            import pr
            integer, intent(in) :: m !! Number of residual equations.
            integer, intent(in) :: n !! Number of unknowns.
            real(pr), intent(in) :: x(n) !! Current unknown vector.
            real(pr), intent(out) :: fvec(m) !! Residual vector.
            integer, intent(in out) :: iflag !! MINPACK callback control/status flag.
         end subroutine fun
      end interface

      real(pr), intent(in) :: tol !! Convergence tolerance passed to MINPACK.
      real(pr), intent(in out) :: X(:) !! Variables vector
      real(pr), intent(out) :: F(:) !! Function values at solved point
      integer, intent(in out) :: info !! MINPACK termination code.

      integer :: m
      integer :: n
      integer :: iwa(size(x))
      integer :: lwa
      real(pr) :: wa(size(F) * size(x) + 5*size(x) + size(f))

      m = size(F)
      n = size(X)
      lwa = size(F) * size(X) + 5 * size(X) + size(F) + 10
      call lmdif1(fun, m, n, X, F, tol, Info, Iwa, Wa, Lwa)
      ! call lmdif(fun, m, n, x, F, tol, tol, Gtol, 1000, Epsfcn=1e-5_pr, Diag, Mode, Factor, Nprint, Info, Nfev, Fjac, Ldfjac, Ipvt, Qtf, Wa1, Wa2, Wa3, Wa4)

   end subroutine levenberg_marquardt

   subroutine powel_hybrid(&
      fun, tol, X, F, info &
      )
      !! # powel_hybrid
      !!
      !! # Description
      !! Solve a square nonlinear system with MINPACK's Powell hybrid method.
      use minpack_module, only: hybrd1
      use yaeos__auxiliar, only: optval
      use yaeos__math_linalg, only: solve_system

      interface
         subroutine fun(n, x, fvec, iflag)
            !! # fun
            !!
            !! # Description
            !! MINPACK residual callback for a square nonlinear system.
            import pr
            integer, intent(in) :: n !! Number of unknowns and residuals.
            real(pr), intent(in) :: x(n) !! Current unknown vector.
            real(pr), intent(out) :: fvec(n) !! Residual vector.
            integer, intent(in out) :: iflag !! MINPACK callback control/status flag.
         end subroutine fun
      end interface

      real(pr), intent(in) :: tol !! Convergence tolerance passed to MINPACK.
      real(pr), intent(in out) :: X(:) !! Variables vector
      real(pr), intent(out) :: F(:) !! Function values at solved point
      integer, intent(out) :: info !! MINPACK termination code.

      integer :: m
      integer :: n
      integer :: iwa(size(x))
      integer :: lwa
      real(pr) :: wa((size(X)*(3*size(X)+13))/2 + 5)

      m = size(F)
      n = size(X)
      lwa = (n*(3*n+13))/2 + 5
      call hybrd1(fun, n, x, F, tol, info, Wa, Lwa)
   end subroutine powel_hybrid

end module yaeos__math
