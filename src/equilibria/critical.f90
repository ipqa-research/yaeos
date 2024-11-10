module yaeos__equilibria_critical
   use yaeos__constants, only: pr
   use yaeos__models, only: ArModel
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState

   implicit none

   type :: CriticalLine
      !! # CriticalLine
      !!
      !! ## Description
      !! This derived type is used to store a critical line between two fluids.
      !! The critical line is calculated using the `critical_line` function. It
      !! uses the continuation method.
      !!
      !! ## Examples
      !! A critical line can be calculated between two fluids using the
      !! `critical_line` function.
      !! In this example we calculate the critical of a binary mixture of
      !! carbon dioxide and bicyclohexyl.
      !!
      !! ```fortran
      !! use yaeos
      !! implicit none
      !! type(CubicEoS) :: model
      !! type(CriticalLine) :: cl
      !! real(pr) :: z0(2), zi(2)
      !!
      !! z0 = [1, 0] ! Pure carbon dioxide
      !! zi = [0, 1] ! Pure bicyclohexyl
      !!
      !! ! Setup the model
      !! tc = [304.21_pr, 727.0_pr]
      !! pc = [73.83_pr, 25.6_pr]
      !! w = [0.223621_pr, 0.427556_pr]
      !! model = PengRobinson76(tc, pc, w)
      !! ! Calculate the critical line
      !! cl = critical_line(model, a0=0.99_pr, z0=z0, zi=zi, dS0=-0.01_pr)
      !! ```
      real(pr), allocatable :: a(:) !! Molar fraction of the second fluid
      real(pr), allocatable :: z0(:) !! Molar fractions of the first fluid
      real(pr), allocatable :: zi(:) !! Molar fractions of the second fluid
      real(pr), allocatable :: P(:) !! Pressure [bar]
      real(pr), allocatable :: V(:) !! Volume [L/mol]
      real(pr), allocatable :: T(:) !! Temperature [K]
   end type CriticalLine

   type, private :: CPSpecs
      integer :: a=1
      integer :: V=2
      integer :: T=3
      integer :: P=4
   end type

   type(CPSpecs), parameter :: CPSpec = CPSpecs()

contains

   type(CriticalLine) function critical_line(model, a0, z0, zi, dS0, max_points, maxP)
      !! # critical_line
      !!
      !! ## Description
      !! Calculates the critical line between two mixtures using the
      !! continuation method. The two mixtures compositions are restricted to
      !! the relation between them, by a parameter \(\alpha\), which represents
      !! the molar fraction of the second fluid with respect to the whole
      !! mixture.
      use yaeos__math_continuation, only: continuation
      use yaeos__math, only: solve_system
      use yaeos__equilibria_equilibrium_state, only: EquilibriumState
      class(ArModel), intent(in) :: model !! Equation of state model
      real(pr), intent(in) :: a0 !! Initial \(\alpha\) value
      real(pr), intent(in) :: z0(:) !! Molar fractions of the first fluid
      real(pr), intent(in) :: zi(:) !! Molar fractions of the second fluid
      real(pr), intent(in) :: dS0 !! Initial step size
      integer, optional, intent(in) :: max_points !! Maximum number of points
      real(pr), optional, intent(in) :: maxP !! Maximum pressure

      real(pr) :: u(size(z0)) !! eigen-vector
      real(pr) :: u_new(size(z0)) !! eigen-vector

      real(pr), allocatable :: XS(:, :) !! Full set of solved points

      real(pr) :: X0(4), T, P, V, z(size(z0))

      integer :: i, j, ns, last_point, npoints

      real(pr) :: max_P

      ! ========================================================================
      ! Handle the input
      ! ------------------------------------------------------------------------

      if (present(max_points)) then
         npoints = max_points
      else
         npoints = 1000
      end if

      if (present(maxP)) then
         max_P = maxP
      else
         max_P = 2500
      end if

      u = (z0 + zi)/sum(z0 + zi)
      z = a0*zi + (1-a0)*z0

      T = sum(model%components%Tc * z)
      P = sum(model%components%Pc * z)
      call model%volume(n=z, P=P, T=T, V=V, root_type="stable")

      X0 = [a0, log([v, T, P])]
      ns = 1

      ! ========================================================================
      ! Calculate the points
      ! ------------------------------------------------------------------------
      XS = continuation(&
         f=foo, X0=X0, ns0=ns, S0=X0(ns), &
         dS0=dS0, max_points=npoints, solver_tol=1e-5_pr, &
         update_specification=update_specification &
         )

      ! Find the last true converged point
      last_point = 0
      do i=1, size(XS, 1)
         if (all(abs(XS(i, :)) < 0.001)) exit
         last_point = i + 1
      end do

      XS = XS(1:last_point-1, :)

      critical_line%z0 = z0
      critical_line%zi = zi
      critical_line%a =  XS(:, 1)
      critical_line%V = exp(XS(:, 2))
      critical_line%T = exp(XS(:, 3))

      allocate(critical_line%P(size(critical_line%a)))
      do i=1, size(critical_line%a)
         z = critical_line%a(i)*zi + (1-critical_line%a(i))*z0

         call model%pressure(&
            n=z, V=critical_line%V(i), T=critical_line%T(i), &
            P=critical_line%P(i))
      end do

   contains

      subroutine foo(X, ns, S, F, dF, dFdS)
         real(pr), intent(in) :: X(:)
         integer, intent(in) :: ns
         real(pr), intent(in) :: S
         real(pr), intent(out) :: F(:)
         real(pr), intent(out) :: dF(:, :)
         real(pr), intent(out) :: dFdS(:)
         real(pr) :: l1

         real(pr) :: z(size(u)), Xsol(3)

         if (X(1) > 1) then
            return
         end if

         F = F_critical(model, X, ns, S, z0, zi, u)
         df = df_critical(model, X, ns, S, z0, zi, u)
         l1 = lambda1(model=model, X=X, s=0.0_pr, z0=z0, zi=zi, u=u, u_new=u_new)
         u = u_new
         dFdS = 0
         dFdS(size(dFdS)) = -1
      end subroutine foo

      subroutine update_specification(X, ns, S, dS, dXdS, iterations)
         real(pr), intent(in out) :: X(:) !! Vector of variables \(X\)
         integer, intent(in out) :: ns !! Position of specified variable
         real(pr), intent(in out) :: S !! Specification variable value
         real(pr), intent(in out) :: dS !! Step of specification in the method
         real(pr), intent(in out) ::  dXdS(:) !! \(\frac{dX}{dS}\)
         integer, intent(in) :: iterations !! Iterations needed to converge point



         ns = maxloc(abs(dXdS), dim=1)
         dS = dXdS(ns)*dS
         dXdS = dXdS/dXdS(ns)
         if (exp(X(4)) > max_P) then
            dS = 0
         end if

      end subroutine update_specification
   end function critical_line

   real(pr) function lambda1(model, X, s, z0, zi, u, u_new, P)
      !! # lambda1
      !!
      !! Calculation of the first restriction of a critical point
      !!
      !! \[
      !!  \lambda_1(s) = \frac{d^2tpd}{ds^2} = 0
      !! \]
      use stdlib_linalg, only: eigh, linalg_state_type
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z0(:) !! Molar fractions of the first fluid
      real(pr), intent(in) :: zi(:) !! Molar fractions of the second fluid
      real(pr), intent(in) :: s
      !! Distance between the two fluids compositions to the main composition
      real(pr), intent(in) :: X(4) !! Vector of variables
      real(pr), intent(in) :: u(:)
      !! Eigen-vector that defines the direction between the two compositions
      !! \[ n_i = z_i + s  u_i \sqrt{z_i} \]
      real(pr), optional, intent(out) :: u_new(:)
      !! Eigen-vector corresponding to the smallest eigenvalue of the matrix
      !! \[ M_{ij} = \sqrt{z_i z_j} \frac{\partial \ln f_i}{\partial n_j} \]
      real(pr), optional, intent(out) :: P !! Pressure of the system [bar]

      real(pr) :: n(size(z0)), V, T
      real(pr) :: dlnf_dn(size(z0), size(z0))
      real(pr) :: lambda(size(z0)), vectors(size(z0), size(z0))

      type(linalg_state_type) :: stat

      integer :: i, j, nc
      real(pr) :: M(size(z0), size(z0)), z(size(z0)), Pin

      nc = size(z0)

      z = X(1) * zi + (1-X(1))*z0
      n = z + s * u * sqrt(z)
      V = exp(X(2))
      T = exp(X(3))

      call model%lnfug_vt(n=n, V=V, T=T, dlnfdn=dlnf_dn, P=Pin)

      do i=1,nc
         do j=1,nc
            M(i, j) = sqrt(z(i)*z(j)) * dlnf_dn(i, j)
         end do
      end do

      call eigh(A=M, lambda=lambda, vectors=vectors, err=stat)
      if (.not. stat%ok()) then
         write(*, *) stat%print_msg()
      end if

      lambda1 = lambda(minloc(abs(lambda), dim=1))
      if (present(u_new)) u_new = vectors(:, minloc(abs(lambda), dim=1))
      if (present(P)) P = Pin
   end function lambda1

   function F_critical(model, X, ns, S, z0, zi, u)
      !! # F_critical
      !!
      !! ## Description
      !! Function that should be equal to zero at a critical point is found.
      !! The second criticality condition is calculated as a numerical
      !! derivative with `eps=1e-10`.
      !!
      !! \[
      !! F = \begin{bmatrix}
      !!   \lambda_1(s) \\
      !!   \frac{\partial \lambda_1(s+\epsilon) - \lambda_1(s-\epsilon)}{2\epsilon} \\
      !!   \ln P - X_4 \\
      !!   X_{ns} - S
      !! \end{bmatrix} = 0
      !! \]
      class(ArModel), intent(in) :: model !! Equation of state model
      real(pr), intent(in) :: X(4) !! Vector of variables
      integer, intent(in) :: ns !! Position of the specification variable
      real(pr), intent(in) :: S !! Specification variable value
      real(pr), intent(in) :: z0(:) !! Molar fractions of the first fluid
      real(pr), intent(in) :: zi(:) !! Molar fractions of the second fluid
      real(pr), intent(in) :: u(:) !! Eigen-vector

      real(pr) :: F_critical(4)
      real(pr) :: z(size(u))

      real(pr) :: V, T, P

      real(pr), parameter :: eps=1e-10_pr

      V = exp(X(2))
      T = exp(X(3))
      z = X(1) * zi + (1-X(1)) * z0

      if(any(z < 0) ) return

      F_critical(1) = lambda1(model=model, X=X, s=0.0_pr, z0=z0, zi=zi, u=u, P=P)
      F_critical(2) = (&
         lambda1(model=model, X=X, s= eps, zi=zi, z0=z0, u=u) &
         - lambda1(model=model, X=X, s=-eps, zi=zi, z0=z0, u=u))/(2*eps)
      F_critical(3) = log(P) - X(4)
      F_critical(4) = X(ns) - S
   end function F_critical

   function df_critical(model, X, ns, S, z0, zi, u)
      !! # df_critical
      !!
      !! ## Description
      !! Calculates the Jacobian of the critical point function `F_critical`.
      class(ArModel), intent(in) :: model !! Equation of state model
      real(pr), intent(in) :: X(4) !! Vector of variables
      integer, intent(in) :: ns !! Position of the specification variable
      real(pr), intent(in) :: S !! Specification variable value
      real(pr), intent(in) :: z0(:) !! Molar fractions of the first fluid
      real(pr), intent(in) :: zi(:) !! Molar fractions of the second fluid
      real(pr), intent(in) :: u(:) !! Eigen-vector
      real(pr) :: df_critical(4, 4) !! Jacobian of the critical point function

      real(pr), parameter :: eps=1e-5_pr

      real(pr) :: dx(4), F1(4), F2(4)

      integer :: i

      df_critical = 0
      do i=1,4
         dx = 0
         dx(i) = eps
         F2 = F_critical(model, X+dx, ns, S, z0, zi, u)
         F1 = F_critical(model, X-dx, ns, S, z0, zi, u)
         df_critical(:, i) = (F2 - F1)/(2*eps)
      end do
   end function df_critical

   type(EquilibriumState) function critical_point(&
      model, z0, zi, spec, S, max_iters, u0, V0, T0, a0 &
      )
      !! # critical_point
      !!
      !! ## Description
      !! Calculates a single critical point of a mixture using a Newton-Raphson
      !! method. It is possible to specify different variables to be fixed with
      !! the `spec` argument. 
      use yaeos__math, only: solve_system
      class(ArModel), intent(in) :: model !! Equation of state model
      real(pr), intent(in) :: z0(:) !! Molar fractions of the first fluid
      real(pr), intent(in) :: zi(:) !! Molar fractions of the second fluid
      integer, intent(in) :: spec !! Specification `[1:"z", 2:"V", 3:"T", 4:"P"]`
      real(pr), intent(in) :: S !! Specification value
      integer, intent(in) :: max_iters !! Maxiumum number of iterations

      real(pr), optional, intent(in) :: V0 !! Initial volume [L/mol].
      real(pr), optional, intent(in) :: T0 !! Initial temperature [K].
      real(pr), optional, intent(in) :: a0 !! Initial \(\alpha\) value
      real(pr), optional, intent(in) :: u0(:) !! Initial eigen-vector

      real(pr) :: X(4)
      integer :: ns
      real(pr) :: F(4), df(4, 4), dX(4), u(size(z0))
      real(pr) :: V, T, P

      real(pr) :: z(size(z0)), u_new(size(z0)), l
      integer :: i

      ! ========================================================================
      ! Handle the input
      ! ------------------------------------------------------------------------
      if (present(a0)) then
         X(1) = a0
      else
         X(1) = 0.0_pr
      end if

      z = X(1)*zi + (1-X(1))*z0

      if (present(u0)) then
         u = u0
      else
         u = (z0 + zi)/sum(z0 + zi)
      end if

      if (present(T0)) then
         X(3) = log(T0)
      else
         X(3) = log(sum(model%components%Tc * z))
      end if

      if (present(V0)) then
         X(2) = log(V0)
      else
         call model%volume(&
            n=z, P=sum(model%components%Pc * z), T=exp(X(3)), V=X(2), root_type="stable")
         X(2) = log(X(2))
      end if

      X(4) = log(sum(model%components%Pc * z))
      ns = spec

      ! ========================================================================
      ! Solve the system of equations
      ! ------------------------------------------------------------------------
      do i=1,max_iters
         F = F_critical(model, X, ns, S, z0, zi, u)
         df = df_critical(model, X, ns, S, z0, zi, u)
         dX = solve_system(A=df, b=-F)

         do while(maxval(abs(dX/X)) > 1e-1)
            dX = dX/10
         end do

         if (maxval(abs(X)) < 1e-5) exit

         X = X + dX
         l = lambda1(model, X, 0.0_pr, z0, zi, u, u_new)
         u = u_new

         critical_point%iters = i
      end do

      critical_point%x = z
      critical_point%y = z
      critical_point%beta = 0
      critical_point%Vx = exp(X(2))
      critical_point%Vy = exp(X(2))
      critical_point%T  = exp(X(3))
      call model%pressure(n=z, V=critical_point%Vx, T=critical_point%T, P=critical_point%P)
      critical_point%kind = "critical"

   end function critical_point
end module yaeos__equilibria_critical
