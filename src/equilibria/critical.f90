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
      integer, allocatable :: ns(:) !! Specified variable
      integer, allocatable :: iters(:) !! Iterations needed for this point
   end type CriticalLine

   type, private :: CPSpecs
      !! Enumerator to handle the possible specifications for a critical point.
      integer :: a=1 !! Specify \( \alpha \)
      integer :: V=2 !! Specify \( V \)
      integer :: T=3 !! Specify \( T \)
      integer :: P=4 !! Specify \( P \)
   end type CPSpecs

   type(CPSpecs), parameter :: spec_CP = CPSpecs()
   !! Specification variables for a critical point or critical line
   !! calculation.

contains

   type(CriticalLine) function critical_line(&
      model, a0, z0, zi, ns0, S0, dS0, max_points, maxP, first_point &
      )
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
      integer, intent(in) :: ns0 !! Position of the specification variable
      real(pr), intent(in) :: S0 !! Specified value
      real(pr), intent(in) :: dS0 !! Initial step size
      integer, optional, intent(in) :: max_points !! Maximum number of points
      real(pr), optional, intent(in) :: maxP !! Maximum pressure
      type(EquilibriumState), optional, intent(in) :: first_point


      real(pr) :: u(size(z0)) !! eigen-vector
      real(pr) :: u_new(size(z0)) !! eigen-vector

      real(pr), allocatable :: XS_i(:), XS(:, :)!! Full set of solved points

      integer :: ns
      real(pr) :: S

      real(pr) :: X0(4), T, P, V, z(size(z0))

      integer :: i, npoints

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
      call model%volume(n=z, P=P, T=T, V=V, root_type="vapor")

      X0 = [a0, log([v, T, P])]

      if (present(first_point)) then
         X0 = [&
            first_point%x(2), &
            log([first_point%Vx, first_point%T, first_point%P])]
      end if

      X0(ns0) = S0
      ns = ns0

      ! ========================================================================
      ! Calculate the points
      ! ------------------------------------------------------------------------
      allocate(critical_line%ns(0), critical_line%iters(0))
      allocate(critical_line%P(0), critical_line%T(0), critical_line%V(0), critical_line%a(0))

      solve_points: block
         use yaeos__math, only: solve_system
         use yaeos__math_continuation, only: full_newton
         real(pr) :: X(4), dX(4), dS, F(4), dF(4,4), dFdS(4), dXdS(4)
         real(pr) :: u_new(size(z0)), l1, Si

         integer :: its, real_its

         X = X0
         dS = dS0
         S = X(ns)

         do i=1,npoints
            dX = 1
            F = 10
            X0 = X
            Si = X0(ns)
            its = 0
            real_its = 0

            do while(&
               maxval(abs(dX)) > 1e-4 &
               .and. maxval(abs(F)) > 1e-5 &
               .and. real_its < 2000)

               its = its + 1
               real_its = real_its + 1

               F = F_critical(model, X, ns, Si, z0, zi, u)
               dF = df_critical(model, X, ns, Si, z0, zi, u)
               l1 = lambda1(model=model, X=X, s=0.0_pr, z0=z0, zi=zi, u=u, u_new=u_new)
               dX = solve_system(dF, -F)
               u = u_new

               do while(abs(dX(1)/X(1)) > 0.1 .and. i < 5)
                  dX = dX/2
               end do

               do while(abs(maxval(dX(:)/X(:))) > 0.1)
                  dX = dX/2
               end do

               X = X + dX
            end do

            ! ==============================================================
            ! Cases where the line must be stopped
            ! --------------------------------------------------------------
            if (real_its == 2000) exit
            if (any(isnan(X))) exit
            if (exp(X(spec_CP%P)) > max_P) exit


            ! ==============================================================
            ! Save point
            ! --------------------------------------------------------------
            critical_line%a = [critical_line%a, X(1)]
            critical_line%V = [critical_line%V, exp(X(2))]
            critical_line%T = [critical_line%T, exp(X(3))]
            critical_line%P = [critical_line%P, exp(X(4))]
            critical_line%ns = [critical_line%ns, ns]
            critical_line%iters = [critical_line%iters, its]

            ! ==============================================================
            ! Determination of new specification
            ! --------------------------------------------------------------
            dFdS = [0, 0, 0, -1]
            dXdS = solve_system(dF, -dFdS)
            ns = maxloc(abs(dXdS), dim=1)
            dS = dXdS(ns)*dS
            dXdS = dXdS/dXdS(ns)

            ! ==============================================================
            ! Avoid big steps in pressure
            ! --------------------------------------------------------------
            ! do while(abs(exp(X(4) + dXdS(4) * dS) - exp(X(4))) > 10)
            !    dS = dS * 0.99
            ! end do

            ! Next step
            X = X + dXdS*dS
         end do
      end block solve_points

      critical_line%z0 = z0
      critical_line%zi = zi
   end function critical_line

   real(pr) function lambda1(model, X, s, z0, zi, u, u_new, P)
      !! # lambda1
      !!
      !! Calculation of the first restriction of a critical point
      !!
      !! \[
      !!  \lambda_1(s) = \frac{d^2tpd}{ds^2} = 0
      !! \]
      use yaeos__math_linalg, only: eigen
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

      ! type(linalg_state_type) :: stat

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

      call eigen(A=M, eigenvalues=lambda, eigenvectors=vectors)

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

      real(pr), parameter :: eps=1e-5_pr

      V = exp(X(2))
      T = exp(X(3))
      z = X(1) * zi + (1-X(1)) * z0

      ! if(any(z < 0) ) return

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

      real(pr) :: eps

      real(pr) :: dx(4), F1(4), F2(4)

      integer :: i

      if (any(X(1)*zi + (1-X(1))*z0 > 0.99)) then
         eps = 1e-3_pr
      else
         eps = 1e-6_pr
      end if

      df_critical = 0
      do i=1,4
         dx = 0
         dx(i) = eps
         F2 = F_critical(model, X+dx, ns, S, z0, zi, u)
         F1 = F_critical(model, X-dx, ns, S, z0, zi, u)
         df_critical(:, i) = (F2 - F1)/(2*eps)
      end do
   end function df_critical

   function F_cep(model, X, ns, S, z0, zi, u, u_new)
      class(ArModel), intent(in) :: model !! Equation of state model
      real(pr), intent(in) :: z0(:) !! Molar fractions of the first fluid
      real(pr), intent(in) :: X(size(z0) + 4) !! Vector of variables
      integer, intent(in) :: ns !! Position of the specification variable
      real(pr), intent(in) :: S !! Specification variable value
      real(pr), intent(in) :: zi(:) !! Molar fractions of the second fluid
      real(pr), intent(in) :: u(:) !! Eigen-vector
      real(pr), optional, intent(out) :: u_new(:) !! updated

      real(pr) :: F_cep(size(z0)+ 4)
      real(pr) :: z(size(u))

      real(pr) :: V, T, P
      real(pr) :: Xcp(4)

      real(pr) :: Vc, Pc, lnf_z(size(z0))
      real(pr) :: y(size(z0)), Vy, Py, lnf_y(size(z0))

      real(pr), parameter :: eps=1e-5_pr

      integer :: nc

      nc = size(z0)

      y = X(:nc)
      Vy = exp(X(nc+1))
      Vc = exp(X(nc+2))
      T = exp(X(nc+3))
      z = X(nc+4) * zi + (1-X(nc+4)) * z0

      if(any(z < 0) ) return

      call model%lnfug_vt(n=y, V=Vy, T=T, P=Py, lnf=lnf_y)
      call model%lnfug_vt(n=z, V=Vc, T=T, P=Pc, lnf=lnf_z)
      Xcp(1) = X(nc+4)
      Xcp(2) = log(Vc)
      Xcp(3) = log(T)
      Xcp(4) = log(Pc)

      F_cep(1) = lambda1(model=model, X=Xcp, s=0.0_pr, z0=z0, zi=zi, u=u, P=Pc, u_new=u_new)
      F_cep(2) = (&
         lambda1(model=model, X=Xcp, s= eps, zi=zi, z0=z0, u=u) &
         - lambda1(model=model, X=Xcp, s=-eps, zi=zi, z0=z0, u=u))/(2*eps)
      F_cep(3) = log(Pc) - log(Py)
      F_cep(4:nc+3) = lnf_y - lnf_z
      F_cep(nc+4) = sum(y) - 1
   end function F_cep

   function df_cep(model, X, ns, S, z0, zi, u)
      !! # df_critical
      !!
      !! ## Description
      !! Calculates the Jacobian of the critical point function `F_critical`.
      class(ArModel), intent(in) :: model !! Equation of state model
      real(pr), intent(in) :: z0(:) !! Molar fractions of the first fluid
      real(pr), intent(in) :: X(size(z0)+4) !! Vector of variables
      integer, intent(in) :: ns !! Position of the specification variable
      real(pr), intent(in) :: S !! Specification variable value
      real(pr), intent(in) :: zi(:) !! Molar fractions of the second fluid
      real(pr), intent(in) :: u(:) !! Eigen-vector
      real(pr) :: df_cep(size(z0)+4, size(z0)+4) !! Jacobian of the critical point function

      real(pr) :: eps

      real(pr) :: dx(size(z0)+4), F1(size(z0)+4), F2(size(z0)+4)

      integer :: i

      if (any(X(1)*zi + (1-X(1))*z0 > 0.99)) then
         eps = 1e-3_pr
      else
         eps = 1e-6_pr
      end if

      df_cep = 0
      do i=1,size(df_cep, 1)
         dx = 0
         dx(i) = eps
         F2 = F_cep(model, X+dx, ns, S, z0, zi, u)
         F1 = F_cep(model, X-dx, ns, S, z0, zi, u)
         df_cep(:, i) = (F2 - F1)/(2*eps)
      end do
   end function df_cep

   type(EquilibriumState) function critical_point(&
      model, z0, zi, spec, S, max_iters, u0, V0, T0, a0, P0 &
      )
      !! # critical_point
      !!
      !! ## Description
      !! Calculates a single critical point of a mixture using a Newton-Raphson
      !! method. It is possible to specify different variables to be fixed with
      !! the `spec` argument, the `spec_CP` variable helps when selecting the
      !! specified variable.
      !!
      !! ## Examples
      !!
      !! ### Default behaviour
      !!
      !! ```fortran
      !!   cp = critical_point(&
      !!        model, z0, zi, S=0.5_pr, spec=spec_CP%a, max_iters=1000)
      !! ```
      !!
      !! ### Specifiying another variable
      !! The natural variables are a, lnV, lnT and lnP. So it is important to
      !! specify the variable in logaritmic scale if that is the case.
      !!
      !! ```fortran
      !!   cp = critical_point(model, z0, zi, S=log(200._pr), spec=spec_CP%P, max_iters=1000)
      !! ```
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
      real(pr), optional, intent(in) :: P0 !! Initial Pressure [bar]
      real(pr), optional, intent(in) :: u0(:) !! Initial eigen-vector

      real(pr) :: X(4)
      integer :: ns
      real(pr) :: F(4), df(4, 4), dX(4), u(size(z0))

      real(pr) :: z(size(z0)), u_new(size(z0)), l
      integer :: i

      ! ========================================================================
      ! Handle the input
      ! ------------------------------------------------------------------------
      if (present(a0)) then
         X(1) = a0
      else if (spec == spec_CP%a) then
         X(1) = S
      else
         X(1) = 0.5_pr
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

      if (ns == spec_CP%P) then
         X(4) = S
      else if (present(P0)) then
         X(4) = log(P0)
      else
         X(4) = log(sum(model%components%Pc * z))
      end if


      if (present(V0)) then
         X(2) = log(V0)
      else
         call model%volume(n=z, P=exp(X(4)), T=exp(X(3)), V=X(2), root_type="stable")

         X(2) = log(X(2))
      end if

      ns = spec
      X(ns) = S

      ! ========================================================================
      ! Solve the system of equations
      ! ------------------------------------------------------------------------
      do i=1,max_iters
         F = F_critical(model, X, ns, S, z0, zi, u)
         df = df_critical(model, X, ns, S, z0, zi, u)
         dX = solve_system(A=df, b=-F)

         do while(maxval(abs(dX)) > 1e-1)
            dX = dX/10
         end do

         if (maxval(abs(F)) < 1e-6) exit

         X = X + dX
         l = lambda1(model, X, 0.0_pr, z0, zi, u, u_new)
         u = u_new

         critical_point%iters = i
      end do

      z = X(1)*zi + (1-X(1))*z0
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
