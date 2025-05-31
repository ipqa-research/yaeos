module yaeos__equilibria_critical
   use yaeos__constants, only: pr
   use yaeos__models, only: ArModel
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState

   implicit none

   private

   public :: critical_line
   public :: critical_point
   public :: CriticalLine
   public :: spec_CP
   public :: lambda1
   public :: F_critical
   public :: df_critical

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
      type(EquilibriumState) :: CEP !! Critical End Point
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
      model, a0, z0, zi, ns0, S0, dS0, &
      v0, t0, p0, &
      max_points, maxP, first_point, &
      stability_analysis &
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
      real(pr), optional, intent(in) :: v0 !! Initial volume [L/mol]
      real(pr), optional, intent(in) :: t0 !! Initial temperature [K]
      real(pr), optional, intent(in) :: p0 !! Initial pressure [bar]
      integer, optional, intent(in) :: max_points !! Maximum number of points
      real(pr), optional, intent(in) :: maxP !! Maximum pressure
      type(EquilibriumState), optional, intent(in) :: first_point
      logical, optional :: stability_analysis


      real(pr) :: u(size(z0)) !! eigen-vector
      real(pr) :: u_new(size(z0)) !! eigen-vector

      real(pr), allocatable :: XS_i(:), XS(:, :)!! Full set of solved points

      integer :: ns
      real(pr) :: S

      real(pr) :: X0(4), T, P, V, a, z(size(z0))

      integer :: i, npoints

      real(pr) :: max_P

      real(pr) :: y_cep(size(z0))
      real(pr) :: V_cep
      logical :: stab_anal
      logical :: found_cep

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

      if (.not. present(stability_analysis)) then
         stab_anal = .false.
      else
         stab_anal = stability_analysis
      end if

      u = (z0 + zi)/sum(z0 + zi)
      z = a0*zi + (1-a0)*z0

      if (present(t0)) then
         T = t0
      else
         T = sum(model%components%Tc * z)
      end if

      if (present(p0)) then
         P = p0
      else
         P = sum(model%components%Pc * z)
      end if

      if (present(v0)) then
         V = v0
      else
         call model%volume(n=z, P=P, T=T, V=V, root_type="vapor")
      end if

      X0 = [set_a(a0), log([v, T, P])]

      if (present(first_point)) then
         X0 = [&
            first_point%x(2), &
            log([first_point%Vx, first_point%T, first_point%P])]
      end if

      if (ns0 == spec_CP%a) then
         X0(ns0) = set_a(S0)
      else
         X0(ns0) = S0
      end if
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

         if (exp(X(4)) > max_P) then
            max_P = exp(X(4)) + 100
         end if


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
               (maxval(abs(dX)) > 1e-5 &
               .or. maxval(abs(F)) > 1e-5) &
               .and. real_its < 500)

               its = its + 1
               real_its = real_its + 1

               F = F_critical(model, X, ns, Si, z0, zi, u)
               dF = df_critical(model, X, ns, Si, z0, zi, u)
               dX = solve_system(dF, -F)

               do while(abs(maxval(dX(:))) > 0.01)
                  dX = dX/2
               end do

               X = X + dX
               l1 = lambda1(model=model, X=X, s=0.0_pr, z0=z0, zi=zi, u=u, u_new=u_new)
               u = u_new
            end do

            ! ==============================================================
            ! Cases where the line must be stopped
            ! --------------------------------------------------------------
            if (real_its == 500) exit
            if (any(isnan(X))) exit
            if (exp(X(spec_CP%P)) > max_P) exit


            a = get_a(X(1))
            V = exp(X(2))
            T = exp(X(3))
            P = exp(X(4))

            ! ==============================================================
            ! Stability analysis
            ! --------------------------------------------------------------
            if (stab_anal) then
               call look_for_cep(model, z0, zi, P, V, T, a, u, found_cep, critical_line%CEP)
               if (found_cep) then
                  exit solve_points
               end if
            end if

            ! ==============================================================
            ! Save point
            ! --------------------------------------------------------------
            critical_line%a = [critical_line%a, a]
            critical_line%V = [critical_line%V, V]
            critical_line%T = [critical_line%T, T]
            critical_line%P = [critical_line%P, P]
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

            ! dS = sign(min(abs(dS * 3./its), dS0), dS)
            ! ==============================================================
            ! Avoid big steps in pressure
            ! --------------------------------------------------------------
            ! do while(abs(exp(X(4) + dXdS(4) * dS) - exp(X(4))) > 20)
            !    dS = dS * 0.9
            ! end do

            ! Next step
            z = a * zi  + (1-a)*z0
            X = X + dXdS*dS
         end do
      end block solve_points

      critical_line%z0 = z0
      critical_line%zi = zi
   end function critical_line

   subroutine look_for_cep(model, z0, zi, Pc, Vc, Tc, a, u, found, CEP)
      use yaeos__math, only: solve_system
      class(ArModel), intent(in) :: model !! Equation of state model
      real(pr), intent(in) :: z0(:) !! Molar fractions of the first fluid
      real(pr), intent(in) :: zi(:) !! Molar fractions of the second fluid
      real(pr), intent(in) :: Pc !! Pressure [bar]
      real(pr), intent(in) :: Vc !! Volume [L/mol]
      real(pr), intent(in) :: Tc !! Temperature [K]
      real(pr), intent(in) :: a !! Molar fraction of the second fluid
      logical, intent(out) :: found !! Found a Critical End Point
      real(pr), intent(in out) :: u(:) !! Eigen-vector
      type(EquilibriumState), intent(out) :: CEP !! Critical End Point

      real(pr) :: y_cep(size(z0)), V_cep
      real(pr) :: Xcep(size(z0)+4),Fcep(size(z0)+4), dFcep(size(z0)+4, size(z0)+4), dXcep(size(z0)+4)

      found = .false.
      y_cep = 0
      V_cep = 0

      call stability_check(model, z0, zi, Pc, Vc, Tc, a, found, y_cep, V_cep)

      if (found) then
         Fcep = 1
         Xcep = [log(y_cep), log(V_cep), log(Vc), log(Tc), set_a(a)]
         do while(maxval(abs(Fcep)) > 1e-5)
            Fcep = F_cep(model, 2, X=Xcep, z0=z0, zi=zi, u=u)
            dFcep = df_cep(model, 2, X=Xcep, z0=z0, zi=zi, u=u)
            dXcep = solve_system(dFcep, -Fcep)

            do while(abs(dXcep(2+4)) > 0.01)
               dXcep(2+4) = dXcep(2+4)/2
            end do
            Xcep = Xcep + dXcep
         end do

         CEP%y = exp(Xcep(:2))
         CEP%Vy = exp(Xcep(3))
         CEP%Vx =  exp(Xcep(4))
         CEP%T = exp(Xcep(5))

         call model%pressure(n=CEP%y, V=CEP%Vy, T=CEP%T, P=CEP%P)
         CEP%x = zi * get_a(Xcep(6)) + (1-get_a(Xcep(6))) * z0
         CEP%kind = "CEP"
         CEP%beta = 0
      end if

   end subroutine look_for_cep

   subroutine stability_check(model, z0, zi, Pc, Vc, Tc, a, unstable, y_other, V_other)
      !! # stability_check
      !!
      !! ## Description
      !! Check the stability of a point in the critical line. The stability is
      !! determined by `tpd` analysis.
      use yaeos__equilibria_stability, only: min_tpd
      class(ArModel), intent(in) :: model !! Equation of state model
      real(pr), intent(in) :: z0(:) !! Molar fractions of the first fluid
      real(pr), intent(in) :: zi(:) !! Molar fractions of the second fluid
      real(pr), intent(in) :: Pc !! Pressure [bar]
      real(pr), intent(in) :: Vc !! Volume [L/mol]
      real(pr), intent(in) :: Tc !! Temperature [K]
      real(pr), intent(in) :: a !! Molar fraction of the second fluid
      logical, intent(out) :: unstable !! Stability of the point)
      real(pr), intent(out) :: V_other !! Volume [L/mol]
      real(pr), intent(out) :: y_other(:) !! Molar fractions of the second fluid

      real(pr) :: z(2)
      real(pr) :: y(2), dy
      real(pr) :: fug_z(2), fug_y(2), P
      integer :: istab, istab0
      real(pr) :: tpd

      logical :: first, possible

      if (size(z0) /= 2) then
         error stop "Stability check only for binary mixtures"
      end if

      z = a*zi + (1-a)*z0
      call model%lnfug_vt(n=z, V=Vc, T=Tc, lnf=fug_z)

      unstable = .false.
      first = .true.

      istab0 = 1

      ! TODO #optimization: Make an adaptative step of this
      possible = .false.
      istab0 = 0

      y(1) = 0
      y(2) = 1

      dy = 0.01_pr
      ! do while (istab0 < 2)
      !    istab = istab0

      do while(y(1) < 1-dy)
         y(1) = y(1) + dy
         y(2) = 1 - y(1)

         call model%volume(n=y, P=Pc, T=Tc, V=V_other, root_type="stable")
         call model%lnfug_vt(n=y, V=V_other, T=Tc, lnf=fug_y, P=P)

         tpd = sum(y * (fug_y - fug_z))
         if (tpd < -1e-2 ) then
            ! TODO: This should be finding the mimima
            unstable = .true.
            y_other = y
            return
         end if
      end do

   end subroutine stability_check

   real(pr) function lambda1(model, X, s, z0, zi, u, u_new, P)
      !! # lambda1
      !!
      !! Calculation of the first restriction of a critical point.
      !!
      !! \[
      !!  \lambda_1(s=0, \mathbf{n}, V, T) = \frac{d^2tpd}{ds^2} = 0
      !! \]
      !!
      !! \(\lambda_1\) is the smallest eigen-value for the matrix:
      !!
      !! \[
      !! M_{ij} = \sqrt{z_i z_j} \frac{d \ln f_i}{dn_j}(\mathbf{n}, V, T)
      !! \]
      !!
      !! Where
      !! \[
      !! \mathbf{n} = \mathbf{z} + s \mathbf{u} \sqrt{\mathbf{z}}
      !! \]
      !! And \( \mathbf{u} \) should be the eigen-vector corresponding to the
      !! smallest eigen-value when \(s = 0\)
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

      real(pr) :: a

      nc = size(z0)

      a = get_a(X(1))

      z = a * zi + (1-a)*z0
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
      if (present(u_new)) then
         u_new = vectors(:, minloc(abs(lambda), dim=1))
         u_new = u_new/sqrt(sum(u_new**2))
      end if
      if (present(P)) P = Pin
   end function lambda1

   function F_critical(model, X, ns, S, z0, zi, u)
      !! # F_critical
      !!
      !! ## Description
      !! Function that should be equal to zero at a critical point is found.
      !! The second criticality condition is calculated as a numerical
      !! derivative with `eps=1e-4`.
      !!
      !! \[
      !! F = \begin{bmatrix}
      !!   \lambda_1(s) \\
      !!   \frac{\partial \lambda_1(s+\epsilon) - \lambda_1(s-\epsilon)}{2\epsilon} \\
      !!   \ln P - X_4 \\
      !!   X_{ns} - S
      !! \end{bmatrix} = 0
      !! \]
      !!
      !! The vector of varibles is
      !!
      !! \[
      !! X = [\alpha, \ln V, \ln T, \ln P]
      !! \]
      !!
      !! Including internally the extra equation:
      !! \[ \mathbf{z} = \alpha \mathbf{z_i} + (1-\alpha) \mathbf{z_0} \]
      class(ArModel), intent(in) :: model !! Equation of state model
      real(pr), intent(in) :: X(4) !! Vector of variables
      integer, intent(in) :: ns !! Position of the specification variable
      real(pr), intent(in) :: S !! Specification variable value
      real(pr), intent(in) :: z0(:) !! Molar fractions of the first fluid
      real(pr), intent(in) :: zi(:) !! Molar fractions of the second fluid
      real(pr), intent(in out) :: u(:) !! Eigen-vector

      real(pr) :: u_new(size(u))

      real(pr) :: F_critical(4)
      real(pr) :: z(size(u))

      real(pr) :: a, V, T, P

      real(pr), parameter :: eps=1e-4_pr

      integer :: i
      real(pr) :: eps_df, F1(4), F2(4), dx(4)

      a = get_a(X(1))
      V = exp(X(2))
      T = exp(X(3))
      z = a * zi + (1-a) * z0

      ! if(any(z < 0) ) return

      F_critical(1) = lambda1(model=model, X=X, s=0.0_pr, z0=z0, zi=zi, u=u, P=P, u_new=u_new)
      if (size(z) == 2) u = u_new
      F_critical(2) = (&
         lambda1(model=model, X=X, s=eps, zi=zi, z0=z0, u=u) &
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
      real(pr), intent(in out) :: u(:) !! Eigen-vector
      real(pr) :: df_critical(4, 4) !! Jacobian of the critical point function

      real(pr) :: eps, a

      real(pr) :: dx(4), F1(4), F2(4)

      integer :: i

      ! if (any(X(1)*zi + (1-X(1))*z0 > 0.99)) then
      !    eps = 1e-3_pr
      ! else
      !    eps = 1e-6_pr
      ! end if

      a = get_a(X(1))

      if (size(zi) == 2) then
         eps = 1e-10
      else
         if (any(a*zi + (1-a)*z0 > 0.99)) then
            eps = 1e-3_pr
         else
            eps = 1e-6_pr
         end if
      end if

      df_critical = 0
      do i=1,4
         dx = 0
         dx(i) = max(abs(eps * X(i)), eps)
         F2 = F_critical(model, X+dx, ns, S, z0, zi, u)
         F1 = F_critical(model, X-dx, ns, S, z0, zi, u)
         df_critical(:, i) = (F2 - F1)/(2*dx(i))
      end do
   end function df_critical

   function F_cep(model, nc, X, z0, zi, u)
      class(ArModel), intent(in) :: model !! Equation of state model
      real(pr), intent(in) :: z0(nc) !! Molar fractions of the first fluid
      real(pr), intent(in) :: X(nc + 4) !! Vector of variables
      real(pr), intent(in) :: zi(nc) !! Molar fractions of the second fluid
      real(pr), intent(in out) :: u(nc) !! Eigen-vector

      real(pr) :: F_cep(nc+4)
      real(pr) :: z(nc)

      real(pr) :: V, T, P
      real(pr) :: Xcp(nc+4)

      real(pr) :: Vc, Pc, lnf_z(nc)
      real(pr) :: y(nc)
      real(pr) :: Vy, Py
      real(pr) :: lnf_y(nc)

      real(pr), parameter :: eps=1e-5_pr

      real(pr) :: a, u_new(nc)

      integer, intent(in) :: nc

      ! nc = size(z0)

      y = exp(X(:nc))
      Vy = exp(X(nc+1))
      Vc = exp(X(nc+2))
      T = exp(X(nc+3))
      a = get_a(X(nc+4))
      z = a * zi + (1-a) * z0

      if(any(z < 0) ) return

      call model%lnfug_vt(n=y, V=Vy, T=T, P=Py, lnf=lnf_y)
      call model%lnfug_vt(n=z, V=Vc, T=T, P=Pc, lnf=lnf_z)
      Xcp(1) = X(nc+4)
      Xcp(2) = log(Vc)
      Xcp(3) = log(T)
      Xcp(4) = log(Pc)

      F_cep(1) = lambda1(model=model, X=Xcp, s=0._pr, z0=z0, zi=zi, u=u, P=Pc, u_new=u_new)
      u = u_new
      F_cep(2) = (&
         lambda1(model=model, X=Xcp, s=eps, zi=zi, z0=z0, u=u) &
         - lambda1(model=model, X=Xcp, s=-eps, zi=zi, z0=z0, u=u))/(2*eps)
      F_cep(3) = log(Pc) - log(Py)
      F_cep(4:nc+3) = lnf_y - lnf_z
      F_cep(nc+4) = sum(y) - 1
   end function F_cep

   function df_cep(model, nc, X, z0, zi, u)
      !! # df_critical
      !!
      !! ## Description
      !! Calculates the Jacobian of the critical point function `F_critical`.
      class(ArModel), intent(in) :: model !! Equation of state model
      integer, intent(in) :: nc
      real(pr), intent(in) :: z0(nc) !! Molar fractions of the first fluid
      real(pr), intent(in) :: X(nc+4) !! Vector of variables
      real(pr), intent(in) :: zi(nc) !! Molar fractions of the second fluid
      real(pr), intent(in out) :: u(nc) !! Eigen-vector
      real(pr) :: df_cep(nc+4, nc+4) !! Jacobian of the critical point function

      real(pr) :: eps

      real(pr) :: dx(nc+4), F1(nc+4), F2(nc+4)

      real(pr) :: a

      integer :: i

      a = get_a(X(1))

      if (any(a*zi + (1-a)*z0 > 0.99)) then
         eps = 1e-3_pr
      else
         eps = 1e-6_pr
      end if

      eps = 1e-10

      df_cep = 0
      do i=1,size(df_cep, 1)
         dx = 0
         dx(i) = eps
         F2 = F_cep(model, nc, X+dx, z0, zi, u)
         F1 = F_cep(model, nc, X-dx, z0, zi, u)
         df_cep(:, i) = (F2 - F1)/(2*eps)
      end do
   end function df_cep

   type(EquilibriumState) function critical_point(&
      model, z0, zi, spec, S, max_iters, V0, T0, a0, P0 &
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

      real(pr) :: X(4)
      integer :: ns
      real(pr) :: F(4), df(4, 4), dX(4), u(size(z0))

      real(pr) :: z(size(z0)), u_new(size(z0)), l, a
      real(pr) :: Sin
      integer :: i


      ! ========================================================================
      ! Handle the input
      ! ------------------------------------------------------------------------
      if (present(a0)) then
         X(1) = set_a(a0)
      else if (spec == spec_CP%a) then
         X(1) = set_a(S)
      else
         X(1) = log(0.5_pr)
      end if

      a = get_a(X(1))
      z = a*zi + (1-a)*z0

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
      if (ns == spec_CP%a) then
         Sin = set_a(S)
         X(ns) = Sin
      else
         Sin = S
         X(ns) = Sin
      end if

      ! ========================================================================
      ! Solve the system of equations
      ! ------------------------------------------------------------------------
      do i=1,max_iters
         l = lambda1(model=model, X=X, s=0.0_pr, z0=z0, zi=zi, u=u, u_new=u_new)
         u = u_new

         F = F_critical(model, X, ns, Sin, z0, zi, u)
         df = df_critical(model, X, ns, Sin, z0, zi, u)
         dX = solve_system(A=df, b=-F)

         do while(maxval(abs(dX)) > 2e-2)
            dX = dX*0.99
         end do

         do while((get_a(X(1) + dX(1)) > 1 .or. get_a(X(1) + dX(1)) < 0) .and. size(z0) == 2)
            dX = dX/2
         end do

         if (maxval(abs(F)) < 1e-8) exit

         X = X + dX
         critical_point%iters = i
      end do

      a = get_a(X(1))
      z = a*zi + (1-a)*z0
      critical_point%x = z
      critical_point%y = z
      critical_point%beta = 0
      critical_point%Vx = exp(X(2))
      critical_point%Vy = exp(X(2))
      critical_point%T  = exp(X(3))
      call model%pressure(n=z, V=critical_point%Vx, T=critical_point%T, P=critical_point%P)
      critical_point%kind = "critical"

   end function critical_point

   real(pr) function get_a(X)
      real(pr), intent(in) :: X
      get_a = X!(X)**2
   end function get_a

   real(pr) function set_a(a)
      real(pr), intent(in) :: a
      set_a = a!sqrt(a)
   end function set_a
end module yaeos__equilibria_critical
