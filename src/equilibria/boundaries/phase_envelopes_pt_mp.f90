module yaeos__equilibria_boundaries_phase_envelopes_mp
   !! Multiphase PT envelope calculation module.
   !!
   !! This module contains the functions to calculate the PT envelope of a
   !! mixture with multiple phases.
   use yaeos__constants, only: pr, R
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState
   use yaeos__equilibria_boundaries_auxiliar, only: detect_critical
   use yaeos__models_ar, only: ArModel
   use yaeos__math, only: solve_system

   implicit none

   private

   public :: PTEnvelMP
   public :: pt_F_NP
   public :: pt_envelope

   type :: PTEnvelMP
      !! Multiphase PT envelope.
      type(MPPoint), allocatable :: points(:) !! Array of converged points.
      real(pr), allocatable :: Tc(:) !! Critical temperatures.
      real(pr), allocatable :: Pc(:) !! Critical pressures.
   contains
      procedure :: write => write_envelope_PT_MP
      procedure, nopass :: solve_point
      procedure, nopass :: get_values_from_X
   end type PTEnvelMP

   type :: MPPoint
      !! Multiphase equilibria point.
      integer :: np !! Number of phases
      integer :: nc !! Number of components
      real(pr) :: beta_w !! Fraction of the reference (incipient) phase.
      real(pr), allocatable :: betas(:) !! Fractions of the main phases.
      real(pr) :: P !! Pressure [bar]
      real(pr) :: T !! Temperature [K]
      real(pr), allocatable :: x_l(:, :) !! Mole fractions of the main phases.
      real(pr), allocatable :: w(:) !! Mole fractions of the incipient phase.
      character(len=14), allocatable :: kinds_x(:) !! Kinds of the main phases.
      character(len=14) :: kind_w !! Kind of the reference phase.
      integer :: iters !! Number of iterations needed to converge the point.
      integer :: ns !! Number of the specified variable.
   end type MPPoint

contains

   type(PTEnvelMP) function pt_envelope(&
      model, z, np, kinds_x, kind_w, x_l0, w0, betas0, P0, T0, ns0, dS0, beta_w, points, &
      max_pressure &
      )
      !! # `pt_envelope`
      !! Calculation of a multiphase PT envelope.
      !!
      !! # Description
      !! Calculates a PT envelope is calculated using the continuation method.
      !! The envelope is calculated by solving the system of equations for each
      !! point of the envelope. The system of equations is solved using the
      !! Newton-Raphson method.
      !!
      !! This function requires the system specification conditions, which are
      !! the fluid composition (\z\), the number of phases that are not
      !! incipient; defined as \(np\), proper intialization values, the
      !! variables that end with `0` are the initial guess; the mole fraction
      !! of the reference phase `beta_w` which when it is equal to 0 means that
      !! we are calculating a phase boundary.
      use yaeos__auxiliar, only: optval
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z(:)
      !! Mixture global composition.
      integer, intent(in) :: np
      !! Number of main phases.
      character(len=14), intent(in) :: kinds_x(np)
      !! Kind of the main phases.
      character(len=14), intent(in) :: kind_w
      !! Kind of the reference phase.
      real(pr), intent(in) :: x_l0(np, size(z))
      !! Initial guess for the mole fractions of each phase. arranged as
      !! an array of size `(np, nc)`, where nc is the number of components
      !! and `np` the number of main phases. Each row correspond to the
      !! composition of each main phaase.
      real(pr), intent(in) :: w0(size(z))
      !! Initial guess for the mole fractions of the
      !! reference/incipient phase.
      real(pr), intent(in) :: betas0(np)
      !! Initial guess for the fractions of the main phases. arranged as
      !! an array of size `(np)`, where `np` is the number of main phases.
      real(pr), intent(in) :: P0 !! Initial guess for the pressure [bar].
      real(pr), intent(in) :: T0 !! Initial guess for the temperature [K].
      integer, intent(in) :: ns0
      !! Number of the specified variable.
      !! The variable to be specified. This is the variable that will be
      !! used to calculate the first point of the envelope. The variable
      !! can be any of the variables in the vector X, but it is recommended
      !! to use the temperature or pressure. The variables are aranged as
      !! follows:
      !!
      !! - `X(1:nc*np) = ln(K_i^l)`: \(\frac{x_i^l}{w_i}\)
      !! - `X(nc*np+1:nc*np+np) = \beta_i^l`: Fraction of each main phase.
      !! - `X(nc*np+np+1) = ln(P)`: Pressure [bar].
      !! - `X(nc*np+np+2) = ln(T)`: Temperature [K].
      real(pr), intent(in) :: dS0
      !! Step size of the specification for the next point.
      !! This is the step size that will be used to calculate the next point.
      !! Inside the algorithm this value is modified to adapt the step size
      !! to facilitate the convergence of each point.
      real(pr), intent(in) :: beta_w
      !! Fraction of the reference (incipient) phase.
      integer, optional, intent(in) :: points
      !! Number of points to calculate.
      real(pr), optional, intent(in) :: max_pressure
      !! Maximum pressure [bar] to calculate.
      !! If the pressure of the point is greater than this value, the
      !! calculation is stopped.
      !! This is useful to avoid calculating envelopes that go to infinite
      !! values of pressure.

      type(MPPoint), allocatable :: env_points(:) !! Array of converged points.
      type(MPPoint) :: point !! Converged point.
      real(pr) :: max_P !! Maximum pressure [bar] to calculate.

      real(pr) :: F(size(z) * np + np + 2) !! Vector of functions valuated.
      real(pr) :: dF(size(z) * np + np + 2, size(z) * np + np + 2)
      !! Jacobian matrix.
      real(pr) :: dXdS(size(z) * np + np + 2)
      !! Sensitivity of the variables wrt the specification.
      real(pr) :: X(size(z) * np + np + 2)
      !! Vector of variables.
      real(pr) :: dX(size(z) * np + np + 2)
      !! Step for next point estimation.

      integer :: nc !! Number of components.

      integer :: its
      !! Number of iterations to solve the current point.
      integer :: max_iterations = 10
      !! Maximum number of iterations to solve the point.
      integer :: number_of_points
      !! Number of points to calculate.


      real(pr) :: x_l(np, size(z)) !! Mole fractions of the main phases.
      real(pr) :: w(size(z)) !! Mole fractions of the incipient phase.
      real(pr) :: betas(np) !! Fractions of the main phases.
      real(pr) :: P !! Pressure [bar].
      real(pr) :: T !! Temperature [K].

      integer :: i !! Point calculation index
      integer :: iT !! Index of the temperature variable.
      integer :: iP !! Index of the pressure variable.
      integer :: lb !! Lower bound, index of the first component of a phase
      integer :: ub !! Upper bound, index of the last component of a phase
      integer :: inner !! Number of times a failed point is retried to converge

      integer :: ns !! Number of the specified variable
      real(pr) :: dS !! Step size of the specification for the next point
      real(pr) :: S !! Specified value

      real(pr) :: X0(size(X)) !! Initial guess for the point
      real(pr) :: X_last_converged(size(X)) !! Last converged point
      real(pr) :: Xc(size(X)) !! Vector of variables at the critical point
      logical :: found_critical !! If true, a critical point was found

      character(len=14) :: x_kinds(np), w_kind

      real(pr) :: Tc !! Critical temperature [K]
      real(pr) :: Pc !! Critical pressure [bar]

      nc = size(z)
      iP = np*nc + np + 1
      iT = np*nc + np + 2

      number_of_points = optval(points, 1000)
      max_P = optval(max_pressure, 2000._pr)

      do i=1,np
         lb = (i-1)*nc + 1
         ub = i*nc
         X(lb:ub) = log(x_l0(i, :)/w0)
      end do

      X(np*nc + 1:np*nc + np) = betas0
      X(np*nc + np + 1) = log(P0)
      X(np*nc + np + 2) = log(T0)

      ns = ns0
      S = X(ns)
      dS = dS0

      x_kinds = kinds_x
      w_kind = kind_w

      allocate(env_points(0), pt_envelope%Tc(0), pt_envelope%Pc(0))

      F = 1
      its = 0
      X0 = X
      call solve_point(&
         model, z, np, beta_w, x_kinds, w_kind, X, ns, S, dXdS, &
         F, dF, its, 1000 &
         )
      do i=1,number_of_points
         X0 = X
         call solve_point(&
            model, z, np, beta_w, x_kinds, w_kind, X, ns, S, dXdS, &
            F, dF, its, max_iterations &
            )

         ! The point might not converge, in this case we try again with an
         ! initial guess closer to the previous (converged) point.
         inner = 0
         do while(i > 1 .and. its >= max_iterations .and. inner < 10)
            inner = inner + 1
            X = X0 - (1 - real(inner, pr) / 10._pr) * dX
            S = X(ns)
            call solve_point(&
               model, z, np, beta_w, x_kinds, w_kind, X, ns, S, dXdS, &
               F, dF, its, max_iterations&
               )
         end do

         ! Convert the values of the vector of variables into human-friendly
         ! variables.
         call get_values_from_X(X, np, z, beta_w, x_l, w, betas, P, T)

         ! Attach the new point to the envelope.
         point = MPPoint(&
            np=np, nc=nc, betas=betas, P=P, T=T, x_l=x_l, w=w, beta_w=beta_w, &
            kinds_x=x_kinds, kind_w=w_kind, iters=its, ns=ns &
            )

         ! Check if the system is close to a critical point, and try to jump
         ! over it.
         call detect_critical(&
            nc, np, i, x_kinds, w_kind, .false., &
            X_last_converged, X, dXdS, ns, dS, S, &
            found_critical, Xc)

         if (found_critical) then
            ! Save critical point
            Tc = exp(Xc(iT))
            Pc = exp(Xc(iP))
            pt_envelope%Tc = [pt_envelope%Tc, Tc]
            pt_envelope%Pc = [pt_envelope%Pc, Pc]
         end if

         ! Update the specification for the next point.
         call update_specification(its, nc, np, X, dF, dXdS, ns, dS)

         ! If the point did not converge, stop the calculation
         if (&
            any(isnan(F)) .or. its > max_iterations &
            .or. exp(X(nc*np+np+1)) < 1e-5 &
            .or. P > max_P  &
            .or. any(betas < -1e-14) .or. any(betas > 1 + 1e-14) &
            .or. abs(dS) <= 1e-14 &
            ) exit

         env_points = [env_points, point]

         ! Next point estimation.
         dX = dXdS * dS

         do while(abs(exp(X(iT))  - exp(X(iT) + dX(iT))) > 7)
            dX = dX/2
         end do

         ! do while(abs(exp(X(iP))  - exp(X(iP) + dX(iP))) > 5)
         !    dX = dX/2
         ! end do

         do while(abs(exp(X(iP))  - exp(X(iP) + dX(iP))) < 3 &
            .and. abs(exp(X(iT))  - exp(X(iT) + dX(iT))) < 3)
            dX = dX*1.1
         end do

         dejavu: block
            !! Getting too close to a local minima in the envelope, which can
            !! be hard to select corrects steps.
            real(pr) :: dPdT_1, dPdT_2, d2PdT2, dT2
            if (i > 4) then
               dPdT_1 = (P - env_points(i-1)%P) / (T - env_points(i-1)%T)
               dPdT_2 = (P - env_points(i-2)%P) / (T - env_points(i-2)%T)
               dT2 = (T - env_points(i-1)%T) * (env_points(i-2)%T - env_points(i-3)%T)

               d2PdT2 = (P - 2*env_points(i-1)%P + env_points(i-2)%P) / dT2

               if (abs(d2PdT2) > 0.05 .and. abs(dPdT_1) < 1.5) then
                  ns = iT
                  dS = dXdS(ns) * dS
                  dXdS = dXdS/dXdS(ns)
                  do while(abs(exp(X(iT) + dXdS(iT) * dS) - T) > 1)
                     dS = 0.9 * dS
                  end do
               end if
            end if
         end block dejavu

         X_last_converged = X
         X = X + dX
         S = X(ns)
      end do

      ! This moves the locally saved points to the output variable.
      call move_alloc(env_points, pt_envelope%points)
   end function pt_envelope

   subroutine pt_F_NP(model, z, np, beta_w, kinds_x, kind_w, X, ns, S, F, dF)
      !! Function to solve at each point of a multi-phase envelope.
      use iso_fortran_env, only: error_unit
      class(ArModel), intent(in) :: model !! Model to use.
      real(pr), intent(in) :: z(:) !! Mixture global composition.
      integer, intent(in) :: np !! Number of main phases.
      real(pr), intent(in) :: beta_w !! Fraction of the reference (incipient) phase.
      character(len=14), intent(in) :: kinds_x(np) !! Kind of the main phases.
      character(len=14), intent(in) :: kind_w !! Kind of the reference phase.
      real(pr), intent(in)  :: X(:) !! Vector of variables.
      integer, intent(in)  :: ns !! Number of specification.
      real(pr), intent(in)  :: S !! Specification value.
      real(pr), intent(out) :: F(size(X)) !! Vector of functions valuated.
      real(pr), intent(out) :: df(size(X), size(X)) !! Jacobian matrix.

      ! X variables
      real(pr) :: K(np, size(z))
      real(pr) :: P
      real(pr) :: T
      real(pr) :: betas(np)

      ! Main phases variables
      real(pr) :: moles(size(z))

      real(pr) :: Vl(np)
      real(pr), dimension(np, size(z)) :: x_l, lnphi_l, dlnphi_dt_l, dlnphi_dp_l
      real(pr), dimension(np, size(z), size(z)) :: dlnphi_dn_l

      real(pr) :: lnphi(size(z)), dlnphi_dt(size(z)), dlnphi_dp(size(z))
      real(pr), dimension(size(z), size(z)) :: dlnphi_dn

      ! Incipient phase variables
      real(pr) :: Vw
      real(pr), dimension(size(z)) :: w, lnphi_w, dlnphi_dt_w, dlnphi_dp_w
      real(pr), dimension(size(z), size(z)) :: dlnphi_dn_w

      ! Derivatives of w wrt beta and K
      real(pr) :: dwdb(np, size(z))
      real(pr) :: dwdlnK(np, size(z))

      real(pr) :: denom(size(z))
      real(pr) :: denomdlnK(np, size(z), size(z))


      real(pr) :: dx_l_dlnK(np, np, size(z))

      integer :: i, j, l, phase, nc
      integer :: lb, ub
      integer :: idx_1, idx_2

      nc = size(z)

      ! ========================================================================
      ! Extract variables from the vector X
      ! ------------------------------------------------------------------------
      do l=1,np
         lb = (l-1)*nc + 1
         ub = l*nc
         K(l, :) = exp(X(lb:ub))
      end do
      betas = X(np*nc + 1:np*nc + np)
      P = exp(X(np*nc + np + 1))
      T = exp(X(np*nc + np + 2))

      denom = 0
      denom = matmul(betas, K) + beta_w
      denomdlnK = 0
      do i=1,nc
         denomdlnK(:, i, i) = betas(:)*K(:, i)
      end do

      w = z/denom

      ! ========================================================================
      ! Calculation of fugacities coeficients and their derivatives
      ! ------------------------------------------------------------------------
      call model%lnphi_pt(&
         w, P, T, V=Vw, root_type=kind_w, lnphi=lnphi_w, &
         dlnphidp=dlnphi_dp_w, dlnphidt=dlnphi_dt_w, dlnphidn=dlnphi_dn_w &
         )

      do l=1,np
         x_l(l, :) = K(l, :)*w
         call model%lnphi_pt(&
            x_l(l, :), P, T, V=Vl(l), root_type=kinds_x(l), lnphi=lnphi, &
            dlnphidp=dlnphi_dp, dlnphidt=dlnphi_dt, dlnphidn=dlnphi_dn &
            )
         lnphi_l(l, :) = lnphi
         dlnphi_dn_l(l, :, :) = dlnphi_dn
         dlnphi_dt_l(l, :) = dlnphi_dt
         dlnphi_dp_l(l, :) = dlnphi_dp
      end do

      ! ========================================================================
      ! Calculation of the system of equations
      ! ------------------------------------------------------------------------
      do l=1,np
         ! Select the limits of the function
         lb = (l-1)*nc + 1
         ub = l*nc

         F(lb:ub) = X(lb:ub) + lnphi_l(l, :) - lnphi_w
         F(nc * np + l) = sum(x_l(l, :) - w)
      end do
      F(nc * np + np + 1) = sum(betas) + beta_w - 1
      F(nc * np + np + 2) = X(ns) - S


      ! ========================================================================
      ! Derivatives and Jacobian Matrix of the whole system
      ! ------------------------------------------------------------------------
      df = 0
      dwdlnK = 0

      do l=1,np
         ! Save the derivatives of w wrt beta and K of the incipient phase
         dwdb(l, :) = -z * K(l, :)/denom**2
         dwdlnK(l, :) = -K(l, :) * betas(l)*z/denom**2
      end do

      do l=1,np
         do phase=1,np
            dx_l_dlnK(phase, l, :) = dwdlnK(l, :) * K(phase, :)
            if (phase == l) then
               dx_l_dlnK(phase, l, :) = dx_l_dlnK(phase, l, :) + w * K(l, :)
            end if
         end do
      end do

      do l=1,np
         ! Derivatives of the isofugacity equations

         ! wrt lnK
         do phase=1,np
            do i=1, nc
               do j=1,nc

                  idx_1 = i + (phase-1)*nc
                  idx_2 = j + (l-1)*nc

                  df(idx_1, idx_2) = &
                     dlnphi_dn_l(phase, i, j) * dx_l_dlnK(phase, l, j) &
                     - dlnphi_dn_w(i, j) * dwdlnK(l, j)

                  if (i == j .and. phase == l) then
                     df(idx_1, idx_2) = df(idx_1, idx_2) + 1
                  end if

               end do
            end do
         end do

         ! wrt betas
         do j=1,np
            lb = (j-1)*nc + 1
            ub = j*nc
            do i=1,nc
               df(lb+i-1, np*nc + l) = &
                  sum(K(j, :) * dlnphi_dn_l(j, i, :)*dwdb(l, :) &
                  - dlnphi_dn_w(i, :)*dwdb(l, :))
            end do
         end do

         ! wrt T,p
         do i=1,nc
            lb = (l-1)*nc + i
            df(lb, nc*np+np+1) = P*(dlnphi_dp_l(l, i) - dlnphi_dp_w(i))
            df(lb, nc*np+np+2) = T*(dlnphi_dt_l(l, i) - dlnphi_dt_w(i))
         end do

         ! ====================================================================
         ! Derivatives of the sum of mole fractions
         ! --------------------------------------------------------------------

         ! wrt lnK
         do phase=1,np
            do j=1,nc
               lb = nc*np + phase
               ub = j + (l-1)*nc
               df(lb, ub) = df(lb, ub) + (dx_l_dlnK(phase, l, j) - dwdlnK(l, j))
            end do
         end do

         ! wrt beta
         do j=1,np
            lb = nc*np + j
            df(lb,np*nc+l) = sum(K(j, :) * dwdb(l, :) - dwdb(l, :))
         end do

         ! Derivatives of sum(beta)==1
         df(nc * np + np + 1, np*nc + l) = 1
      end do

      df(nc * np + np + 2, ns) = 1
   end subroutine pt_F_NP

   subroutine solve_point(model, z, np, beta_w, kinds_x, kind_w, X, ns, S, dXdS, F, dF, iters, max_iterations)
      use iso_fortran_env, only: error_unit
      use yaeos__math, only: solve_system
      class(ArModel), intent(in) :: model !! Model to use.
      real(pr), intent(in) :: z(:) !! Mixture global composition.
      integer, intent(in) :: np !! Number of main phases
      real(pr), intent(in) :: beta_w !! Fraction of the reference (incipient) phase
      character(len=14), intent(in) :: kinds_x(np) !! Kind of the main phases
      character(len=14), intent(in) :: kind_w !! Kind of the reference phase
      real(pr), intent(in out)  :: X(:) !! Vector of variables
      integer, intent(in)  :: ns !! Number of specification
      real(pr), intent(in)  :: S !! Specification value
      real(pr), intent(in) :: dXdS(size(X))
      real(pr), intent(out) :: F(size(X)) !! Vector of functions valuated
      real(pr), intent(out) :: df(size(X), size(X)) !! Jacobian matrix
      integer, intent(in) :: max_iterations
      !! Maximum number of iterations to solve the point
      integer, intent(out) :: iters
      !! Number of iterations to solve the current point


      integer :: i, l
      integer :: iT
      integer :: iP
      integer :: iBetas(np)
      integer :: nc

      real(pr) :: X0(size(X))
      real(pr) :: dX(size(X))

      logical :: can_solve

      nc = size(z)
      iP = np*nc + np + 1
      iT = np*nc + np + 2

      X0 = X

      can_solve = .true.

      iBetas = [(i, i=np*nc+1, np*nc+np)]

      do iters=1,max_iterations
         call pt_F_NP(model, z, np, beta_w, kinds_x, kind_w, X, ns, S, F, dF)

         if (any(isnan(F)) .and. can_solve) then
            X = X - 0.9 * dX
            can_solve = .false.
            cycle
         end if

         dX = solve_system(dF, -F)

         do l=1,np
            if (maxval(abs(X(l:nc*l))) < 1e-1) then
               do while(maxval(abs(dX(l:nc*l))) > 1e-1)
                  dX = dX/2
               end do
            end if
         end do

         do while(abs(dX(iT)) > 0.5)
            dX = dX/2
         end do

         do while(abs(dX(iP)) > 0.25)
            dX = dX/2
         end do

         do i=1,np
            if (X(iBetas(i)) /= 0) then
               do while(abs(dX(iBetas(i))/X(iBetas(i))) > 0.5)
                  dX = dX/2
               end do
            end if
         end do

         if (maxval(abs(F)) < 1e-9_pr) exit
         if (iters < 3) then
            X = X + 0.1 * dX
         else
            X = X + dX
         end if
      end do
   end subroutine solve_point

   subroutine update_specification(its, nc, np, X, dF, dXdS, ns, dS)
      !! # update_specification
      !! Change the specified variable for the next step.
      !!
      !! # Description
      !! Using the information of a converged point and the Jacobian matrix of
      !! the function. It is possible to determine the sensitivity of the
      !! variables with respect to the specification. This information is used
      !! to update the specification for the next point. Choosing the variable
      !! with the highest sensitivity.
      !! This can be done by solving the system of equations:
      !!
      !! \[
      !! J \frac{dX}{dS} + \frac{dF}{dS} = 0
      !! \]
      !!
      !! for the \( \frac{dX}{dS} \) vector. The variable with the highest value
      !! of \( \frac{dX}{dS} \) is chosen as the new specification.
      !!
      !! # References
      !!
      integer, intent(in) :: its
      !! Iterations to solve the current point.
      integer, intent(in) :: nc
      !! Number of components in the mixture.
      integer, intent(in) :: np
      !! Number of main phases.
      real(pr), intent(in out) :: X(:)
      !! Vector of variables.
      real(pr), intent(in out) :: dF(:, :)
      !! Jacobian matrix.
      real(pr), intent(in out) :: dXdS(:)
      !! Sensitivity of the variables wrt the specification.
      integer, intent(in out) :: ns
      !! Number of the specified variable.
      real(pr), intent(in out) :: dS
      !! Step size of the specification for the next point.

      real(pr) :: dFdS(size(X))
      !! Sensitivity of the functions wrt the specification.

      integer :: i
      integer :: lb !! Lower bound of each phase
      integer :: ub !! Upper bound of each phase

      integer :: iT
      integer :: iP
      integer :: iBetas(np)

      real(pr) :: dT, dP

      iBetas = [(i, i=np*nc+1, np*nc+np)]
      iP = size(X) - 1
      iT = size(X)

      dFdS = 0
      dFdS(size(X)) = -1

      dXdS = solve_system(dF, -dFdS)

      ns = maxloc(abs(dXdS), dim=1)

      ! ========================================================================
      ! For each phase, check if the mole fractions are too low.
      ! this can be related to criticality and it is useful to force the
      ! specification of compositions.
      ! ------------------------------------------------------------------------
      do i=1,np
         lb = (i-1)*nc + 1
         ub = i*nc

         if (maxval(abs(X(lb:ub))) < 0.1) then
            ns = lb + maxloc(abs(X(lb:ub)), dim=1) - 1
            dS = dXdS(ns) * dS * 0.1
            dXdS = dXdS/dXdS(ns)
            exit
         end if
      end do

      dS = dXdS(ns)*dS
      dXdS = dXdS/dXdS(ns)

      dS = sign(min(dS, 0.01_pr, sqrt(abs((X(ns)+1e-15)/5))), dS)
   end subroutine update_specification

   subroutine get_values_from_X(X, np, z, beta_w, x_l, w, betas, P, T)
      !! # get_values_from_X
      !! Extract the values of the variables from the vector X.
      !!
      real(pr), intent(in) :: X(:) !! Vector of variables.
      integer, intent(in) :: np !! Number of main phases.
      real(pr), intent(in) :: z(:) !! Mixture composition.
      real(pr), intent(in) :: beta_w !! Fraction of the reference phase.
      real(pr), intent(out) :: x_l(np, size(z)) !! Mole fractions of the main phases.
      real(pr), intent(out) :: w(size(z)) !! Mole fractions of the incipient phase.
      real(pr), intent(out) :: betas(np) !! Fractions of the main phases.
      real(pr), intent(out) :: P !! Pressure [bar].
      real(pr), intent(out) :: T !! Temperature [K].

      integer :: nc !! Number of components.
      integer :: i !! Loop index.
      integer :: lb !! Lower bound of each phase.
      integer :: ub !! Upper bound of each phase.

      nc = size(z)

      betas = X(np*nc + 1:np*nc + np)
      P = exp(X(np*nc + np + 1))
      T = exp(X(np*nc + np + 2))

      ! Extract the K values from the vector of variables
      do i=1,np
         lb = (i-1)*nc + 1
         ub = i*nc
         x_l(i, :) = exp(X(lb:ub))
      end do

      ! Calculate the mole fractions of the incipient phase
      w = z/(matmul(betas, x_l) + beta_w)

      ! Calculate the mole fractions of the main phases with the previously
      ! calculated K values
      do i=1,np
         x_l(i, :) = x_l(i, :)*w
      end do
   end subroutine get_values_from_X

   subroutine write_envelope_PT_MP(env, unit)
      class(PTEnvelMP), intent(in) :: env
      integer, intent(in) :: unit

      integer :: i, j
      integer :: np, nc
      real(pr) :: P, T
      real(pr), allocatable :: betas(:)
      real(pr), allocatable :: w(:)
      real(pr), allocatable :: x_l(:, :)

      np = size(env%points)
      nc = size(env%points(1)%w)

      do i=1,np
         P = env%points(i)%P
         T = env%points(i)%T
         betas = env%points(i)%betas
         w = env%points(i)%w
         x_l = env%points(i)%x_l
         write(unit, "(*(E15.5,2x))") P, T, betas, w, (x_l(j, :), j=1, size(x_l,dim=1))
      end do
   end subroutine write_envelope_PT_MP
end module yaeos__equilibria_boundaries_phase_envelopes_mp
