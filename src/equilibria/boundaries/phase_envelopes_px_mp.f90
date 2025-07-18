module yaeos__equilibria_boundaries_phase_envelopes_mp_px
   !! Multiphase Px envelope calculation module.
   !!
   !! This module contains the functions to calculate the PT envelope of a
   !! mixture with multiple phases.
   use yaeos__constants, only: pr, R
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState
   use yaeos__models_ar, only: ArModel
   use yaeos__math, only: solve_system
   use yaeos__equilibria_boundaries_auxiliar, only: get_z, detect_critical

   implicit none

   private

   public :: PXEnvelMP
   public :: px_F_NP
   public :: px_envelope

   type :: PXEnvelMP
      !! Multiphase PX envelope.
      type(MPPoint), allocatable :: points(:) !! Array of converged points.
      real(pr), allocatable :: alpha(:) !! Molar relation between two mixtures.
      real(pr), allocatable :: z0(:) !! Original mixture mole fractions.
      real(pr), allocatable :: zi(:) !! Other mixture mole fractions
      real(pr), allocatable :: Pc(:) !! Critical pressure [bar]
      real(pr), allocatable :: ac(:) !! Critical alphas.
   contains
      procedure :: write => write_envelope_Px_MP
      procedure, nopass :: solve_point
      procedure, nopass :: get_values_from_X
   end type PXEnvelMP

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
      integer :: iters !! Number of iterations needed to converge the point.
      integer :: ns !! Number of the specified variable.
      character(len=14), allocatable :: kinds_x(:) !! Kinds of the main phases.
      character(len=14) :: kind_w !! Kind of the reference phase.
   end type MPPoint

contains
   type(PXEnvelMP) function px_envelope(&
      model, z0, zi, np, T, x_l0, w0, betas0, P0, alpha0, ns0, dS0, beta_w, points, kinds_x, kind_w &
      )
      !! # `px_envelope`
      !! Calculation of a multiphase Px envelope.
      !!
      !! # Description
      !! Calculates a phase envelope at costant temperature, using a numerical
      !! continuation method.
      use yaeos__auxiliar, only: optval
      class(ArModel), intent(in) :: model !! Model to use.
      real(pr), intent(in) :: z0(:) !! Original fluid composition.
      real(pr), intent(in) :: zi(:) !! Other fluid compostion.
      integer, intent(in) :: np
      !! Number of phases, without including the reference phaes
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: x_l0(np, size(z0))
      !! Initial guess for composition of phases.
      real(pr), intent(in) :: w0(size(z0))
      !! Initial guess for composition of reference phase.
      real(pr), intent(in) :: betas0(np)
      !! Mole fractions of each phase. Excluding the reference phase.
      real(pr), intent(in) :: P0 !! Initial guess for pressure [bar]
      real(pr), intent(in) :: alpha0
      !! Initial guess for relation between two fluids \(\alpha\)
      integer, intent(in) :: ns0
      !! First specified variable.
      !!
      !! The first `nc*np` values correspond to
      !! the different \(\ln K_i^l\), which are sorted like
      !! \([\ln K_1^1, \ln K_2^1, \dots \ln K_1^2, \dots, ln K_{nc}^{np}]\).
      !!
      !! From `nc*np+1` to `nc*np + np`, the different \(\beta^l\) values.
      !!
      !! `nc*np+np+1` and `cp*np+np+2` are \(P\) and \(\alpha\),
      !! respectively.
      real(pr), intent(in) :: dS0
      !! First step to extrapolate for next point calculation. After that
      !! It will use an adaptive algorithm.
      real(pr), intent(in) :: beta_w
      !! Fraction of the reference (incipient) phase.
      integer, optional, intent(in) :: points
      !! Maximum number of points to calculate.
      character(len=14), optional, intent(in) :: kinds_x(np)
      character(len=14), optional, intent(in) :: kind_w
      ! integer, optional, intent(in) :: kinds_x(np)
      ! integer, optional, intent(in) :: kind_w

      type(MPPoint), allocatable :: env_points(:)
      real(pr), allocatable :: alphas(:)
      type(MPPoint) :: point

      real(pr) :: F(size(z0) * np + np + 2)
      real(pr) :: dF(size(z0) * np + np + 2, size(z0) * np + np + 2)
      real(pr) :: dXdS(size(z0) * np + np + 2)
      real(pr) :: X(size(z0) * np + np + 2), dX(size(z0) * np + np + 2)
      real(pr) :: z(size(z0))

      integer :: nc

      integer :: its
      integer :: max_iterations = 100
      integer :: number_of_points

      real(pr) :: x_l(np, size(z0)), w(size(z0)), betas(np), P, alpha

      character(len=14) :: x_kinds(np)
      character(len=14) :: w_kind

      integer :: i !! Point calculation index
      integer :: lb !! Lower bound, index of the first component of a phase
      integer :: ub !! Upper bound, index of the last component of a phase
      integer :: inner !! Number of times a failed point is retried to converge

      integer :: ns !! Number of the specified variable
      real(pr) :: dS !! Step size of the specification for the next point
      real(pr) :: S !! Specified value

      real(pr) :: X0(size(X)) !! Initial guess for the point

      real(pr) :: Xold(size(X)) !! Old vector of variables
      real(pr) :: X_last_converged(size(X)) !! Last converged point
      real(pr) :: Xc(size(X)) !! Critical variables
      logical :: found_critical
      real(pr) :: ac, Pc

      integer :: ia, iP

      nc = size(z0)
      ia = nc*np+np+2
      iP = nc*np+np+1

      number_of_points = optval(points, 1000)

      if (present(kinds_x) .and. present(kind_w)) then
         x_kinds = kinds_x
         w_kind = kind_w
      else
         x_kinds = "stable"
         w_kind = "stable"
      end if

      do i=1,np
         lb = (i-1)*nc + 1
         ub = i*nc
         X(lb:ub) = log(x_l0(i, :)/w0)
      end do

      X(np*nc + 1:np*nc + np) = betas0
      X(np*nc + np + 1) = log(P0)
      X(np*nc + np + 2) = alpha0

      ns = ns0
      S = X(ns)
      dS = dS0

      allocate(env_points(0), alphas(0), px_envelope%ac(0), px_envelope%Pc(0))
      call solve_point(&
         model, z0, zi, np, T, beta_w, x_kinds, w_kind, &
         X, ns, S, dXdS, &
         F, dF, its, 1000 &
         )
      do i=1,number_of_points
         X0 = X
         call solve_point(&
            model=model, z0=z0, zi=zi, np=np, T=T, beta_w=beta_w, &
            x_kinds=x_kinds, w_kind=w_kind, &
            X=X, ns=ns, S=S, dXdS=dXdS, &
            F=F, dF=dF, iters=its, max_iterations=max_iterations &
            )

         ! The point might not converge, in this case we try again with an
         ! initial guess closer to the previous (converged) point.
         inner = 0
         do while(i > 1 .and. its >= max_iterations .and. inner < 10)
            inner = inner + 1
            X = X0 - (1 - real(inner, pr) / 10._pr) * dX
            S = X(ns)
            call solve_point(&
               model=model, z0=z0, zi=zi, np=np, T=T, beta_w=beta_w, &
               x_kinds=x_kinds, w_kind=w_kind, &
               X=X, ns=ns, S=S, dXdS=dXdS, &
               F=F, dF=dF, iters=its, max_iterations=5 &
               )
         end do

         ! If the point did not converge, stop the calculation
         if (any(isnan(F)) .or. its > max_iterations .or. abs(dS) < 1e-14) exit

         ! Save the information of the converged point
         call get_values_from_X(X, np, z0, zi, beta_w, x_l, w, betas, P, alpha)
         point = MPPoint(&
            np=np, nc=nc, betas=betas, P=P, T=T, x_l=x_l, w=w, beta_w=beta_w, &
            iters=its, ns=ns, kinds_x=x_kinds, kind_w=w_kind &
            )
         env_points = [env_points, point]
         alphas = [alphas, alpha]

         ! Update the specification for the next point.
         call update_specification(its, nc, np, X, dF, dXdS, ns, dS)

         ! Check if the system is close to a critical point, and try to jump
         ! over it.
         call detect_critical(nc, np, i, x_kinds, w_kind, .true., X_last_converged, X, dXdS, ns, dS, S, found_critical, Xc)
         if (found_critical) then
            ac = Xc(ia)
            Pc = exp(Xc(iP))
            px_envelope%Pc = [px_envelope%Pc, Pc]
            px_envelope%ac = [px_envelope%ac, ac]
         end if


         if (nc == 2) then
            alpha = X(ia) + dXdS(ia)*dS
            z = alpha * zi + (1- alpha) * z0

            do while((any(z < 0) .or. any(z > 1)) .and. abs(dS) > 0)
               dS = dS/2
               alpha = X(ia) + dXdS(ia)*dS
               z = alpha * zi + (1- alpha) * z0
            end do
         end if

         ! Next point estimation.
         dX = dXdS * dS
         
         X_last_converged = X
         X = X + dX
         S = X(ns)
      end do

      ! This moves the locally saved points to the output variable.
      call move_alloc(env_points, px_envelope%points)
      call move_alloc(alphas, px_envelope%alpha)
   end function px_envelope

   subroutine px_F_NP(model, z0, zi, np, T, beta_w, x_kinds, w_kind, X, ns, S, F, dF)
      !! # `px_F_NP`
      !! System of equations to solve a multiphase-point at constant
      !! temperature.
      !!
      !! # Description
      !! A multiphase equilibria point between `np+1` phases and `nc`
      !! components, where the `np+1` phase is a phase taken as reference for
      !! the calculation of equilibria rations
      !! \(K_i^l = \frac{\mathbf{x}_i^l}{\mathbf{w}_i}\),
      !! can be defined by the system of equations:
      !!
      !! \[
      !! \begin{bmatrix}
      !! \ln K_i^{l} + \ln \phi_i^{l}(\mathbf{x}^l, P, T) - \ln \phi_i^{l}(\mathbf{w}, P, T) \\
      !! \sum_i{\mathbf{x}^l_i - \mathbf{w}_i}
      !! \sum^l{\beta^l} + \beta^{np+1} - 1
      !! \end{bmatrix}
      !! \]
      ! ------------------------------------------------------------------------
      use iso_fortran_env, only: error_unit
      class(ArModel), intent(in) :: model !! Model to use.
      real(pr), intent(in) :: z0(:) !! First mixture composition.
      real(pr), intent(in) :: zi(:) !! Second mixture composition.
      integer, intent(in) :: np !! Number of main phases.
      real(pr), intent(in) :: T !! Temperature [K].
      real(pr), intent(in) :: beta_w !! Fraction of the reference (incipient) phase.
      real(pr), intent(in)  :: X(:) !! Vector of variables.
      integer, intent(in)  :: ns !! Number of specification.
      character(len=14), intent(in) :: x_kinds(np)
      character(len=14), intent(in) :: w_kind
      real(pr), intent(in)  :: S !! Specification value.
      real(pr), intent(out) :: F(size(X)) !! Vector of functions valuated.
      real(pr), intent(out) :: df(size(X), size(X)) !! Jacobian matrix.

      ! X variables
      real(pr) :: K(np, size(z0))
      real(pr) :: P
      real(pr) :: betas(np)

      real(pr) :: z(size(z0)), alpha, dzda(size(z0))

      ! Main phases variables
      real(pr) :: moles(size(z0))

      real(pr) :: Vl(np)
      real(pr), dimension(np, size(z0)) :: x_l, lnphi_l, dlnphi_dt_l, dlnphi_dp_l
      real(pr), dimension(np, size(z0), size(z0)) :: dlnphi_dn_l

      real(pr) :: lnphi(size(z0)), dlnphi_dt(size(z0)), dlnphi_dp(size(z0))
      real(pr), dimension(size(z0), size(z0)) :: dlnphi_dn

      ! Incipient phase variables
      real(pr) :: Vw
      real(pr), dimension(size(z0)) :: w, lnphi_w, dlnphi_dt_w, dlnphi_dp_w
      real(pr), dimension(size(z0), size(z0)) :: dlnphi_dn_w
      real(pr) :: dwda(size(z0))

      ! Derivatives of w wrt beta and K
      real(pr) :: dwdb(np, size(z0))
      real(pr) :: dwdlnK(np, size(z0))

      real(pr) :: denom(size(z0))
      real(pr) :: denomdlnK(np, size(z0), size(z0))


      real(pr) :: dx_l_dlnK(np, np, size(z0))

      integer :: i, j, l, phase, nc
      integer :: lb, ub
      integer :: idx_1, idx_2

      nc = size(z0)


      ! ========================================================================
      ! Extract variables from the vector X
      ! ------------------------------------------------------------------------
      P = exp(X(np*nc + np + 1))
      alpha = X(np*nc + np + 2)
      do l=1,np
         lb = (l-1)*nc + 1
         ub = l*nc
         K(l, :) = exp(X(lb:ub))
      end do
      betas = X(np*nc + 1:np*nc + np)

      call get_z(alpha, z0, zi, z, dzda)

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
         w, P, T, V=Vw, root_type=w_kind, lnphip=lnphi_w, &
         dlnphidp=dlnphi_dp_w, dlnphidt=dlnphi_dt_w, dlnphidn=dlnphi_dn_w &
         )

      do l=1,np
         x_l(l, :) = K(l, :)*w
         call model%lnphi_pt(&
            x_l(l, :), P, T, V=Vl(l), root_type=x_kinds(l), lnphip=lnphi, &
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
      dwda = dzda/denom

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
         ! disofug wrt P, alpha
         do i=1,nc
            lb = (l-1)*nc + i
            df(lb, nc*np+np+1) = P*(dlnphi_dp_l(l, i) - dlnphi_dp_w(i))
            df(lb, nc*np+np+2) = sum(&
               dwda * K(l, :) * dlnphi_dn_l(l, i, :) - dwda*dlnphi_dn_w(i, :) &
               )
         end do

         ! Derivatives of the sum of mole fractions

         ! wrt lnK
         do phase=1,np
            do j=1,nc
               lb = nc*np + phase
               ub = j + (l-1)*nc
               df(lb, ub) = df(lb, ub) + (dx_l_dlnK(phase, l, j) - dwdlnK(l, j))
            end do
         end do

         ! wrt beta, alpha
         do j=1,np
            lb = nc*np + j
            df(lb,np*nc+l) = sum(K(j, :) * dwdb(l, :) - dwdb(l, :))
         end do
         df(nc*np+l, nc*np+np+2) = sum(K(l, :) * dwda - dwda)

         ! Derivatives of sum(beta) + beta_w == 1
         df(nc * np + np + 1, np*nc + l) = 1
      end do



      df(nc * np + np + 2, ns) = 1
   end subroutine px_F_NP

   subroutine solve_point(model, z0, zi, np, T, beta_w, x_kinds, w_kind, X, ns, S, dXdS, F, dF, iters, max_iterations)
      !! # `solve_point`
      !! Solve the system of equations for a multiphase point.
      !!
      !! # Description
      !! Solves the point of a multiphase system using the Newton-Raphson
      !! method. The system of equations is defined in [[px_F_NP(procedure)]]
      use iso_fortran_env, only: error_unit
      use yaeos__math, only: solve_system
      class(ArModel), intent(in) :: model !! Model to use.
      real(pr), intent(in) :: z0(:) !! First mixture composition.
      real(pr), intent(in) :: zi(:) !! Second mixture composition.
      integer, intent(in) :: np !! Number of main phases
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: beta_w !! Fraction of the reference (incipient) phase
      character(len=14), intent(in) :: x_kinds(np)
      character(len=14), intent(in) :: w_kind
      real(pr), intent(in out)  :: X(:) !! Vector of variables
      integer, intent(in)  :: ns !! Number of specification
      real(pr), intent(in)  :: S !! Specification value
      real(pr), intent(in) :: dXdS(size(X))
      !! Sensitivity of the variables wrt the specification
      real(pr), intent(out) :: F(size(X)) !! Vector of functions valuated
      real(pr), intent(out) :: df(size(X), size(X)) !! Jacobian matrix
      integer, intent(in) :: max_iterations
      !! Maximum number of iterations to solve the point
      integer, intent(out) :: iters
      !! Number of iterations needed to converge the point


      integer :: ia
      integer :: iP
      integer :: nc

      real(pr) :: X0(size(X))
      real(pr) :: dX(size(X))

      nc = size(z0)
      iP = np*nc + np + 1
      ia = np*nc + np + 2

      X0 = X
      F = 1
      dX = 1

      do iters=1,max_iterations
         call px_F_NP(&
            model=model, z0=z0, zi=zi, np=np, T=T, beta_w=beta_w, &
            x_kinds=x_kinds, w_kind=w_kind, &
            X=X, ns=ns, S=S, F=F, dF=dF)

         if (any(isnan(F))) then
            X = X - 0.9 * dX
            cycle
         end if

         dX = solve_system(dF, -F)

         if (maxval(abs(F)) < 1e-9_pr) exit

         do while(abs( exp(X(iP))  - exp(X(iP) + dX(iP))) > 10)
            dX = dX*0.9
         end do

         X = X + dX
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

      integer :: ia !! Index of the \(\alpha\) variable
      integer :: iP !! Index of the \(\ln P\) variable
      integer :: lb !! Lower bound of each phase
      integer :: ub !! Upper bound of each phase

      dFdS = 0
      dFdS(size(X)) = -1

      dXdS = solve_system(dF, -dFdS)

      ns = maxloc(abs(dXdS), dim=1)

      iP = nc*np + np + 1
      ia = nc*np + np + 2

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
            exit
         end if
      end do

      dS = dXdS(ns) * dS
      dXdS = dXdS/dXdS(ns)

      ! We adapt the step size to the number of iterations, the desired number
      ! of iterations for each point is around 3.
      ! dS = dS * 3._pr/its

      do while(maxval(abs(dXdS(:nc*np)*dS)) > 0.1_pr)
         dS = dS/2
      end do

      do while(minval(abs(dXdS(:nc*np)*dS)) < 1e-5_pr)
         dS = dS*1.1
      end do

      do while(abs(dXdS(ia)*dS) < 1e-3 .and. abs(dXdS(iP)*dS) < 1e-2)
         dS = dS*1.1
      end do
   end subroutine update_specification

   subroutine get_values_from_X(X, np, z0, zi, beta_w, x_l, w, betas, P, alpha)
      !! # get_values_from_X
      !! Extract the values of the variables from the vector X.
      !!
      real(pr), intent(in) :: X(:) !! Vector of variables.
      integer, intent(in) :: np !! Number of main phases.
      real(pr), intent(in) :: z0(:) !! Initial mixture composition.
      real(pr), intent(in) :: zi(:) !! Second mixture composition.
      real(pr), intent(in) :: beta_w !! Reference phase beta.
      real(pr), intent(out) :: x_l(np, size(z0)) !! Mole fractions of the main phases.
      real(pr), intent(out) :: w(size(z0)) !! Mole fractions of the incipient phase.
      real(pr), intent(out) :: betas(np) !! Fractions of the main phases.
      real(pr), intent(out) :: P !! Pressure [bar].
      real(pr), intent(out) :: alpha !! \(alpha\).

      real(pr) :: z(size(z0))
      real(pr) :: K(np, size(z0))

      integer :: nc !! Number of components.
      integer :: i !! Loop index.
      integer :: l !! Phase index.
      integer :: lb !! Lower bound of each phase.
      integer :: ub !! Upper bound of each phase.

      nc = size(z0)

      ! ========================================================================
      ! Extract variables from the vector X
      ! ------------------------------------------------------------------------
      alpha = X(np*nc + np + 2)
      call get_z(alpha, z0, zi, z)
      P = exp(X(np*nc + np + 1))
      do l=1,np
         lb = (l-1)*nc + 1
         ub = l*nc
         K(l, :) = exp(X(lb:ub))
      end do
      betas = X(np*nc + 1:np*nc + np)

      w = z/(matmul(betas, K) + beta_w)

      do l=1,np
         x_l(l, :) = K(l, :) * w
      end do
   end subroutine get_values_from_X

   subroutine write_envelope_PX_MP(env, unit)
      class(PXEnvelMP), intent(in) :: env
      integer, intent(in) :: unit

      integer :: i, j
      integer :: np, nc
      real(pr) :: P, T, alpha
      real(pr), allocatable :: betas(:)
      real(pr), allocatable :: w(:)
      real(pr), allocatable :: x_l(:, :)

      np = size(env%points)
      nc = size(env%points(1)%w)

      do i=1,np
         alpha = env%alpha(i)
         P = env%points(i)%P
         T = env%points(i)%T
         betas = env%points(i)%betas
         w = env%points(i)%w
         x_l = env%points(i)%x_l
         write(unit, "(*(E15.5,2x))") alpha, P, T, betas, w, (x_l(j, :), j=1, size(x_l,dim=1))
      end do
   end subroutine write_envelope_PX_MP
end module yaeos__equilibria_boundaries_phase_envelopes_mp_px
