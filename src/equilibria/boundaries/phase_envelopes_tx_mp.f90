module yaeos__equilibria_boundaries_phase_envelopes_mp_tx
   !! Multiphase Px envelope calculation module.
   !!
   !! This module contains the functions to calculate the PT envelope of a
   !! mixture with multiple phases.
   use yaeos__constants, only: pr, R
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState
   use yaeos__models_ar, only: ArModel
   use yaeos__math, only: solve_system
   use yaeos__equilibria_boundaries_auxiliar, only: get_z

   implicit none

   private

   public :: TXEnvelMP
   public :: tx_F_NP
   public :: tx_envelope

   type :: TXEnvelMP
      !! Multiphase PT envelope.
      type(MPPoint), allocatable :: points(:) !! Array of converged points.
      real(pr), allocatable :: alpha(:)
      real(pr), allocatable :: z0(:)
      real(pr), allocatable :: zi(:)
   contains
      procedure :: write => write_envelope_Tx_MP
      procedure, nopass :: solve_point
      procedure, nopass :: get_values_from_X
   end type TXEnvelMP

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
   end type MPPoint

contains
   type(TXEnvelMP) function tx_envelope(&
      model, z0, zi, np, P, kinds_x, kind_w, &
      x_l0, w0, betas0, T0, alpha0, ns0, dS0, beta_w, points &
      )
      use yaeos__auxiliar, only: optval
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z0(:)
      real(pr), intent(in) :: zi(:)
      integer, intent(in) :: np
      real(pr), intent(in) :: P
      character(len=14), intent(in) :: kinds_x(np)
      character(len=14), intent(in) :: kind_w
      real(pr), intent(in) :: x_l0(np, size(z0))
      real(pr), intent(in) :: w0(size(z0))
      real(pr), intent(in) :: betas0(np)
      real(pr), intent(in) :: T0
      real(pr), intent(in) :: alpha0
      integer, intent(in) :: ns0
      real(pr), intent(in) :: dS0
      real(pr), intent(in) :: beta_w
      !! Fraction of the reference (incipient) phase.
      integer, optional, intent(in) :: points

      type(MPPoint), allocatable :: env_points(:)
      real(pr), allocatable :: alphas(:)
      type(MPPoint) :: point

      real(pr) :: F(size(z0) * np + np + 2)
      real(pr) :: dF(size(z0) * np + np + 2, size(z0) * np + np + 2)
      real(pr) :: dXdS(size(z0) * np + np + 2)
      real(pr) :: X(size(z0) * np + np + 2), dX(size(z0) * np + np + 2)

      integer :: nc

      integer :: its
      integer :: max_iterations = 100
      integer :: number_of_points


      real(pr) :: x_l(np, size(z0)), w(size(z0)), betas(np), T, alpha

      integer :: i !! Point calculation index
      integer :: lb !! Lower bound, index of the first component of a phase
      integer :: ub !! Upper bound, index of the last component of a phase
      integer :: inner !! Number of times a failed point is retried to converge

      integer :: ns !! Number of the specified variable
      real(pr) :: dS !! Step size of the specification for the next point
      real(pr) :: S !! Specified value

      real(pr) :: X0(size(X)) !! Initial guess for the point

      integer :: ia
      real(pr) :: z(size(z0))

      character(len=14) :: x_kinds(np)
      character(len=14) :: w_kind

      nc = size(z0)
      ia = np*nc + np + 2

      number_of_points = optval(points, 1000)

      do i=1,np
         lb = (i-1)*nc + 1
         ub = i*nc
         X(lb:ub) = log(x_l0(i, :)/w0)
      end do

      X(np*nc + 1:np*nc + np) = betas0
      X(np*nc + np + 1) = log(T0)
      X(np*nc + np + 2) = alpha0

      ns = ns0
      S = X(ns)
      dS = dS0

      x_kinds = kinds_x
      w_kind = kind_w

      allocate(env_points(0), alphas(0))
      call solve_point(&
         model=model, z0=z0, zi=zi, np=np, P=P, beta_w=beta_w, kinds_x=x_kinds, kind_w=w_kind, &
         X=X, ns=ns, S=S, dXdS=dXdS, &
         F=F, dF=dF, iters=its, max_iterations=1000 &
         )
      do i=1,number_of_points
         X0 = X
         call solve_point(&
            model=model, z0=z0, zi=zi, np=np, P=P, beta_w=beta_w, kinds_x=x_kinds, kind_w=w_kind, &
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
               model=model, z0=z0, zi=zi, np=np, P=P, beta_w=beta_w, kinds_x=x_kinds, kind_w=w_kind, &
               X=X, ns=ns, S=S, dXdS=dXdS, &
               F=F, dF=dF, iters=its, max_iterations=5 &
               )
         end do

         ! If the point did not converge, stop the calculation
         if (any(isnan(F)) .or. its > max_iterations) exit

         ! Save the information of the converged point
         call get_values_from_X(X, np, z0, zi, beta_w, x_l, w, betas, T, alpha)
         point = MPPoint(&
            np=np, nc=nc, betas=betas, P=P, T=T, x_l=x_l, w=w, beta_w=beta_w, &
            iters=its, ns=ns &
            )
         env_points = [env_points, point]
         alphas = [alphas, alpha]

         ! Update the specification for the next point.
         call update_specification(its, nc, np, X, dF, dXdS, ns, dS)

         ! Check if the system is close to a critical point, and try to jump
         ! over it.
         call detect_critical(&
            nc=nc, np=np, kinds_x=x_kinds, kind_w=w_kind, &
            X=X, dXdS=dXdS, ns=ns, dS=dS, S=S)

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
         X = X + dX
         S = X(ns)
      end do

      ! This moves the locally saved points to the output variable.
      call move_alloc(env_points, tx_envelope%points)
      call move_alloc(alphas, tx_envelope%alpha)
   end function tx_envelope

   subroutine tx_F_NP(model, z0, zi, np, P, beta_w, kinds_x, kind_w, &
      X, ns, S, F, dF)
      !! Function to solve at each point of a multi-phase envelope.
      use iso_fortran_env, only: error_unit
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z0(:)
      real(pr), intent(in) :: zi(:)
      integer, intent(in) :: np !! Number of main phases.
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: beta_w !! Fraction of the reference (incipient) phase.
      character(len=14), intent(in) :: kinds_x(np) !! Kinds of the main phases.
      character(len=14), intent(in) :: kind_w !! Kind of the reference phase.
      real(pr), intent(in)  :: X(:) !! Vector of variables.
      integer, intent(in)  :: ns !! Number of specification.
      real(pr), intent(in)  :: S !! Specification value.
      real(pr), intent(out) :: F(size(X)) !! Vector of functions valuated.
      real(pr), intent(out) :: df(size(X), size(X)) !! Jacobian matrix.

      ! X variables
      real(pr) :: K(np, size(z0))
      real(pr) :: T
      real(pr) :: betas(np)

      real(pr) :: z(size(z0)), alpha, dzda(size(z0))

      ! Main phases variables
      real(pr) :: moles(size(z0))

      real(pr) :: Vl(np)
      real(pr), dimension(np, size(z0)) :: x_l, lnphi_l, dlnphi_dt_l
      real(pr), dimension(np, size(z0), size(z0)) :: dlnphi_dn_l

      real(pr) :: lnphi(size(z0)), dlnphi_dt(size(z0))
      real(pr), dimension(size(z0), size(z0)) :: dlnphi_dn

      ! Incipient phase variables
      real(pr) :: Vw
      real(pr), dimension(size(z0)) :: w, lnphi_w, dlnphi_dt_w
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
      T = exp(X(np*nc + np + 1))
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
         w, P, T, V=Vw, root_type=kind_w, lnphi=lnphi_w, &
         dlnphidt=dlnphi_dt_w, dlnphidn=dlnphi_dn_w &
         )

      do l=1,np
         x_l(l, :) = K(l, :)*w
         call model%lnphi_pt(&
            x_l(l, :), P, T, V=Vl(l), root_type=kinds_x(l), lnphi=lnphi, &
            dlnphidt=dlnphi_dt, dlnphidn=dlnphi_dn &
            )
         lnphi_l(l, :) = lnphi
         dlnphi_dn_l(l, :, :) = dlnphi_dn
         dlnphi_dt_l(l, :) = dlnphi_dt
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
            df(lb, nc*np+np+1) = T*(dlnphi_dt_l(l, i) - dlnphi_dt_w(i))
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
   end subroutine tx_F_NP

   subroutine solve_point(&
      model, z0, zi, np, P, beta_w, kinds_x, kind_w, &
      X, ns, S, dXdS, F, dF, iters, max_iterations)
      use iso_fortran_env, only: error_unit
      use yaeos__math, only: solve_system
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z0(:)
      real(pr), intent(in) :: zi(:)
      integer, intent(in) :: np !! Number of main phases
      real(pr), intent(in) :: P
      real(pr), intent(in) :: beta_w !! Fraction of the reference (incipient) phase
      character(len=14), intent(in) :: kinds_x(np) !! Kinds of the main phases
      character(len=14), intent(in) :: kind_w !! Kind of the reference phase
      real(pr), intent(in out)  :: X(:) !! Vector of variables
      integer, intent(in)  :: ns !! Number of specification
      real(pr), intent(in)  :: S !! Specification value
      real(pr), intent(in) :: dXdS(size(X))
      real(pr), intent(out) :: F(size(X)) !! Vector of functions valuated
      real(pr), intent(out) :: df(size(X), size(X)) !! Jacobian matrix
      integer, intent(in) :: max_iterations
      integer, intent(out) :: iters


      integer :: ia
      integer :: iT
      integer :: nc

      real(pr) :: X0(size(X))
      real(pr) :: dX(size(X))

      nc = size(z0)
      iT = np*nc + np + 1
      ia = np*nc + np + 2

      X0 = X

      do iters=1,max_iterations
         call tx_F_NP(&
            model=model, z0=z0, zi=zi, np=np, P=P, beta_w=beta_w, &
            kinds_x=kinds_x, kind_w=kind_w, X=X, ns=ns, S=S, F=F, dF=dF)

         if (any(isnan(F))) then
            X = X - 0.9 * dX
            cycle
         end if

         dX = solve_system(dF, -F)

         if (maxval(abs(F)) < 1e-9_pr) exit

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
      integer :: lb !! Lower bound of each phase
      integer :: ub !! Upper bound of each phase

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

         if (maxval(abs(X(lb:ub))) < 0.3) then
            ns = lb + maxloc(abs(X(lb:ub)), dim=1) - 1
            exit
         end if
      end do
      dS = dXdS(ns) * dS
      dXdS = dXdS/dXdS(ns)

      ! We adapt the step size to the number of iterations, the desired number
      ! of iterations for each point is around 3.
      dS = dS * 3._pr/its
   end subroutine update_specification

   subroutine detect_critical(nc, np, kinds_x, kind_w, X, dXdS, ns, dS, S)
      !! # detect_critical
      !! Detect if the system is close to a critical point.
      !!
      !! # Description
      !! When the system is close to a critical point, the \(\ln K_i^l\) values
      !! are close to zero, since the composition of the incipient phase and the
      !! \(l\) phase are similar (equal in the critical point). This can be used
      !! to detect if the system is close to a critical point and force a jump
      !! above it.
      !!
      !! # References
      !!
      integer, intent(in) :: nc
      !! Number of components in the mixture.
      integer, intent(in) :: np
      !! Number of main phases.
      character(len=14), intent(in out) :: kinds_x(np)
      !! Kinds of the main phases.
      character(len=14), intent(in out) :: kind_w
      !! Kind of the reference phase.
      real(pr), intent(in out) :: X(:)
      !! Vector of variables.
      real(pr), intent(in out) :: dXdS(:)
      !! Sensitivity of the variables wrt the specification.
      integer, intent(in out) :: ns
      !! Number of the specified variable.
      real(pr), intent(in out) :: dS
      !! Step size of the specification for the next point.
      real(pr), intent(in out) :: S
      !! Specification value.

      integer :: i, lb, ub
      character(len=14) :: incipient_kind

      real(pr) :: Xold(nc*np+np+2)

      do i=1,np
         lb = (i-1)*nc + 1
         ub = i*nc

         do while(maxval(abs(X(lb:ub))) < 0.01)
            X = X + dXdS * dS
         end do

         if (all(Xold(lb:ub) * (X(lb:ub) + dXdS(lb:ub)*dS) < 0)) then
            incipient_kind = kind_w
            kind_w = kinds_x(i)
            kinds_x(i) = incipient_kind
            ! Interpolate here

            if (nc == 2) then
               dS=0
               return
            end if
         end if

      end do
   end subroutine detect_critical

   subroutine get_values_from_X(X, np, z0, zi, beta_w, x_l, w, betas, T, alpha)
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
      real(pr), intent(out) :: T !! Pressure [bar].
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
      T = exp(X(np*nc + np + 1))
      alpha = X(np*nc + np + 2)
      do l=1,np
         lb = (l-1)*nc + 1
         ub = l*nc
         K(l, :) = exp(X(lb:ub))
      end do
      betas = X(np*nc + 1:np*nc + np)


      w = z/matmul(betas, K) + beta_w

      do l=1,np
         x_l(l, :) = K(l, :) * w
      end do
   end subroutine get_values_from_X

   subroutine write_envelope_TX_MP(env, unit)
      class(TXEnvelMP), intent(in) :: env
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
   end subroutine write_envelope_TX_MP
end module yaeos__equilibria_boundaries_phase_envelopes_mp_tx
