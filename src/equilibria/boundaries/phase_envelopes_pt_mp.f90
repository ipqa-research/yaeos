module yaeos__equilibria_boundaries_phase_envelopes_mp
   !! Multiphase PT envelope calculation module.
   !!
   !! This module contains the functions to calculate the PT envelope of a
   !! mixture with multiple phases.
   use yaeos__constants, only: pr, R
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState
   use yaeos__models_ar, only: ArModel
   use yaeos__math, only: solve_system

   implicit none

   private

   public :: PTEnvelMP
   public :: pt_F_NP
   public :: solve_point
   public :: pt_envelope
   public :: get_values_from_X

   real(pr), public :: beta_w = 0

   type :: PTEnvelMP
      !! Multiphase PT envelope.
      type(MPPoint), allocatable :: points(:) !! Array of converged points.
   contains
      procedure :: write => write_envelope_PT_MP
   end type PTEnvelMP

   type :: MPPoint
      !! Multiphase equilibria point.
      integer :: np !! Number of phases
      integer :: nc !! Number of components
      real(pr), allocatable :: betas(:) !! Fractions of the main phases.
      real(pr) :: P !! Pressure [bar]
      real(pr) :: T !! Temperature [K]
      real(pr), allocatable :: x_l(:, :) !! Mole fractions of the main phases.
      real(pr), allocatable :: w(:) !! Mole fractions of the incipient phase.
   end type MPPoint

contains

   type(PTEnvelMP) function pt_envelope(&
      model, z, np, x_l0, w0, betas0, P0, T0, ns0, dS0, points&
      )
      use yaeos__auxiliar, only: optval
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z(:)
      integer, intent(in) :: np
      real(pr), intent(in) :: x_l0(np, size(z))
      real(pr), intent(in) :: w0(size(z))
      real(pr), intent(in) :: betas0(np)
      real(pr), intent(in) :: P0
      real(pr), intent(in) :: T0
      integer, intent(in) :: ns0
      real(pr), intent(in) :: dS0
      integer, optional, intent(in) :: points

      type(MPPoint), allocatable :: env_points(:)
      type(MPPoint) :: point

      real(pr) :: F(size(z) * np + np + 2)
      real(pr) :: dF(size(z) * np + np + 2, size(z) * np + np + 2)
      real(pr) :: dXdS(size(z) * np + np + 2)
      real(pr) :: X(size(z) * np + np + 2), dX(size(z) * np + np + 2)

      integer :: nc

      integer :: its
      integer :: max_iterations = 1000
      integer :: number_of_points


      real(pr) :: x_l(np, size(z)), w(size(z)), betas(np), P, T

      integer :: i, lb, ub

      integer :: ns
      real(pr) :: dS, S

      real(pr) :: X0(size(X))

      nc = size(z)

      number_of_points = optval(points, 1000)

      do i=1,np
         lb = (i-1)*nc + 1
         ub = i*nc
         X(lb:ub) = log(x_l0(i, :)/w0)
      end do

      X(np*nc + 1:np*nc + np) = betas0
      X(np*nc + np + 2) = log(T0)
      X(np*nc + np + 1) = log(P0)

      ns = ns0
      S = X(ns)
      dS = dS0

      allocate(env_points(0))
      do i=1,number_of_points
         X0 = X
         call solve_point(model, z, np, X, ns, S, dXdS, i, F, dF, its, max_iterations)

         ! The point might not converge, in this case we try again with an
         ! initial guess closer to the previous (converged) point.
         if (i < 1 .and. its > max_iterations) then
            X = X0 - 0.5 * dXdS * dS
            S = X(ns)
            call solve_point(model, z, np, X, ns, S, dXdS, i, F, dF, its, max_iterations)
         end if

         ! If the point did not converge, stop the calculation
         if (any(isnan(F)) .or. its > max_iterations) exit

         ! Save the information of the converged point
         call get_values_from_X(X, np, z, x_l, w, betas, P, T)
         point = MPPoint(np=np, nc=nc, betas=betas, P=P, T=T, x_l=x_l, w=w)
         env_points = [env_points, point]

         ! Update the specification for the next point.
         call update_specification(its, nc, np, X, dF, dXdS, ns, dS)

         ! Check if the system is close to a critical point, and try to jump
         ! over it.
         call detect_critical(nc, np, X, dXdS, ns, dS, S)

         ! Next point estimation.
         dX = dXdS * dS
         X = X + dX
         S = X(ns)
      end do

      call move_alloc(env_points, pt_envelope%points)
   end function pt_envelope

   subroutine pt_F_NP(model, z, np, x, ns, S, F, dF)
      !! Function to solve at each point of a multi-phase envelope.
      use iso_fortran_env, only: error_unit
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z(:) !! Mixture global composition.
      integer, intent(in) :: np !! Number of main phases.
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
      T = exp(X(np*nc + np + 1))
      P = exp(X(np*nc + np + 2))
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
         w, P, T, V=Vw, root_type="stable", lnphi=lnphi_w, &
         dlnphidp=dlnphi_dp_w, dlnphidt=dlnphi_dt_w, dlnphidn=dlnphi_dn_w &
         )

      do l=1,np
         x_l(l, :) = K(l, :)*w
         call model%lnphi_pt(&
            x_l(l, :), P, T, V=Vl(l), root_type="stable", lnphi=lnphi, &
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

                  df(idx_1, idx_2) = dlnphi_dn_l(phase, i, j) * dx_l_dlnK(phase, l, j) - dlnphi_dn_w(i, j) * dwdlnK(l, j)

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
         do phase=1,np
            do i=1,nc
               lb = (phase-1)*nc + i
               df(lb, nc*np+np+1) = P*(dlnphi_dp_l(phase, i) - dlnphi_dp_w(i))
               df(lb, nc*np+np+2) = T*(dlnphi_dt_l(phase, i) - dlnphi_dt_w(i))
            end do
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

   subroutine solve_point(model, z, np, X, ns, S, dXdS, point, F, dF, iters, max_iterations)
      use iso_fortran_env, only: error_unit
      use yaeos__math, only: solve_system
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z(:)
      integer, intent(in) :: np !! Number of main phases
      real(pr), intent(in out)  :: X(:) !! Vector of variables
      integer, intent(in)  :: ns !! Number of specification
      real(pr), intent(in)  :: S !! Specification value
      real(pr), intent(in) :: dXdS(size(X))
      integer, intent(in) :: point !! Number of point being calculated
      real(pr), intent(out) :: F(size(X)) !! Vector of functions valuated
      real(pr), intent(out) :: df(size(X), size(X)) !! Jacobian matrix
      integer, intent(in) :: max_iterations
      integer, intent(out) :: iters


      real(pr) :: X0(size(X))
      real(pr) :: dX(size(X))

      X0 = X

      do iters=1,max_iterations
         call pt_F_NP(model, z, np, x, ns, S, F, dF)

         dX = solve_system(dF, -F)
         ! print *, maxval(abs(F)), maxval(abs(dX)), exp(X(size(X)-1:))

         do while(maxval(abs(dX)) > 5)
            dX = dX/2
         end do

         if (maxval(abs(F)) < 1e-7_pr) exit

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

         if (maxval(abs(X(lb:ub))) < 0.05) then
            ns = lb + maxloc(abs(X(lb:ub)), dim=1)
            exit
         end if
      end do


      dS = dXdS(ns)*dS
      dXdS = dXdS/dXdS(ns)

      ! We adapt the step size to the number of iterations, the desired number
      ! of iterations for each point is around 3.
      dS = dS * 3._pr/its
   end subroutine update_specification

   subroutine detect_critical(nc, np, X, dXdS, ns, dS, S)
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

      do i=1,np
         lb = (i-1)*nc + 1
         ub = i*nc

         do while(maxval(abs(X(lb:ub))) < 0.01)
            X = X + 0.1 * dXdS * dS
         end do
      end do
   end subroutine detect_critical

   subroutine get_values_from_X(X, np, z, x_l, w, betas, P, T)
      !! # get_values_from_X
      !! Extract the values of the variables from the vector X.
      !!
      real(pr), intent(in) :: X(:) !! Vector of variables.
      integer, intent(in) :: np !! Number of main phases.
      real(pr), intent(in) :: z(:) !! Mixture composition.
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
      w = z/matmul(betas, x_l)

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
