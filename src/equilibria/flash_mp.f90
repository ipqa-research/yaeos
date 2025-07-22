module yaeos__equilibria_multiphase_flash
   !! # `yaeos__equilibria_multiphase_flash`
   !! Module for multiphase flash calculations.
   !!
   !! # Description
   !! This module contains routines and functions to perform multiphase flash
   !! calculations.
   !!
   !! # Examples
   !!
   !! ```fortran
   !!  type(MPEquilibriumState) :: mpfr
   !!  type(ArModel) :: model
   !!  real(pr) :: z(3), P, T, Tc(3), Pc(3), w(3)
   !!  Tc = [374, 31, -83] + 273
   !!  Pc = [221, 74, 46]
   !!  w = [0.344, 0.293, 0.011]
   !!  z = [0.03_pr, 1-0.13_pr, 0.1_pr]
   !!  P = 45.6_pr
   !!  T = 190._pr
   !!  mpfr = pt_mp_flash(model, z, P, T)
   !!  print *, "Number of phases:", mpfr%np
   !!  print *, "Phase compositions:"
   !!  do i=1, mpfr%np
   !!     print *, "Phase", i, "composition:", mpfr%x_l(i, :)
   !!  end do
   !!  print *, "Reference phase composition:", mpfr%w(:)
   !! ```
   !!
   !! # References
   !!
   use yaeos__constants, only: pr, R
   use yaeos__models, only: ArModel
   implicit none

   type :: MPEquilibriumState
      !! # `MPEquilibriumState`
      !! Type to hold the state of a multiphase equilibrium calculation.
      !!
      !! # Description
      !! This type holds the results of a multiphase equilibrium calculation,
      !! including phase compositions, pressures, and temperatures.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! ```
      !!
      !! # References
      !!
      real(pr), allocatable :: z(:) !! Global composition
      real(pr) :: P !! Pressure
      real(pr) :: T !! Temperature
      integer :: np !! Number of phases
      real(pr), allocatable :: x_l(:,:)  !! Mole fractions of the main phases
      real(pr), allocatable :: w(:)  !! Mole fractions of the reference phase
      real(pr), allocatable :: betas(:)  !! Mole fractions of each phase
      character(len=14), allocatable :: kinds_x(:)  !! Kinds of the main phases
      character(len=14) :: kind_w !! Kind of the reference phase
   end type MPEquilibriumState

contains

   type(MPEquilibriumState) function pt_mp_flash(model, z, P, T)
      !! # `pt_mp_flash`
      !! Perform a multiphase flash calculation at constant zPT.
      !!
      !! # Description
      !! This method will do stability analysis to detect the possibility of
      !! new phases. For each new phase detected it will calculate a multiphase
      !! flash and repeat stability analysis until no new phases are detected.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! ```
      !!
      !! # References
      !!
      use yaeos__equilibria_stability, only: min_tpd
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z(:)
      real(pr), intent(in) :: P, T
      integer, parameter :: max_phases = 4
      integer :: np

      real(pr) :: mintpd, all_minima(size(z), size(z)), w(size(z)), w_stab(size(z))
      real(pr) :: mintpd_xl1, mintpd_w
      real(pr) :: x_l(max_phases, size(z))
      real(pr) :: K(max_phases, size(z))

      character(len=14) :: kinds_x(max_phases)
      character(len=14) :: kind_w
      logical :: less_phases

      integer :: beta_0_index, iters, ns1, ns2, nc
      real(pr) :: S1, S2
      integer :: max_iters

      real(pr), allocatable :: X(:), F(:), betas(:)
      real(pr) :: beta0

      max_iters = 1000
      kinds_x = "liquid"
      kind_w = "vapor"

      nc = size(z)
      np = 0
      S1 = log(P)
      S2 = log(T)

      call min_tpd(model, z, P, T, mintpd, w_stab)
      if (mintpd < -0.01) then
         np = np + 1

         ns1 = np*nc + np + 1 + 1
         ns2 = np*nc + np + 1 + 2

         x_l(1, :) = z
         K(1, :) = x_l(1, :) / w_stab

         beta0 = z(maxloc(w_stab, dim=1))
         beta0 = 0.0001
         betas = [1-beta0, beta0]

         X = [log(K(1, :)), betas, log(P), log(T)]
         F = X
         call solve_mp_flash_point(&
            model, z, np, kinds_x, kind_w, X, ns1, S1, ns2, S2, max_iters, F, &
            less_phases, beta_0_index, iters &
            )

         K(1, :) = exp(X(:nc))

         betas = X(np*nc+1 : np*nc+np+1)

         w = z/(matmul(betas(:np), K(:np, :)) + betas(np+1))
         x_l(1, :) = K(1, :) * w

         call min_tpd(model, w, P, T, mintpd_w, w_stab)
         call min_tpd(model, x_l(1, :), P, T, mintpd_xl1, w_stab)
         if (mintpd_w < -0.001) then
            np = np + 1

            ns1 = np*nc + np + 1 + 1
            ns2 = np*nc + np + 1 + 2

            x_l(2, :) = w
            w = w_stab

            K(1, :) = x_l(1, :) / w
            K(2, :) = x_l(2, :) / w

            betas = [betas, 0.5_pr]

            X = [log(K(1, :)), log(K(2, :)), betas, log(P), log(T)]
            F = X
            call solve_mp_flash_point(&
               model, z, np, kinds_x, kind_w, X, ns1, S1, ns2, S2, max_iters, F, &
               less_phases, beta_0_index, iters &
               )

            K(1, :) = exp(X(:nc))
            K(2, :) = exp(X(nc+1:2*nc))

            betas = X(np*nc+1 : np*nc+np+1)
            w = z/(matmul(betas(:np), K(:np, :)) + betas(np+1))
            x_l(1, :) = K(1, :) * w
            x_l(2, :) = K(2, :) * w
         end if
      end if

      pt_mp_flash%z = z
      pt_mp_flash%P = P
      pt_mp_flash%T = T
      pt_mp_flash%np = np
      pt_mp_flash%x_l = x_l(1:np, :)
      pt_mp_flash%w = w
      pt_mp_flash%kinds_x = kinds_x(1:np)
      pt_mp_flash%kind_w = kind_w
      pt_mp_flash%betas = betas


   end function pt_mp_flash

   subroutine pt_F_NP(model, z, np, kinds_x, kind_w, X, ns1, S1, ns2, S2, F, dF)
      !! Function to solve at each point of a multi-phase envelope.
      use iso_fortran_env, only: error_unit
      class(ArModel), intent(in) :: model !! Model to use.
      real(pr), intent(in) :: z(:) !! Mixture global composition.
      integer, intent(in) :: np !! Number of main phases.
      character(len=14), intent(in) :: kinds_x(np) !! Kind of the main phases.
      character(len=14), intent(in) :: kind_w !! Kind of the reference phase.
      real(pr), intent(in)  :: X(:) !! Vector of variables.
      integer, intent(in)  :: ns1 !! Number of first specification.
      real(pr), intent(in)  :: S1 !! First specification value.
      integer, intent(in)  :: ns2 !! Number of second specification.
      real(pr), intent(in)  :: S2 !! Second specification value.
      real(pr), intent(out) :: F(size(X)) !! Vector of functions valuated.
      real(pr), intent(out) :: df(size(X), size(X)) !! Jacobian matrix.

      ! X variables
      real(pr) :: K(np, size(z))
      real(pr) :: P
      real(pr) :: T
      real(pr) :: betas(np), beta_w

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
      real(pr) :: dwdb(np, size(z)), dwdbw(size(z))
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
      beta_w = X(np*nc + np + 1)
      P = exp(X(np*nc + np + 2))
      T = exp(X(np*nc + np + 3))

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
      F(nc * np + np + 2) = X(ns1) - S1
      F(nc * np + np + 3) = X(ns2) - S2


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

      dwdbw = -z / denom**2

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
            df(lb, nc*np+np+2) = P*(dlnphi_dp_l(l, i) - dlnphi_dp_w(i))
            df(lb, nc*np+np+3) = T*(dlnphi_dt_l(l, i) - dlnphi_dt_w(i))
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

      do j=1,np
         lb = nc*np + j
         df(lb,np*nc+np+1) = sum(K(j, :) * dwdbw(:) - dwdbw(:))
      end do

      do j=1,np
         lb = (j-1)*nc + 1
         ub = j*nc
         do i=1,nc
            df(lb+i-1, np*nc + np + 1) = &
               sum(K(j, :) * dlnphi_dn_l(j, i, :)*dwdbw(:) &
               - dlnphi_dn_w(i, :)*dwdbw(:))
         end do
      end do


      df(nc * np + np + 1, np*nc + np + 1) = 1

      df(nc * np + np + 2, ns1) = 1
      df(nc * np + np + 3, ns2) = 1
   end subroutine pt_F_NP

   subroutine solve_mp_flash_point(&
      model, z, np, kinds_x, kind_w, X, ns1, S1, ns2, S2, max_iters, F, &
      less_phases, beta_0_index, iters &
      )
      !! Function to solve the multiphase flash problem.
      use yaeos__math, only: solve_system
      class(ArModel), intent(in) :: model !! Model to use.
      real(pr), intent(in) :: z(:) !! Mixture global composition.
      integer, intent(in) :: np !! Number of x phases.
      character(len=14), intent(in) :: kinds_x(np) !! Kind of the x phases.
      character(len=14), intent(in) :: kind_w !! Kind of the w phase.
      real(pr), intent(inout) :: X(:) !! Vector of variables.
      integer, intent(in) :: ns1 !! Number of first specification.
      real(pr), intent(in) :: S1 !! First specification value.
      integer, intent(in) :: ns2 !! Number of second specification.
      real(pr), intent(in) :: S2 !! Second specification value.
      integer, intent(in) :: max_iters !! Maximum number of iterations.
      real(pr), intent(out) :: F(size(X)) !! Vector of functions valuated.
      logical, intent(out) :: less_phases !! True if the solution has less phases than expected.
      integer, intent(out) :: beta_0_index !! Index of beta that equals zero.
      integer, intent(out) :: iters !! Number of iterations performed.

      real(pr), dimension(size(X), size(X)) :: df !! Jacobian matrix.
      real(pr), dimension(size(X)) :: dX !! Newton step vector.

      integer :: nc !! Number of components
      integer :: iBetas(np+1) !! Index of the betas in the vector X
      integer :: i

      nc = size(z)
      iBetas = [(i, i=nc*np+1,nc*np+np+1)]

      less_phases = .false.

      do iters=1,max_iters
         call pt_F_NP(model, z, np, kinds_x, kind_w, X, ns1, S1, ns2, S2, F, df)

         if (maxval(abs(F)) < 1e-10) exit

         dX = solve_system(df, -F)

         do while(maxval(abs(dX)) > 1)
            dX = dX/2
         end do

         do while(any(abs(X(iBetas) + dX(iBetas)) > 1))
            X(iBetas) = X(iBetas) / 2
         end do

         X = X + dX
      end do
   end subroutine solve_mp_flash_point
end module yaeos__equilibria_multiphase_flash
