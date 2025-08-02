module yaeos__equilibria_boundaries_generalized_isopleths
   !! Calculation of isoplethic phase equilibria lines.
   !!
   !! This module contains the subroutines to calculate any kind of phase
   !! equilibria lines with constant composition.
   use yaeos__constants, only: pr, R
   use yaeos__models, only: ArModel
   use yaeos__equilibria_equilibrium_state, only: MPEquilibriumState, MPEquilibriumState_from_X
   use yaeos__equilibria_stability, only: tm
   use yaeos__math, only: solve_system
   implicit none

   type :: GeneralizedIsoZLine
      type(MPEquilibriumState), allocatable :: points(:)
      logical :: did_stability = .false.
      logical :: found_unstability = .false.
      real(pr), allocatable :: w_more_stable(:)
   end type GeneralizedIsoZLine

contains
   type(GeneralizedIsoZLine) function create_generalized_isoz_line(&
      model, nc, np, nstab, kinds_x, kind_w, z, x_l0, w0, betas0, P0, T0, &
      spec_variable, spec_variable_value, ns0, S0, dS0, &
      ws_stab, max_points &
      )
      !! Create a new generalized line.
      !! This function initializes a new instance of the GeneralizedLine type.
      class(ArModel), intent(in) :: model
      integer, intent(in) :: nc
      integer, intent(in) :: np
      integer, intent(in) :: nstab
      character(len=14), intent(in) :: kinds_x(np)
      character(len=14), intent(in) :: kind_w
      real(pr), intent(in) :: z(nc)
      real(pr), intent(in) :: x_l0(np, nc)
      real(pr), intent(in) :: w0(nc)
      real(pr), intent(in) :: betas0(np+1)
      real(pr), intent(in) :: P0
      real(pr), intent(in) :: T0
      integer, intent(in) :: spec_variable
      real(pr), intent(in) :: spec_variable_value
      integer, intent(in) :: ns0
      real(pr), intent(in) :: S0
      real(pr), intent(in) :: dS0
      real(pr), optional, intent(in) :: ws_stab(nstab, nc)
      integer, optional, intent(in) :: max_points

      real(pr) :: tms(nstab)

      integer :: ns, ns1
      real(pr) :: S, S1, dS
      real(pr) :: X(nc*np + np + 3)
      real(pr) :: F(nc*np + np + 3), df(nc*np + np + 3, nc*np + np + 3)
      real(pr) :: dX(nc*np + np + 3), dXdS(nc*np+np+3), dFdS(nc*np+np+3)

      integer :: lb, ub, i, l, i_point
      integer :: iT, iP, iBetas(np+1)
      character(len=14) :: x_kinds(np), w_kind
      integer :: iters, max_iters = 1000
      type(MPEquilibriumState) :: point

      logical :: found_unstability

      ns = ns0
      S = S0
      dS = dS0

      found_unstability = .false.

      ibetas = [(i, i=nc*np+1, nc*np+np+1)]

      dFdS = 0
      dFdS(nc*np+np+3) = -1

      X = [log([(x_l0(i, :)/w0, i=1,np)]), betas0, log(P0), log(T0)]

      i_point = 0
      allocate(create_generalized_isoz_line%points(0))
      do while(iters < max_iters .and. .not. found_unstability)
         i_point = i_point + 1
         call solve_generalized_point(model, z, np, kinds_x, kind_w, &
            X, spec_variable, spec_variable_value, ns, S, max_iters, F, dF, &
            iters)

         dXdS = solve_system(dF, -dFdS)
         ns = maxloc(dXdS, dim=1)

         dS = dXdS(ns) * dS * 3./iters
         dXdS = dXdS/dXdS(ns)

         point = MPEquilibriumState_from_X(nc, np, z, kinds_x, kind_w, X)
         
         do i=1,nstab
            tms(i) = tm(model, point%w, ws_stab(i, :), point%P, point%T)
         end do

         create_generalized_isoz_line%points = [&
            create_generalized_isoz_line%points, point &
         ]

         if (any(tms < -0.01)) then
            i = minloc(tms, dim=1)
            create_generalized_isoz_line%w_more_stable = ws_stab(i,:)
            found_unstability = .true.
            create_generalized_isoz_line%found_unstability = found_unstability
         end if

         if (present(max_points)) then
            if (i_point > max_points) exit
         end if

         dX = dXdS * dS
         X = X + dX
         S = X(ns)
      end do

   end function create_generalized_isoz_line

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

   subroutine solve_generalized_point( &
      model, z, np, kinds_x, kind_w, X, ns1, S1, ns2, S2, max_iters, F, dF, &
      iters &
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
      real(pr), intent(out), dimension(size(X), size(X)) :: df !! Jacobian matrix.
      integer, intent(out) :: iters !! Number of iterations performed.

      real(pr), dimension(size(X)) :: dX !! Newton step vector.

      integer :: nc !! Number of components
      integer :: iBetas(np+1) !! Index of the betas in the vector X
      integer :: i

      nc = size(z)
      iBetas = [(i, i=nc*np+1,nc*np+np+1)]

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
   end subroutine solve_generalized_point
end module yaeos__equilibria_boundaries_generalized_isopleths
