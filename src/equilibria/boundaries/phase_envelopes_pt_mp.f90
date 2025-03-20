module yaeos__equilibria_boundaries_phase_envelopes_mp
   use yaeos__constants, only: pr, R
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState
   use yaeos__models_ar, only: ArModel
   use yaeos__math, only: solve_system

   implicit none

   private

   public :: PTEnvelMP
   public :: pt_F_NP
   ! public :: pt_envelope
   ! public :: solve_point
   ! public :: get_values_from_X

   type :: PTEnvelMP
      integer :: NP !! Number of phases, besides the incipient phase
      integer :: nc !! Number of components
      integer, allocatable :: its(:) !! Number of needed iterations
      real(pr), allocatable :: beta(:, :) !! Phases mole fractions
      real(pr), allocatable :: P(:) !! Pressures [bar]
      real(pr), allocatable :: T(:) !! Temperatures [K]
      integer, allocatable :: ns(:) !! Number of specified variable
      real(pr), allocatable :: S(:) !! Value of specification
   end type PTEnvelMP

   real(pr), parameter :: lnK_min = 2.0_pr

contains

   subroutine pt_F_NP(model, z, np, x, ns, S, F, dF)
      !! Function to solve at each point of a multi-phase envelope.
      use iso_fortran_env, only: error_unit
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z(:)
      integer, intent(in) :: np !! Number of main phases
      real(pr), intent(in)  :: X(:) !! Vector of variables
      integer, intent(in)  :: ns !! Number of specification
      real(pr), intent(in)  :: S !! Specification value
      real(pr), intent(out) :: F(size(X)) !! Vector of functions valuated
      real(pr), intent(out) :: df(size(X), size(X)) !! Jacobian matrix

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
      denom = matmul(betas, K)
      denomdlnK = 0
      do i=1,nc
         denomdlnK(:, i, i) = betas(:)*K(:, i)
      end do

      w = z/denom

      ! autodiff:block
      !    use hyperdual_mod
      !    type(hyperdual) :: hd_w(size(z)), hd_denom, lnKs(np, size(z)), hd_betas(np)
      !    type(hyperdual) :: hd_xl(np, size(z))

      !    do i=1,nc


      !       do l=1,np
      !          hd_betas = betas
      !          lnKs = log(K)
      !          lnKs(l, i)%f1 = 1

      !          hd_denom = 0._pr
      !          do j=1,np
      !             hd_denom = hd_denom + sum(hd_betas(l) * exp(lnKs(l, :)))
      !          end do
      !          hd_w = z/hd_denom

      !          hd_xl(l, :) = exp(lnKs(l, :)) * hd_w


      !          dwdlnK(l, :, i) = hd_w%f1
      !       end do
      !    end do

      ! end block autodiff

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
      F(nc * np + np + 1) = sum(betas) - 1
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
end module yaeos__equilibria_boundaries_phase_envelopes_mp
