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

      ! Incipient phase variables
      real(pr) :: Vw
      real(pr), dimension(size(z)) :: w, lnphi_w, dlnphi_dt_w, dlnphi_dp_w
      real(pr), dimension(size(z), size(z)) :: dlnphi_dn_w

      ! Derivative of w wrt beta
      real(pr) :: dwdb((Size(X)-3)/2)

      real(pr) :: dwdKx((Size(X)-3)/2), dxdKx((Size(X)-3)/2), dydKx((Size(X)-3)/2)

      real(pr) :: denom

      integer :: i, j, l, nc
      integer :: lb, ub

      nc = size(z)

      T = exp(X(np*nc + np + 1))
      P = exp(X(np*nc + np + 2))

      do l=1,np
         lb = (l-1)*nc + 1
         ub = l*nc
         K(l, :) = exp(X(lb:ub))
      end do

      betas = X(np*nc + 1:np*nc + np)

      denom = 0
      do i=1,np
         denom = denom + sum(betas(i)*K(i, :))
      end do


      P = exp(X(np*nc + np + 1))
      T = exp(X(np*nc + np + 2))

      w = z/denom

      call model%lnphi_pt(&
         w, P, T, V=Vw, root_type="stable", lnphi=lnphi_w, &
         dlnphidp=dlnphi_dp_w, dlnphidt=dlnphi_dt_w, dlnphidn=dlnphi_dn_w &
         )

      do l=1,np
         x_l(l, :) = K(l, :)*w

         call model%lnphi_pt(&
            x_l(l, :), P, T, V=Vl(l), root_type="stable", lnphi=lnphi_l(l, :), &
            dlnphidp=dlnphi_dp_l(l, :), dlnphidt=dlnphi_dt_l, dlnphidn=dlnphi_dn_l(l, :, :) &
            )
      end do

      F = 0
      df = 0

      do l=1,np
         lb = (l-1)*nc + 1
         ub = l*nc
         F(lb:ub) = X(lb:ub) + lnphi_l(l, :) - lnphi_w

         F(nc * np + l) = sum(x_l(l, :) - w )
      end do

      F(nc * np + np + 1) = sum(betas) - 1

      F(nc * np + np + 2) = X(ns) - S
   end subroutine pt_F_NP
end module yaeos__equilibria_boundaries_phase_envelopes_mp
