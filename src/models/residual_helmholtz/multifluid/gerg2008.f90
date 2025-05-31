module yaeos__models_ar_gerg2008
   use yaeos__constants, only: Ryaeos => R, pr !! Ideal gas constants used on yaeos
   use yaeos__adiff_hyperdual_ar_api, only: ArModelAdiff
   use yaeos__models_ar_cubic_implementations, only: SoaveRedlichKwong
   use yaeos__models_ar_genericcubic, only: CubicEoS
   use yaeos__models_ar_multifluid_parameters_gerg2008, only: Gerg2008Binary, Gerg2008Pure

   use hyperdual_mod

   implicit none

   type, extends(ArModelAdiff) :: Gerg2008
      type(Gerg2008Pure), allocatable :: pures(:)
      type(Gerg2008Binary), allocatable :: binaries(:, :)
      type(CubicEoS) :: srk
   contains
      procedure :: ar => arfun
      procedure :: get_v0 => volume_initalizer
   end type Gerg2008

   type, private :: GERG2008Selector
      integer :: methane=1
      integer :: nitrogen=2
      integer :: carbon_dioxide=3
      integer :: ethane=4
      integer :: propane=5
      integer :: nbutane=6
      integer :: isobutane=7
      integer :: npentane=8
      integer :: isopentane=9
      integer :: nhexane=10
      integer :: nheptane=11
      integer :: noctane=12
      integer :: nonane=13
      integer :: decane=14
      integer :: hydrogen=15
      integer :: oxygen=16
      integer :: carbon_monoxide=17
      integer :: water=18
      integer :: hydrogen_sulfide=19
      integer :: helium=20
      integer :: argon=21
   end type GERG2008Selector

   type(GERG2008Selector) :: G2008Components

contains

   type(Gerg2008) function gerg_2008(ids)
      use yaeos__models_ar_multifluid_parameters_gerg2008, only: get_original_parameters
      integer, intent(in) :: ids(:)
      type(Gerg2008Pure) :: pures(size(ids))
      type(Gerg2008Binary) :: binaries(size(ids), size(ids))

      call get_original_parameters(ids, pures, binaries, gerg_2008%components)
      gerg_2008%pures = pures
      gerg_2008%binaries = binaries
      gerg_2008%srk =SoaveRedlichKwong(gerg_2008%components%Tc, gerg_2008%components%Pc, gerg_2008%components%w)
   end function gerg_2008

   subroutine reducing_functions(self, n, Vr, Tr)
      class(Gerg2008), intent(in) :: self
      type(hyperdual), intent(in) :: n(:)
      type(hyperdual), intent(out) :: Vr
      type(hyperdual), intent(out) :: Tr

      type(hyperdual) :: X(size(n))

      real(8) :: Vc(size(n)), Tc(size(n)), rho_c(size(n))

      real(8) :: Bv(size(n), size(n)), Gv(size(n), size(n))
      real(8) :: Bt(size(n), size(n)), Gt(size(n), size(n))

      integer :: i, j, nc

      Vc = self%components%Vc
      Tc = self%components%Tc
      Bv = self%binaries%Bv
      Gv = self%binaries%Gv
      Bt = self%binaries%Bt
      Gt = self%binaries%Gt

      rho_c = 1/Vc
      X = n / sum(n)
      nc = size(n)

      Vr = sum(X ** 2 * Vc)
      Tr = sum(X ** 2 * Tc)

      do i=1,nc
         do j=i+1,nc
            Vr = Vr + &
               2 * X(i) * X(j) * Bv(i, j) * Gv(i, j) &
               * (X(i) + X(j)) / (Bv(i, j) ** 2 * X(i) + X(j)) &
               * 1._pr / 8._pr * (rho_c(i) ** (- 1._pr / 3._pr) &
               + rho_c(j) ** (- 1._pr / 3)) ** 3

            Tr = Tr + &
               2 * X(i) * X(j) * Bt(i, j) * Gt(i, j) &
               * (X(i) + X(j)) / (Bt(i, j) ** 2 * X(i) + X(j)) &
               * sqrt((Tc(i) * Tc(j)))
         end do
      end do
   end subroutine reducing_functions

   subroutine ar_pure(pure, delta, tau, ar)
      type(Gerg2008Pure), intent(in) :: pure
      type(hyperdual), intent(in) :: delta
      type(hyperdual), intent(in) :: tau
      type(hyperdual), intent(out) :: ar

      integer :: i, Kpol, Kexp

      real(8) :: n_pol(pure%Kpol), d_pol(pure%Kpol), t_pol(pure%Kpol)
      real(8) :: n_exp(pure%Kexp), d_exp(pure%Kexp), t_exp(pure%Kexp)
      real(8) :: c_exp(pure%Kexp)

      Kpol = pure%Kpol
      Kexp = pure%Kexp

      n_pol = pure%n(1:Kpol)
      d_pol = pure%d(1:Kpol)
      t_pol = pure%t(1:Kpol)
      n_exp = pure%n(Kpol+1:Kpol+Kexp)
      d_exp = pure%d(Kpol+1:Kpol+Kexp)
      t_exp = pure%t(Kpol+1:Kpol+Kexp)

      c_exp = pure%c

      ar = sum(n_pol * delta ** d_pol * tau ** t_pol) + &
         sum(n_exp * delta**d_exp * tau**t_exp * exp(-delta**c_exp))
   end subroutine ar_pure

   subroutine ar_ij(delta, tau, binary, aij)
      type(hyperdual), intent(in) :: delta
      type(hyperdual), intent(in) :: tau
      type(Gerg2008Binary), intent(in) :: binary
      type(hyperdual), intent(out) :: aij

      integer :: idx_poly, idx_exp

      real(8) :: n_pol(binary%Kpolij), d_pol(binary%Kpolij), t_pol(binary%Kpolij)
      real(8) :: n_exp(binary%Kexpij), d_exp(binary%Kexpij), t_exp(binary%Kexpij)
      real(8) :: etha(binary%Kexpij), eps(binary%Kexpij), beta(binary%Kexpij), gama(binary%Kexpij)

      idx_poly = binary%Kpolij
      idx_exp = binary%Kexpij + idx_poly

      n_pol = binary%nij(1:idx_poly)
      d_pol = binary%dij(1:idx_poly)
      t_pol = binary%tij(1:idx_poly)

      n_exp = binary%nij(idx_poly+1:idx_exp)
      d_exp = binary%dij(idx_poly+1:idx_exp)
      t_exp = binary%tij(idx_poly+1:idx_exp)

      etha = binary%ethaij(1:binary%Kexpij)
      eps = binary%epsij(1:binary%Kexpij)
      beta = binary%betaij(1:binary%Kexpij)
      gama = binary%gammaij(1:binary%Kexpij)

      aij = sum(n_pol * delta ** d_pol * tau ** t_pol) + &
         sum(n_exp * delta**d_exp * tau**t_exp * exp(-etha * (delta - eps) ** 2 - beta * (delta - gama)))
   end subroutine ar_ij

   function arfun(self, n, v, t) result(arval)
      class(Gerg2008) :: self
      type(hyperdual), intent(in) :: n(:), v, t
      type(hyperdual) :: arval

      type(hyperdual) :: Vr, Tr, X(size(n)), rho_r
      type(hyperdual) :: delta, tau
      type(hyperdual) :: aij
      type(hyperdual) :: ar_pures(size(n))

      type(Gerg2008Pure) :: pures(size(n))
      type(Gerg2008Binary) :: binary

      real(pr) :: rho_c(size(n))
      real(pr) :: Fij(size(n), size(n))
      integer :: i, j, nc

      Fij = self%binaries%Fij

      pures = self%pures

      nc = size(n)
      X = n / sum(n)
      call reducing_functions(self, n, Vr, Tr)

      rho_r = 1._pr/Vr

      delta = (1._pr/(V/sum(n)))/rho_r
      tau = Tr/T

      do i=1,nc
         call ar_pure(pures(i), delta, tau, ar_pures(i))
      end do

      arval = sum(x * ar_pures)

      do i=1,nc
         do j=1,nc!i+1,nc
            if (Fij(i, j) == 0._pr) cycle
            binary = self%binaries(i, j)
            call ar_ij(delta, tau, binary, aij)
            arval = arval + X(i) * X(j) * Fij(i, j) * aij
         end do
      end do
      arval = arval * (sum(n) * ryaeos * t)
   end function arfun

   function volume_initalizer(self, n, p, t) result(v0)
      class(Gerg2008), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: p
      real(pr), intent(in) :: t
      real(pr) :: v0
      v0 = self%srk%get_v0(n, p, t)
   end function volume_initalizer
end module yaeos__models_ar_gerg2008