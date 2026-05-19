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
      integer, allocatable :: ids(:)
      type(CubicEoS) :: srk
   contains
      procedure :: ar => arfun
      procedure :: get_v0 => volume_initalizer
      procedure :: volume => volume
      procedure :: ln_activity_coefficient => ln_activity_coefficient
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
      gerg_2008%ids = ids
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

      real(pr) :: n_pol(binary%Kpolij), d_pol(binary%Kpolij), t_pol(binary%Kpolij)
      real(pr) :: n_exp(binary%Kexpij), d_exp(binary%Kexpij), t_exp(binary%Kexpij)
      real(pr) :: etha(binary%Kexpij), eps(binary%Kexpij), beta(binary%Kexpij), gama(binary%Kexpij)

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

   subroutine volume(eos, n, P, T, V, root_type)
      !! Volume solver routine for the GERG2008.
      !!
      !! Solves volume roots using newton method. Given pressure and
      !! temperature. It will use the SRK equation of state to initialize the
      !! volume values.
      !!
      !! # Description
      !! This subroutine solves the volume using a newton method. The variable
      !! `root_type` is used to specify the desired root to solve. The options
      !! are: `["liquid", "vapor", "stable"]`
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! P = 1.0_pr
      !!
      !! call eos%volume(n, P, T, V, root_type="liquid")
      !! call eos%volume(n, P, T, V, root_type="vapor")
      !! call eos%volume(n, P, T, V, root_type="stable")
      !! ```
      use yaeos__constants, only: pr, R
      use yaeos__math, only: newton

      class(GERG2008), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(out) :: V !! Volume [L]
      character(len=*), intent(in) :: root_type
      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]`

      integer :: max_iters=500
      real(pr) :: tol=1e-8

      real(pr) :: totnRT, GrL, GrV, Gr
      real(pr) :: Vliq, Vvap
      logical :: failed

      GrL = HUGE(GrL)
      GrV = HUGE(GrV)

      totnRT = sum(n) * R * T
      select case(root_type)
       case("liquid")
         call eos%srk%volume(n, P=P, T=T, V=Vliq, root_type="liquid")
         Vliq = log(Vliq)
         call newton(foo, Vliq, tol=tol, max_iters=max_iters, failed=failed)
         Vliq = exp(Vliq)
         GrL = Gr
       case("vapor")
         call eos%srk%volume(n=n, P=P, T=T, V=Vvap, root_type="vapor")
         Vvap = log(Vvap)
         call newton(foo, Vvap, tol=tol, max_iters=max_iters, failed=failed)
         Vvap = exp(Vvap)
         GrV = Gr
       case("stable")
         call eos%srk%volume(n, P=P, T=T, V=Vliq, root_type="liquid")
         Vliq = log(Vliq)
         call newton(foo, Vliq, tol=tol, max_iters=max_iters, failed=failed)
         Vliq = exp(Vliq)
         GrL = Gr

         call eos%srk%volume(n, P=P, T=T, V=Vvap, root_type="vapor")
         Vvap = log(Vvap)
         call newton(foo, Vvap, tol=tol, max_iters=max_iters, failed=failed)
         Vvap = exp(Vvap)
         GrV = Gr
      end select

      if (GrL < GrV) then
         V = Vliq
      else
         V = Vvap
      end if

      if (failed) V = -1

   contains
      subroutine foo(x, f, df)
         real(pr), intent(in) :: x
         real(pr), intent(out) :: f, df
         real(pr) :: Ar, ArV, ArV2, Pcalc, dPcalcdV, Vin
         Vin = exp(x)
         call eos%residual_helmholtz(n, Vin, T, Ar=Ar, ArV=ArV, ArV2=ArV2)
         Pcalc = totnRT / Vin - ArV
         dPcalcdV = -totnRT / Vin**2 - ArV2
         f = Pcalc - P
         df = Vin * dPcalcdV
         Gr = Ar + P * Vin - totnRT - totnRT * log(P*Vin/(R*T))
      end subroutine foo
   end subroutine volume

   function volume_initalizer(self, n, p, t) result(v0)
      class(Gerg2008), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: p
      real(pr), intent(in) :: t
      real(pr) :: v0
      v0 = self%srk%get_v0(n, p, t)
   end function volume_initalizer

   subroutine ln_activity_coefficient(&
      eos, n, P, T, root_type, lngamma, dlngammadP, dlngammadT, dlngammadn &
      )
      !! Calculate natural logarithm of activity coefficients and its
      !! derivatives given pressure and temperature.
      !!
      !! # Examples
      !!
      !! ```fortran ! eos = PengRobinson76(Tc, Pc, w)
      !!
      !! n = [1.0_pr, 1.0_pr] ! T = 300.0_pr ! P = 1.0_pr
      !!
      !! call eos%ln_activity_coefficient(&
      !!    n, P, T, root_type="stable", &
      !!    lngamma=lngamma, dlngammadP=dlngammadP, &
      !!    dlngammadT=dlngammadT, dlngammadn=dlngammadn &
      !!    )
      !! ```
      class(GERG2008), intent(in) :: eos !! Model
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      character(len=*), intent(in) :: root_type
      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]`
      real(pr), intent(out), optional :: lngamma(size(n))
      !! Natural logarithm of activity coefficient
      real(pr), intent(out), optional :: dlngammadP(size(n))
      !! \(\frac{d\ln\gamma}}{dP}\)
      real(pr), intent(out), optional :: dlngammadT(size(n))
      !! \(\frac{d\ln\gamma}}{dT}\)
      real(pr), intent(out), optional :: dlngammadn(size(n),size(n))
      !! \(\frac{d\ln\gamma}}{dn}\)

      ! Mixture properties
      real(pr) :: lnPhi(size(n)), dlnPhidT(size(n))
      real(pr) :: dlnPhidn(size(n),size(n))
      real(pr) :: dPdV, dPdn(size(n)), dVdn(size(n))

      ! Pure properties
      real(pr) :: vi(size(n)), vi_temp, npure(size(n))
      real(pr) :: lnPhi_i(size(n)), dlnPhi_i_dT(size(n))
      real(pr) :: lnPhi_i_temp(1), dlnPhi_i_dT_temp(1)
      type(GERG2008) :: eos_temp

      integer :: i

      logical :: gam, dp, dt, dn, present_derivs

      gam = present(lngamma)
      dp = present(dlngammadP)
      dt = present(dlngammadT)
      dn = present(dlngammadn)
      present_derivs = dp .or. dt .or. dn


      ! Call to mixture fugacity coefficient at PT
      ! TODO: maybe more efficient later?
      if (.not. present_derivs) then
         call eos%lnphi_pt(n=n, P=P, T=T, root_type=root_type, lnPhi=lnPhi)
      else
         call eos%lnphi_pt(&
            n=n, P=P, T=T, root_type=root_type, &
            lnPhi=lnPhi, dlnPhidT=dlnPhidT, dlnPhidn=dlnPhidn, &
            dPdV=dPdV, dPdn=dPdn &
            )

         dVdn = -dPdn / dPdV
      end if

      ! Pure components calls
      vi = 0.0_pr
      lnPhi_i = 0.0_pr

      if (dt) then
         dlnPhi_i_dT = 0.0_pr
      end if

      do i=1, size(n)
         eos_temp = gerg_2008([eos%ids(i)])

         if (.not. dt) then
            call eos_temp%lnphi_pt(&
               n=[1.0_pr], P=P, T=T, V=vi_temp, root_type="stable", &
               lnPhi=lnPhi_i_temp &
               )
         else
            call eos_temp%lnphi_pt(&
               n=[1.0_pr], P=P, T=T, V=vi_temp, root_type="stable", &
               lnPhi=lnPhi_i_temp, dlnPhidT=dlnPhi_i_dT_temp &
               )
         end if

         vi(i) = vi_temp
         lnPhi_i(i) = lnPhi_i_temp(1)

         if (dt) then
            dlnPhi_i_dT(i) = dlnPhi_i_dT_temp(1)
         end if
      end do

      ! returns
      if (gam) lngamma = lnPhi - lnPhi_i
      if (dt) dlngammadT = dlnPhidT - dlnPhi_i_dT
      if (dp) dlngammadP = (dVdn - vi) / (Ryaeos * T)
      if (dn) dlngammadn = dlnPhidn
   end subroutine ln_activity_coefficient
end module yaeos__models_ar_gerg2008
