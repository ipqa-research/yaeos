module yaeos__models_ar_gerg2008
   !! # GERG2008 equation of state
   !!
   !! # Description
   !! Multifluid equation of state for natural-gas mixtures based on the
   !! reduced Helmholtz energy. The model is constructed with
   !! [[gerg_2008]] from the component identifiers in
   !! [[G2008Components]].
   !! The model uses the composition-dependent reducing functions and binary
   !! departure terms of the original 21-component formulation.
   !!
   !! # Examples
   !! ```fortran
   !! model = gerg_2008([G2008Components%methane, G2008Components%nitrogen])
   !! call model%pressure(n, V, T, P)
   !! ```
   use yaeos__constants, only: Ryaeos => R, pr !! Ideal gas constants used on yaeos
   use yaeos__adiff_hyperdual_ar_api, only: ArModelAdiff
   use yaeos__models_ar_cubic_implementations, only: SoaveRedlichKwong
   use yaeos__models_ar_genericcubic, only: CubicEoS
   use yaeos__models_ar_multifluid_parameters_gerg2008, only: Gerg2008Binary, Gerg2008Pure

   use hyperdual_mod

   implicit none

   type, extends(ArModelAdiff) :: Gerg2008
      !! # Gerg2008
      !!
      !! GERG2008 multifluid equation of state. Holds the pure-fluid and
      !! binary parameters of the selected components and an SRK model used
      !! to initialize the volume solver.
      !!
      !! # Description
      !! Provides the standard residual Helmholtz properties through the
      !! inherited hyperdual adapter and specializes volume solving with an
      !! SRK initial estimate.
      !!
      !! # Examples
      !! ```fortran
      !! model = gerg_2008([G2008Components%methane, G2008Components%nitrogen])
      !! call model%pressure(n, V, T, P)
      !! ```
      type(Gerg2008Pure), allocatable :: pures(:) !! Pure-fluid Helmholtz parameters, in model component order.
      type(Gerg2008Binary), allocatable :: binaries(:, :) !! Pair interaction data indexed by component pair.
      integer, allocatable :: ids(:) !! GERG2008 identifiers, in model component order.
      type(CubicEoS) :: srk !! SRK helper used to initialize the volume root.
   contains
      procedure :: ar => arfun
      procedure :: get_v0 => volume_initalizer
      procedure :: volume => volume
      procedure :: ln_activity_coefficient => ln_activity_coefficient
   end type Gerg2008

   type, private :: GERG2008Selector
      !! # GERG2008Selector
      !! Internal integer codes for the 21 components supported by GERG2008.
      !!
      !! # Description
      !! The public `G2008Components` object exposes these code values as
      !! named component selectors for model construction.
      !!
      !! # Examples
      !! ```fortran
      !! ids = [G2008Components%methane, G2008Components%nitrogen]
      !! ```
      integer :: methane=1 !! Methane.
      integer :: nitrogen=2 !! Nitrogen.
      integer :: carbon_dioxide=3 !! Carbon dioxide.
      integer :: ethane=4 !! Ethane.
      integer :: propane=5 !! Propane.
      integer :: nbutane=6 !! n-Butane.
      integer :: isobutane=7 !! Isobutane.
      integer :: npentane=8 !! n-Pentane.
      integer :: isopentane=9 !! Isopentane.
      integer :: nhexane=10 !! n-Hexane.
      integer :: nheptane=11 !! n-Heptane.
      integer :: noctane=12 !! n-Octane.
      integer :: nonane=13 !! n-Nonane.
      integer :: decane=14 !! n-Decane.
      integer :: hydrogen=15 !! Hydrogen.
      integer :: oxygen=16 !! Oxygen.
      integer :: carbon_monoxide=17 !! Carbon monoxide.
      integer :: water=18 !! Water.
      integer :: hydrogen_sulfide=19 !! Hydrogen sulfide.
      integer :: helium=20 !! Helium.
      integer :: argon=21 !! Argon.
   end type GERG2008Selector

   type(GERG2008Selector) :: G2008Components
   !! Component identifiers for the GERG2008 model.

contains

   type(Gerg2008) function gerg_2008(ids)
      !! Create a GERG2008 model from a list of component identifiers.
      !!
      !! The identifiers are the fields of [[G2008Components]], for example
      !! `G2008Components%methane`. The component order of the model is the
      !! order of `ids`.
      use yaeos__models_ar_multifluid_parameters_gerg2008, only: get_original_parameters
      integer, intent(in) :: ids(:) !! GERG2008 component identifiers
      type(Gerg2008Pure) :: pures(size(ids))
      type(Gerg2008Binary) :: binaries(size(ids), size(ids))

      call get_original_parameters(ids, pures, binaries, gerg_2008%components)
      gerg_2008%pures = pures
      gerg_2008%binaries = binaries
      gerg_2008%ids = ids
      gerg_2008%srk =SoaveRedlichKwong(gerg_2008%components%Tc, gerg_2008%components%Pc, gerg_2008%components%w)
   end function gerg_2008

   subroutine reducing_functions(self, n, Vr, Tr)
      !! # reducing_functions
      !! Evaluate composition-dependent GERG reducing volume and temperature.
      !!
      !! # Description
      !! Normalizes the component amounts to mole fractions and combines pure
      !! critical properties with binary reducing parameters. Hyperdual inputs
      !! preserve composition derivatives for the mixture Helmholtz function.
      class(Gerg2008), intent(in) :: self !! GERG model containing binary reducing data.
      type(hyperdual), intent(in) :: n(:) !! Component mole amounts [mol].
      type(hyperdual), intent(out) :: Vr !! Mixture reducing molar volume [L/mol].
      type(hyperdual), intent(out) :: Tr !! Mixture reducing temperature [K].

      type(hyperdual) :: X(size(n)) !! Mole fractions.

      real(8) :: Vc(size(n))
      real(8) :: Tc(size(n))
      real(8) :: rho_c(size(n))

      real(8) :: Bv(size(n), size(n))
      real(8) :: Gv(size(n), size(n))
      real(8) :: Bt(size(n), size(n))
      real(8) :: Gt(size(n), size(n))

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.
      integer :: nc !! Number of components.

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
      !! # ar_pure
      !! Evaluate one GERG pure-fluid residual Helmholtz contribution.
      !!
      !! # Description
      !! Sums the tabulated polynomial terms and exponentially damped terms in
      !! reduced density and inverse reduced temperature.
      type(Gerg2008Pure), intent(in) :: pure !! Pure-fluid coefficient arrays.
      type(hyperdual), intent(in) :: delta !! Reduced density.
      type(hyperdual), intent(in) :: tau !! Reduced inverse temperature.
      type(hyperdual), intent(out) :: ar !! Pure-fluid residual Helmholtz contribution.

      integer :: i !! Loop/index variable.
      integer :: Kpol
      integer :: Kexp

      real(8) :: n_pol(pure%Kpol)
      real(8) :: d_pol(pure%Kpol)
      real(8) :: t_pol(pure%Kpol)
      real(8) :: n_exp(pure%Kexp)
      real(8) :: d_exp(pure%Kexp)
      real(8) :: t_exp(pure%Kexp)
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
      !! # ar_ij
      !! Evaluate the GERG binary departure contribution for one component pair.
      !!
      !! # Description
      !! Sums the pair's polynomial and exponentially damped departure terms
      !! using its tabulated binary coefficients.
      type(hyperdual), intent(in) :: delta !! Mixture reduced density.
      type(hyperdual), intent(in) :: tau !! Mixture reduced inverse temperature.
      type(Gerg2008Binary), intent(in) :: binary !! Binary departure coefficients for the pair.
      type(hyperdual), intent(out) :: aij !! Pair departure contribution.

      integer :: idx_poly
      integer :: idx_exp

      real(pr) :: n_pol(binary%Kpolij)
      real(pr) :: d_pol(binary%Kpolij)
      real(pr) :: t_pol(binary%Kpolij)
      real(pr) :: n_exp(binary%Kexpij)
      real(pr) :: d_exp(binary%Kexpij)
      real(pr) :: t_exp(binary%Kexpij)
      real(pr) :: etha(binary%Kexpij)
      real(pr) :: eps(binary%Kexpij) !! Density offsets \(\varepsilon\) of the quadratic damping term.
      real(pr) :: beta(binary%Kexpij) !! Coefficients \(\beta\) of the linear damping term.
      real(pr) :: gama(binary%Kexpij)

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
      !! # arfun
      !! Evaluate the GERG2008 mixture residual Helmholtz energy.
      !!
      !! # Description
      !! Forms mixture reducing functions, computes reduced density and
      !! temperature, then sums composition-weighted pure-fluid contributions
      !! and binary departure terms. Hyperdual arithmetic supplies derivatives.
      class(Gerg2008) :: self !! Model with component and binary parameter data.
      type(hyperdual), intent(in) :: n(:) !! Component mole amounts [mol].
      type(hyperdual), intent(in) :: v !! Total volume [L].
      type(hyperdual), intent(in) :: t !! Temperature [K].
      type(hyperdual) :: arval !! Mixture residual Helmholtz energy [bar L].

      type(hyperdual) :: Vr
      type(hyperdual) :: Tr
      type(hyperdual) :: X(size(n)) !! Mole fractions.
      type(hyperdual) :: rho_r !! Mixture reducing density [mol/L].
      type(hyperdual) :: delta !! Reduced density.
      type(hyperdual) :: tau !! Reduced inverse temperature.
      type(hyperdual) :: aij
      type(hyperdual) :: ar_pures(size(n))

      type(Gerg2008Pure) :: pures(size(n))
      type(Gerg2008Binary) :: binary

      real(pr) :: rho_c(size(n))
      real(pr) :: Fij(size(n), size(n))
      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.
      integer :: nc !! Number of components.

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
         do j=1,nc !i+1,nc
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
      !! are: `["liquid", "vapor", "stable"]`. The Newton iterations are done
      !! in \(\ln V\); `"stable"` returns the root with lower residual Gibbs
      !! energy. `V` is set to `-1` if the last Newton solve fails.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = gerg_2008([G2008Components%methane, G2008Components%ethane])
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
      real(pr) :: tol=1e-8 !! Numerical convergence tolerance.

      real(pr) :: totnRT
      real(pr) :: GrL
      real(pr) :: GrV
      real(pr) :: Gr
      real(pr) :: Vliq
      real(pr) :: Vvap
      logical :: failed !! True when the solver did not converge.

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
         !! Pressure residual callback used by the GERG volume Newton solve.
         real(pr), intent(in) :: x !! Logarithm of the trial volume.
         real(pr), intent(out) :: f !! Pressure residual [bar].
         real(pr), intent(out) :: df !! Pressure-residual derivative with respect to log-volume.
         real(pr) :: Ar
         real(pr) :: ArV
         real(pr) :: ArV2
         real(pr) :: Pcalc
         real(pr) :: dPcalcdV
         real(pr) :: Vin
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
      !! # volume_initalizer
      !! Return the SRK covolume used as the initial volume estimate for GERG.
      !!
      !! # Description
      !! GERG has no simple analytic covolume, so its embedded SRK model supplies
      !! the initial volume estimate expected by the `ArModel` interface.
      class(Gerg2008), intent(in) :: self !! GERG model containing its SRK helper.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
      real(pr), intent(in) :: p !! Pressure [bar].
      real(pr), intent(in) :: t !! Temperature [K].
      real(pr) :: v0 !! Initial total volume [L].
      v0 = self%srk%get_v0(n, p, t)
   end function volume_initalizer

   subroutine ln_activity_coefficient(&
      eos, n, P, T, root_type, lngamma, dlngammadP, dlngammadT, dlngammadn &
      )
      !! Calculate natural logarithm of activity coefficients and its
      !! derivatives given pressure and temperature.
      !!
      !! # Description
      !! The activity coefficient is \(\ln\gamma_i = \ln\hat\phi_i - \ln\phi_i^0\),
      !! where the mixture uses `root_type` and each pure fluid \(\phi_i^0\)
      !! uses its stable root at the same \(P\) and \(T\).
      !!
      !! # Examples
      !!
      !! ```fortran
      !! eos = gerg_2008([G2008Components%methane, G2008Components%ethane])
      !!
      !! n = [1.0_pr, 1.0_pr]
      !! T = 300.0_pr
      !! P = 1.0_pr
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
      !! \(\frac{d\ln\gamma}{dP}\)
      real(pr), intent(out), optional :: dlngammadT(size(n))
      !! \(\frac{d\ln\gamma}{dT}\)
      real(pr), intent(out), optional :: dlngammadn(size(n),size(n))
      !! \(\frac{d\ln\gamma}{dn}\)

      ! Mixture properties
      real(pr) :: lnPhi(size(n)) !! Log fugacity coefficients.
      real(pr) :: dlnPhidT(size(n))
      real(pr) :: dlnPhidn(size(n),size(n))
      real(pr) :: dPdV
      real(pr) :: dPdn(size(n))
      real(pr) :: dVdn(size(n))

      ! Pure properties
      real(pr) :: vi(size(n))
      real(pr) :: vi_temp
      real(pr) :: npure(size(n))
      real(pr) :: lnPhi_i(size(n))
      real(pr) :: dlnPhi_i_dT(size(n))
      real(pr) :: lnPhi_i_temp(1)
      real(pr) :: dlnPhi_i_dT_temp(1)
      type(GERG2008) :: eos_temp

      integer :: i !! Loop/index variable.

      logical :: gam
      logical :: dp
      logical :: dt
      logical :: dn
      logical :: present_derivs

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
