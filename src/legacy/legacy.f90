module legacy_ar_models
   !! # Legacy cubic EOS procedures
   !! Compatibility layer for older applications that use global EOS state and procedure pointers.
   !!
   !! # Description
   !! This module stores model configuration, component data, and mixing
   !! parameters in module variables, then dispatches properties through the
   !! `ar_fun` and `vinit` pointers from `ar_interface`. It is retained for
   !! compatibility only: it is not reentrant or thread-safe, and new code
   !! should use the object-oriented model types exported by `yaeos`.
   !!
   !! # Examples
   !! ```fortran
   !! call setup(nc, nmodel, ntdep, ncomb)
   !! call PR78_factory(z, ac, b, Tc, Pc, omega, k)
   !! ```
   use yaeos__constants, only: pr, R
   use ar_interface, only: ar_fun, vinit
   implicit none

   ! Model settings
   integer :: thermo_model !! Legacy integer selector for the cubic EOS alpha/covolume formulation.
   integer :: tdep !! Legacy selector for temperature-dependent binary interactions.
   integer :: mixing_rule !! Legacy integer selector for the mixing rule.
   integer :: nc !! Number of active components.

   ! Mole fractions
   real(pr), allocatable :: z(:) !! Mole fractions vector

   ! ==========================================================================
   !  Cubic EoS Possible parameters
   ! --------------------------------------------------------------------------
   ! Critical constants
   real(pr), allocatable :: tc(:) !! Critical temperature [K]
   real(pr), allocatable :: pc(:) !! Critical pressure [bar]
   real(pr), allocatable :: dc(:) !! Critical density [mol/L]
   real(pr), allocatable :: w(:) !! Acentric factor

   ! Model parameters
   real(pr), allocatable :: ac(:) !! Critical attractive parameter [bar (L/mol)^2]
   real(pr), allocatable :: b(:) !! repulsive parameter [L]
   real(pr), allocatable :: del1(:) !! $$\delta_1$$ parameter
   real(pr), allocatable :: k(:) !! Attractive parameter constant

   ! Classic VdW mixing rules parameters
   real(pr), allocatable :: kij(:, :) !! Attractive BIP
   real(pr), allocatable :: lij(:, :) !! Repulsive BIP
   real(pr), allocatable :: bij(:, :) !! Pair covolume matrix used by quadratic mixing.

   ! T dependant mixing rule parameters
   real(pr), allocatable :: kij0(:, :) !! Low-temperature binary interaction coefficients.
   real(pr), allocatable :: kinf(:, :) !! Asymptotic binary interaction coefficients.
   real(pr), allocatable :: tstar(:, :) !! Reference temperatures for the binary interaction law.
   ! ==========================================================================

contains
   ! ==========================================================================
   !  Initializer routines
   ! --------------------------------------------------------------------------
   subroutine setup(n, nmodel, ntdep, ncomb)
      !! # setup
      !! Initialize the global configuration and allocate legacy model arrays.
      !!
      !! # Description
      !! Sets component count, EOS selector, temperature-dependence selector,
      !! and mixing-rule selector, then reallocates the shared parameter arrays.
      !! This routine resets the global legacy state and is not thread-safe.
      !!
      !! # Examples
      !! ```fortran
      !! call setup(nc, model_code, temperature_dependence_code, mixing_code)
      !! ```
      integer, intent(in) :: n !! Number of components
      integer, intent(in) :: nmodel !! Number of model
      integer, intent(in) :: ntdep !! Kij dependant of temperature
      integer, intent(in) :: ncomb !! Combining rule

      thermo_model = nmodel
      tdep = ntdep
      mixing_rule = ncomb
      nc = n
      
      if (allocated(tc)) deallocate(tc)
      if (allocated(pc)) deallocate(pc)
      if (allocated(dc)) deallocate(dc)
      if (allocated(w)) deallocate(w)
      if (allocated(ac)) deallocate(ac)
      if (allocated(b)) deallocate(b)
      if (allocated(del1)) deallocate(del1)
      if (allocated(k)) deallocate(k)
      if (allocated(kij)) deallocate(kij)
      if (allocated(lij)) deallocate(lij)
      if (allocated(kinf)) deallocate(kinf)
      if (allocated(tstar)) deallocate(tstar)
      if (allocated(bij)) deallocate(bij)

      allocate(tc(n))
      allocate(pc(n))
      allocate(dc(n))
      allocate(w(n))
      allocate(ac(n))
      allocate(b(n))
      allocate(del1(n))
      allocate(k(n))
      allocate(kij(n, n))
      allocate(lij(n, n))
      allocate(kinf(n, n))
      allocate(tstar(n, n))
      allocate(bij(n, n))
   end subroutine setup

   subroutine PR78_factory(moles_in, ac_in, b_in, tc_in, pc_in, w_in, k_in)
        !! # PR78_factory
        !!
        !! # Description
        !! PengRobinson 78 factory
        !!
        !! Takes either the critical parameters or the fitted model parameters
        !! and gets ones in base of the others
        real(pr), intent(in) :: moles_in(nc) !! Component mole amounts or feed composition used to initialize the legacy model.
        real(pr), optional, intent(in) :: ac_in(nc) !! Optional precomputed critical attractive parameters.
        real(pr), optional, intent(in) :: b_in(nc) !! Optional precomputed covolumes.
        real(pr), optional, intent(in) :: tc_in(nc) !! Optional critical temperatures [K].
        real(pr), optional, intent(in) :: pc_in(nc) !! Optional critical pressures [bar].
        real(pr), optional, intent(in) :: w_in(nc) !! Optional acentric factors.
        real(pr), optional, intent(in) :: k_in(nc) !! Optional alpha-function parameters.

        integer :: i !! Loop/index variable.

        logical :: params_spec
        logical :: critical_spec
        real(pr) :: zc(nc)
        real(pr) :: oma(nc)
        real(pr) :: omb(nc)
        real(pr) :: vceos(nc)
        real(pr) :: al
        real(pr) :: be
        real(pr) :: ga(nc)
        real(pr) :: RTc(nc)

        ar_fun => ar_srkpr
        vinit => cubic_v0

        del1 = 1 + sqrt(2.0_pr)
        z = moles_in

        params_spec = (present(ac_in) .and. present(b_in) .and. present(k_in))
        critical_spec = (present(tc_in) .and. present(pc_in) .and. present(w_in))

        if (params_spec) then
            ac = ac_in
            b = b_in
            k = k_in

            call get_Zc_OMa_OMb(del1, zc, oma, omb)
            Tc = OMb * ac / (OMa * R* b)
            RTc = R * Tc
            Pc = OMb * RTc / b
            Vceos = Zc * RTc / Pc
            al = -0.26992
            be = 1.54226
            ga = 0.37464 - k
            w = 0.5 * (-be + sqrt(be**2 - 4 * al * ga)) / al
        else if (critical_spec) then
            tc = tc_in
            pc = pc_in
            w = w_in
            RTc = R*Tc

            call get_Zc_OMa_OMb(del1, Zc, OMa, OMb)

            ac = OMa * RTc**2 / Pc
            b = OMb * RTc / Pc
            Vceos = Zc * RTc / Pc
            ! k (or m) constant to calculate attractive parameter depending on temperature
            do i=1,nc
               if (w(i) <= 0.491) then
                  ! m from PR
                  k(i) = 0.37464 + 1.54226 * w(i) - 0.26992 * w(i)**2
               else
                  ! PR78
                  k(i) = 0.379642 + 1.48503 * w(i) - 0.164423 * w(i)**2 + 0.016666 * w(i)**3
               end if
            end do
        end if
   end subroutine

   subroutine PR76_factory(moles_in, ac_in, b_in, tc_in, pc_in, w_in, k_in)
      !! # PR76_factory
      !!
      !! # Description
      !! PengRobinson 76 factory
      !!
      !! Takes either the critical parameters or the fitted model parameters
      !! and gets ones in base of the others
      real(pr), intent(in) :: moles_in(nc) !! Component mole amounts or feed composition used to initialize the legacy model.
      real(pr), optional, intent(in) :: ac_in(nc) !! Optional precomputed critical attractive parameters.
      real(pr), optional, intent(in) :: b_in(nc) !! Optional precomputed covolumes.
      real(pr), optional, intent(in) :: tc_in(nc) !! Optional critical temperatures [K].
      real(pr), optional, intent(in) :: pc_in(nc) !! Optional critical pressures [bar].
      real(pr), optional, intent(in) :: w_in(nc) !! Optional acentric factors.
      real(pr), optional, intent(in) :: k_in(nc) !! Optional alpha-function parameters.

      integer :: i !! Loop/index variable.

      logical :: params_spec
      logical :: critical_spec
      real(pr) :: zc(nc)
      real(pr) :: oma(nc)
      real(pr) :: omb(nc)
      real(pr) :: vceos(nc)
      real(pr) :: al
      real(pr) :: be
      real(pr) :: ga(nc)
      real(pr) :: RTc(nc)

      ar_fun => ar_srkpr
      vinit => cubic_v0

      del1 = 1 + sqrt(2.0_pr)
      z = moles_in

      params_spec = (present(ac_in) .and. present(b_in) .and. present(k_in))
      critical_spec = (present(tc_in) .and. present(pc_in) .and. present(w_in))

      if (params_spec) then
         ac = ac_in
         b = b_in
         k = k_in

         call get_Zc_OMa_OMb(del1, zc, oma, omb)
         Tc = OMb * ac / (OMa * R* b)
         RTc = R * Tc
         Pc = OMb * RTc / b
         Vceos = Zc * RTc / Pc
         al = -0.26992
         be = 1.54226
         ga = 0.37464 - k
         w = 0.5 * (-be + sqrt(be**2 - 4 * al * ga)) / al
      else if (critical_spec) then
         tc = tc_in
         pc = pc_in
         w = w_in
         RTc = R*Tc

         call get_Zc_OMa_OMb(del1, Zc, OMa, OMb)

         ac = OMa * RTc**2 / Pc
         b = OMb * RTc / Pc
         Vceos = Zc * RTc / Pc
         ! k (or m) constant to calculate attractive parameter depending on temperature
         do i=1,nc
            k(i) = 0.37464 + 1.54226 * w(i) - 0.26992 * w(i)**2
         end do
      end if
      ! ac = 0.45723553_pr * R**2 * tc**2 / pc
      ! b = 0.07779607_pr * R * tc/pc
      ! k = 0.37464_pr + 1.54226_pr * w - 0.26993_pr * w**2
   end subroutine

   subroutine SRK_factory(moles_in, ac_in, b_in, tc_in, pc_in, w_in, k_in)
        !! # SRK_factory
        !!
        !! # Description
        !! SoaveRedlichKwong factory
        !!
        !! Takes either the critical parameters or the fitted model parameters
        !! and gets ones in base of the others
        real(pr), intent(in) :: moles_in(nc) !! Component mole amounts or feed composition used to initialize the legacy model.
        real(pr), optional, intent(in) :: ac_in(nc) !! Optional precomputed critical attractive parameters.
        real(pr), optional, intent(in) :: b_in(nc) !! Optional precomputed covolumes.
        real(pr), optional, intent(in) :: tc_in(nc) !! Optional critical temperatures [K].
        real(pr), optional, intent(in) :: pc_in(nc) !! Optional critical pressures [bar].
        real(pr), optional, intent(in) :: w_in(nc) !! Optional acentric factors.
        real(pr), optional, intent(in) :: k_in(nc) !! Optional alpha-function parameters.

        logical :: params_spec
        logical :: critical_spec
        real(pr) :: zc(nc)
        real(pr) :: oma(nc)
        real(pr) :: omb(nc)
        real(pr) :: vceos(nc)
        real(pr) :: al
        real(pr) :: be
        real(pr) :: ga(nc)
        real(pr) :: RTc(nc)

        integer :: i !! Loop index `i`.
        integer :: j !! Loop index `j`.

        ar_fun => ar_srkpr
        vinit => cubic_v0

        del1 = 1
        z = moles_in

        params_spec = (present(ac_in) .and. present(b_in) .and. present(k_in))
        critical_spec = (present(tc_in) .and. present(pc_in) .and. present(w_in))

        if (params_spec) then
            ac = ac_in
            b = b_in
            k = k_in

            call get_Zc_OMa_OMb(del1, zc, oma, omb)
            Tc = OMb * ac / (OMa * R* b)
            RTc = R * Tc
            Pc = OMb * RTc / b
            Vceos = Zc * RTc / Pc
            dc = 1/vceos
            al = -0.26992
            be = 1.54226
            ga = 0.37464 - k
            w = 0.5 * (-be + sqrt(be**2 - 4 * al * ga)) / al
        else if (critical_spec) then
            tc = tc_in
            pc = pc_in
            w = w_in
            RTc = R * Tc

            call get_Zc_OMa_OMb(del1, Zc, OMa, OMb)

            ac = OMa * RTc**2 / Pc
            b = OMb * RTc / Pc
            Vceos = Zc * RTc / Pc

            k = 0.48 + 1.574 * w - 0.175 * w**2
        end if
   end subroutine

   subroutine get_Zc_OMa_OMb(del1, Zc, OMa, OMb)
      !! # get_Zc_OMa_OMb
      !!
      !! # Description
      !! Calculate Zc, OMa and OMb from the delta_1 parameter.
      real(pr), intent(in)  :: del1(:) !! delta_1 parameter
      real(pr), intent(out) :: Zc(:) !! Critical compressibility factor
      real(pr), intent(out) :: OMa(:) !! OMa
      real(pr), intent(out) :: OMb(:) !! OMb

      real(pr) :: d1(size(del1))
      real(pr) :: y(size(del1)) !! Auxiliary root \(y\) of the cubic critical condition.

      d1 = (1._pr + del1**2._pr)/(1._pr + del1)
      y = 1._pr + (2._pr*(1._pr + del1))**(1.0_pr/3._pr) + (4._pr/(1._pr + del1))**(1.0_pr/3)
      OMa = (3._pr*y*y + 3._pr*y*d1 + d1**2._pr + d1 - 1.0_pr)/(3._pr*y + d1 - 1.0_pr)**2._pr
      OMb = 1._pr/(3._pr*y + d1 - 1.0_pr)
      Zc = y/(3._pr*y + d1 - 1.0_pr)
   end subroutine get_Zc_OMa_OMb
   ! ==========================================================================

   ! ==========================================================================
   !  Ar Functions
   ! --------------------------------------------------------------------------
   subroutine ar_srkpr(z, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
      !! # ar_srkpr
      !!
      !! # Description
      !! Wrapper subroutine to the SRK/PR Residula Helmholtz function to
      !! use the general interface
      real(pr), intent(in) :: z(:) !! Number of moles
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(in) :: t !! Temperature [K]

      real(pr), intent(out) :: ar !! Residual Helmholtz
      real(pr), intent(out) :: arv !! dAr/dV
      real(pr), intent(out) :: artv !! dAr2/dTV
      real(pr), intent(out) :: arv2 !! dAr2/dV2
      real(pr), intent(out) :: Arn(size(z)) !! dAr/dn
      real(pr), intent(out) :: ArVn(size(z)) !! dAr2/dVn
      real(pr), intent(out) :: ArTn(size(z)) !! dAr2/dTn
      real(pr), intent(out) :: Arn2(size(z), size(z)) !! dAr2/dn2

      integer :: nd !! Compositional derivatives
      integer :: nt !! Temperature derivatives

      nd = 2
      nt = 1
      call HelmSRKPR(size(z), nd, nt, z, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
   end subroutine
   
   subroutine ar_rkpr(z, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
      !! # ar_rkpr
      !! Evaluate residual Helmholtz derivatives for the legacy RKPR model.
      !!
      !! # Description
      !! Dispatches the requested properties to the legacy RKPR Helmholtz
      !! expression using the current module-global component and mixing data.
      real(pr), intent(in) :: z(:) !! Number of moles
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(in) :: t !! Temperature [K]

      real(pr), intent(out) :: ar !! Residual Helmholtz
      real(pr), intent(out) :: arv !! dAr/dV
      real(pr), intent(out) :: artv !! dAr2/dTV
      real(pr), intent(out) :: arv2 !! dAr2/dV2
      real(pr), intent(out) :: Arn(size(z)) !! dAr/dn
      real(pr), intent(out) :: ArVn(size(z)) !! dAr2/dVn
      real(pr), intent(out) :: ArTn(size(z)) !! dAr2/dTn
      real(pr), intent(out) :: Arn2(size(z), size(z)) !! dAr2/dn2
      
      integer :: nd !! Compositional derivatives
      integer :: nt !! Temperature derivatives

      nd = 2
      nt = 1
      call HelmRKPR(size(z), nd, nt, z, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
   end subroutine
   
   subroutine HelmSRKPR(nc, ND, NT, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
      !! # HelmSRKPR
      !! Evaluate the SRK/PR-family residual Helmholtz expression and derivatives.
      !!
      !! # Description
      !! Combines residual Helmholtz contributions and requested derivatives
      !! for the legacy cubic-model selectors.
      integer, intent(in) :: nc !! Number of components
      integer, intent(in) :: nd !! Compositional derivatives
      integer, intent(in) :: nt !! Temperature derivatives
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), intent(in) :: rn(nc) !! Number of moles

      real(pr), intent(out) :: ar !! Residual Helmholtz
      real(pr), intent(out) :: arv !! dAr/dV
      real(pr), intent(out) :: artv !! dAr2/dTV
      real(pr), intent(out) :: arv2 !! dAr2/dV2
      real(pr), intent(out) :: Arn(nc) !! dAr/dn
      real(pr), intent(out) :: ArVn(nc) !! dAr2/dVn
      real(pr), intent(out) :: ArTn(nc) !! dAr2/dTn
      real(pr), intent(out) :: Arn2(nc, nc) !! dAr2/dn2
      real(pr) :: ArT
      real(pr) :: ArTT

      real(pr) :: Bmix !! Mixture covolume.
      real(pr) :: dBi(nc) !! Mole-number gradient of the mixture covolume.
      real(pr) :: dBij(nc, nc) !! Mole-number Hessian of the mixture covolume.
      real(pr) :: D !! Attractive mixture parameter.
      real(pr) :: dDi(nc) !! Mole-number gradient of the attractive mixture parameter.
      real(pr) :: dDij(nc, nc) !! Mole-number Hessian of the attractive mixture parameter.
      real(pr) :: dDiT(nc)
      real(pr) :: dDdT !! Temperature derivative of the attractive mixture parameter.
      real(pr) :: dDdT2 !! Second temperature derivative of the attractive mixture parameter.

      real(pr) :: totn
      real(pr) :: d1
      real(pr) :: d2

      real(pr) :: f
      real(pr) :: g
      real(pr) :: fv
      real(pr) :: fB
      real(pr) :: gv
      real(pr) :: fv2
      real(pr) :: gv2
      real(pr) :: AUX
      real(pr) :: FFB
      real(pr) :: FFBV
      real(pr) :: FFBB


      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.

      real(pr) :: b_v
      real(pr) :: a


      TOTN = sum(rn)
      D1 = del1(1)
      D2 = (1._pr - D1)/(1._pr + D1)

      if (mixing_rule .lt. 2) then
         call Bnder(nc, rn, Bmix, dBi, dBij)
         call DandTnder(NT, nc, T, rn, D, dDi, dDiT, dDij, dDdT, dDdT2)
      end if

      ! The f's and g's used here are for Ar, not F (reduced Ar)
      ! This requires to multiply by R all g, f and its derivatives as defined by Mollerup
      f = log((V + D1*Bmix)/(V + D2*Bmix))/Bmix/(D1 - D2)
      g = R*log(1 - Bmix/V)
      fv = -1/((V + D1*Bmix)*(V + D2*Bmix))
      fB = -(f + V*fv)/Bmix
      gv = R*Bmix/(V*(V - Bmix))
      fv2 = (-1/(V + D1*Bmix)**2 + 1/(V + D2*Bmix)**2)/Bmix/(D1 - D2)
      gv2 = R*(1/V**2 - 1/(V - Bmix)**2)

      ! Reduced Helmholtz Energy and derivatives
      Ar = -TOTN*g*T - D*f
      ArV = -TOTN*gv*T - D*fv
      ArV2 = -TOTN*gv2*T - D*fv2

      AUX = R*T/(V - Bmix)
      FFB = TOTN*AUX - D*fB
      FFBV = -TOTN*AUX/(V - Bmix) + D*(2*fv + V*fv2)/Bmix
      FFBB = TOTN*AUX/(V - Bmix) - D*(2*f + 4*V*fv + V**2*fv2)/Bmix**2

      do i = 1, nc
         Arn(i) = -g*T + FFB*dBi(i) - f*dDi(i)
         ArVn(i) = -gv*T + FFBV*dBi(i) - fv*dDi(i)
         if (ND .eq. 2) then
            do j = 1, i
               Arn2(i, j) = AUX*(dBi(i) + dBi(j)) - fB*(dBi(i)*dDi(j) + dBi(j)*dDi(i)) &
                            + FFB*dBij(i, j) + FFBB*dBi(i)*dBi(j) - f*dDij(i, j)
               Arn2(j, i) = Arn2(i, j)
            end do
         end if
      end do

      ! TEMPERATURE DERIVATIVES
      if (NT .eq. 1) then
         ArT = -TOTN*g - dDdT*f
         ArTV = -TOTN*gv - dDdT*fV
         ArTT = -dDdT2*f
         do i = 1, nc
            ArTn(i) = -g + (TOTN*AUX/T - dDdT*fB)*dBi(i) - f*dDiT(i)
         end do
      end if
   end subroutine HelmSRKPR

   subroutine HelmRKPR(nco, NDE, NTD, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
      !! # HelmRKPR
      !!
      !! # Description
      !! Calculate the reduced residual Helmholtz Energy and it's derivatives with the RKPR EOS
      integer, intent(in) :: nco !! Number of components in the mixture.
      integer, intent(in) :: NDE !! Compositional derivative order; `2` also computes `Arn2`.
      integer, intent(in) :: NTD !! Flag selecting temperature-dependent derivative calculations.
      real(pr), intent(in) :: rn(nco) !! Component mole amounts used in mixture rules.
      real(pr), intent(in) :: V !! Total volume [L].
      real(pr), intent(in) :: T !! Temperature [K].

      real(pr), intent(out) :: Ar !! Residual Helmholtz energy.
      real(pr), intent(out) :: ArV !! First volume derivative of residual Helmholtz energy.
      real(pr), intent(out) :: ArTV !! Mixed volume-temperature derivative of residual Helmholtz energy.
      real(pr), intent(out) :: ArV2 !! Second volume derivative of residual Helmholtz energy.
      real(pr), intent(out) :: Arn(nco) !! Mole-number gradient of residual Helmholtz energy.
      real(pr), intent(out) :: ArVn(nco) !! Mixed volume/mole-number derivative of residual Helmholtz energy.
      real(pr), intent(out) :: ArTn(nco) !! Mixed temperature/mole-number derivative of residual Helmholtz energy.
      real(pr), intent(out) :: Arn2(nco, nco) !! Mole-number Hessian of residual Helmholtz energy.

      real(pr) :: totn !! Total number of moles.
      real(pr) :: Bmix !! Mixture covolume.
      real(pr) :: dBi(nco) !! Mole-number gradient of the mixture covolume.
      real(pr) :: dBij(nco, nco) !! Mole-number Hessian of the mixture covolume.
      real(pr) :: dD1i(nco) !! Mole-number gradient of the cubic denominator parameter.
      real(pr) :: dD1ij(nco, nco) !! Mole-number Hessian of the cubic denominator parameter.
      real(pr) :: D !! Attractive mixture parameter.
      real(pr) :: dDi(nco) !! Mole-number gradient of the attractive mixture parameter.
      real(pr) :: dDij(nco, nco) !! Mole-number Hessian of the attractive mixture parameter.
      real(pr) :: dDiT(nco)
      real(pr) :: dDdT !! Temperature derivative of the attractive mixture parameter.
      real(pr) :: dDdT2 !! Second temperature derivative of the attractive mixture parameter.
      real(pr) :: D1
      real(pr) :: D2

      ! Auxiliar functions for Ar
      real(pr) :: f
      real(pr) :: g
      real(pr) :: fv
      real(pr) :: fB
      real(pr) :: gv
      real(pr) :: fv2
      real(pr) :: gv2
      real(pr) :: AUX
      real(pr) :: FFB
      real(pr) :: FFBV
      real(pr) :: FFBB
      ! Extra auxiliar functions for RKPR
      real(pr) :: auxD2
      real(pr) :: fD1
      real(pr) :: fBD1
      real(pr) :: fVD1
      real(pr) :: fD1D1

      real(pr) :: ArT
      real(pr) :: ArTT

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.

      nc = nco
      TOTN = sum(rn)

      call DELTAnder(nc, rn, D1, dD1i, dD1ij)

      D2 = (1 - D1)/(1 + D1)

      if (mixing_rule .lt. 2) then
         call Bnder(nc, rn, Bmix, dBi, dBij)
         call DandTnder(NTD, nc, T, rn, D, dDi, dDiT, dDij, dDdT, dDdT2)
      else
         ! call Bcubicnder(nc,rn,Bmix,dBi,dBij)
         ! call DCubicandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
      end if

      !  The f's and g's used here are for Ar, not F (reduced Ar)
      !  This requires to multiply by R all g, f and its derivatives as defined by Mollerup
      f = log((V + D1*Bmix)/(V + D2*Bmix))/Bmix/(D1 - D2)
      g = R * log(1 - Bmix/V)
      fv = -1/((V + D1*Bmix)*(V + D2*Bmix))
      fB = -(f + V*fv)/Bmix
      gv = R*Bmix/(V*(V - Bmix))
      fv2 = (-1/(V + D1*Bmix)**2 + 1/(V + D2*Bmix)**2)/Bmix/(D1 - D2)
      gv2 = R*(1/V**2 - 1/(V - Bmix)**2)

      ! DERIVATIVES OF f WITH RESPECT TO DELTA1
      auxD2 = (1 + 2/(1 + D1)**2)
      fD1 = (1/(V + D1*Bmix) + 2/(V + D2*Bmix)/(1 + D1)**2) - f*auxD2
      fD1 = fD1/(D1 - D2)
      fBD1 = -(fB*auxD2 + D1/(V + D1*Bmix)**2 + 2*D2/(V + D2*Bmix)**2/(1 + D1)**2)
      fBD1 = fBD1/(D1 - D2)
      fVD1 = -(fV*auxD2 + 1/(V + D1*Bmix)**2 + 2/(V + D2*Bmix)**2/(1 + D1)**2)/(D1 - D2)
      fD1D1 = 4*(f - 1/(V + D2*Bmix))/(1 + D1)**3 + Bmix*(-1/(V + D1*Bmix)**2 &
                                                          + 4/(V + D2*Bmix)**2/(1 + D1)**4) - 2*fD1*(1 + 2/(1 + D1)**2)
      fD1D1 = fD1D1/(D1 - D2)

      ! Reduced Helmholtz Energy and derivatives
      Ar = -TOTN*g*T - D*f
      ArV = -TOTN*gv*T - D*fv
      ArV2 = -TOTN*gv2*T - D*fv2

      AUX = R*T/(V - Bmix)
      FFB = TOTN*AUX - D*fB
      FFBV = -TOTN*AUX/(V - Bmix) + D*(2*fv + V*fv2)/Bmix
      FFBB = TOTN*AUX/(V - Bmix) - D*(2*f + 4*V*fv + V**2*fv2)/Bmix**2

      do i = 1, nc
         Arn(i) = -g*T + FFB*dBi(i) - f*dDi(i) - D*fD1*dD1i(i)
         ArVn(i) = -gv*T + FFBV*dBi(i) - fv*dDi(i) - D*fVD1*dD1i(i)
         if (NDE .eq. 2) then
            do j = 1, i
               Arn2(i, j) = AUX*(dBi(i) + dBi(j)) - fB*(dBi(i)*dDi(j) + dBi(j)*dDi(i)) &
                            + FFB*dBij(i, j) + FFBB*dBi(i)*dBi(j) - f*dDij(i, j)
               Arn2(i, j) = Arn2(i, j) - D*fBD1*(dBi(i)*dD1i(j) + dBi(j)*dD1i(i)) &
                            - fD1*(dDi(i)*dD1i(j) + dDi(j)*dD1i(i)) &
                            - D*fD1*dD1ij(i, j) - D*fD1D1*dD1i(i)*dD1i(j)
               Arn2(j, i) = Arn2(i, j)
            end do
         end if
      end do

      ! TEMPERATURE DERIVATIVES
      if (NTD .eq. 1) then
         ArT = -TOTN*g - dDdT*f
         ArTV = -TOTN*gv - dDdT*fV
         ArTT = -dDdT2*f
         do i = 1, nc
            ArTn(i) = -g + (TOTN*AUX/T - dDdT*fB)*dBi(i) - f*dDiT(i) - dDdT*fD1*dD1i(i)
         end do
      end if
   end subroutine HelmRKPR

   subroutine ArVnder(nc, NDER, NTEMP, z, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
      !! # ArVnder
      !! Dispatch the configured legacy residual Helmholtz implementation.
      !!
      !! # Description
      !! Selects the active residual Helmholtz callback through `ar_fun` and
      !! returns the requested residual energy derivatives. `NDER` and `NTEMP`
      !! control the derivative sets supported by the selected implementation.
      integer, intent(in) :: nc !! Number of active components.
      integer, intent(in) :: nder !! Get compositional derivatives
      integer, intent(in) :: ntemp !! Get temperature derivatives
      real(pr), intent(in) :: z(nc) !! Component mole amounts/composition vector.
      real(pr), intent(in) :: V !! Total volume [L].
      real(pr), intent(in) :: T !! Temperature [K].
      real(pr), intent(out) :: ar !! Residual Helmholtz energy.
      real(pr), intent(out) :: arv !! First volume derivative of residual Helmholtz energy.
      real(pr), intent(out) :: artv !! Mixed volume-temperature derivative of residual Helmholtz energy.
      real(pr), intent(out) :: arv2 !! Second volume derivative of residual Helmholtz energy.
      real(pr), dimension(size(z)), intent(out) :: Arn !! Composition gradient of residual Helmholtz energy.
      real(pr), dimension(size(z)), intent(out) :: ArVn !! Mixed volume-composition derivatives.
      real(pr), dimension(size(z)), intent(out) :: ArTn !! Mixed temperature-composition derivatives.
      real(pr), intent(out) :: Arn2(size(z),size(z)) !! Composition Hessian of residual Helmholtz energy.

      vinit => cubic_v0

      call ar_fun(z, v, t, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
   end subroutine ArVnder
   ! ==========================================================================

   ! ==========================================================================
   !  Attractive parameter routines
   ! --------------------------------------------------------------------------
   subroutine aTder(ac, Tc, k, T, a, dadT, dadT2)
      !! # aTder
      !! Evaluate the configured cubic alpha correction and its temperature derivatives.
      !!
      !! # Description
      !! Uses the global EOS selector to choose the Soave-family or RKPR alpha
      !! correlation, then returns the attractive parameter at `T` and its
      !! first two temperature derivatives.
      real(pr), intent(in) :: ac !! Pure-component critical attractive parameter.
      real(pr), intent(in) :: Tc !! Pure-component critical temperature [K].
      real(pr), intent(in) :: k !! Alpha-function parameter.
      real(pr), intent(in) :: T !! Temperature [K].
      real(pr), intent(out) :: a !! Attractive parameter at the state.
      real(pr), intent(out) :: dadT !! First temperature derivative of the attractive parameter.
      real(pr), intent(out) :: dadT2 !! Second temperature derivative of the attractive parameter.

      real(pr) :: Tr

      Tr = T/Tc

      if (thermo_model .le. 3) then
         a = ac*(1 + k*(1 - sqrt(Tr)))**2
         dadT = ac*k*(k - (k + 1)/sqrt(Tr))/Tc
         dadT2 = ac*k*(k + 1)/(2*Tc**2*Tr**1.5D0)
      else if (thermo_model == 4) then
         a = ac*(3/(2 + Tr))**k
         dadT = -k*a/Tc/(2 + Tr)
         dadT2 = -(k + 1)*dadT/Tc/(2 + Tr)
      end if
   end subroutine aTder

   subroutine aijTder(NTD, nc, T, aij, daijdT, daijdT2)
      !! # aijTder
      !! Build pair attractive parameters and their temperature derivatives.
      !!
      !! # Description
      !! Evaluates each pure attractive parameter, applies the configured
      !! binary interaction law, and symmetrizes the resulting pair matrix.
      integer, intent(in) :: ntd !! Flag selecting temperature-dependent derivative calculations.
      integer, intent(in) :: nc !! Number of active components.
      real(pr), intent(in) :: T !! Temperature [K].
      real(pr), intent(out) :: aij(nc, nc) !! Pair attractive-parameter matrix.
      real(pr), intent(out) :: daijdT(nc, nc) !! Temperature derivative of the pair attractive parameters.
      real(pr), intent(out) :: daijdT2(nc, nc) !! Second temperature derivative of the pair attractive parameters.

      real(pr) :: ai(nc) !! Pure-component attractive parameters.
      real(pr) :: daidT(nc)
      real(pr) :: daidT2(nc)

      real(pr) :: aux(nc, nc) !! Intermediate accumulator.
      real(pr) :: ratK(nc, nc)
      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.

      if (tdep .ge. 1) then
         Kij = 0.0D0
         do i = 1, nc
            Kij(:i - 1, i) = Kinf(:i - 1, i) + Kij0(:i - 1, i)*exp(-T/Tstar(:i - 1, i))
         end do
      end if

      do i = 1, nc
         call aTder(ac(i), Tc(i), k(i), T, ai(i), daidT(i), daidT2(i))
         aij(i, i) = ai(i)
         daijdT(i, i) = daidT(i)
         daijdT2(i, i) = daidT2(i)
         if (i .gt. 1) then
            do j = 1, i - 1
               aij(j, i) = sqrt(ai(i)*ai(j))*(1 - Kij(j, i))
               aij(i, j) = aij(j, i)
               if (NTD .eq. 1) then
                  daijdT(j, i) = (1 - Kij(j, i))*(sqrt(ai(i)/ai(j))*daidT(j) &
                                  + sqrt(ai(j)/ai(i))*daidT(i))/2
                  daijdT2(j, i) = (1 - Kij(j, i))*(daidT(j)*daidT(i)/sqrt(ai(i)*ai(j)) &
                                 + sqrt(ai(i)/ai(j))*(daidT2(j) - daidT(j)**2/(2*ai(j))) &
                                 + sqrt(ai(j)/ai(i))*(daidT2(i) - daidT(i)**2/(2*ai(i))))/2
                  daijdT(i, j) = daijdT(j, i)
                  daijdT2(i, j) = daijdT2(j, i)
               end if
            end do
         end if
      end do
   end subroutine aijTder

   subroutine DandTnder(NTD, nc, T, rn, D, dDi, dDiT, dDij, dDdT, dDdT2)
      !! # DandTnder
      !! Evaluate the legacy mixture attractive parameter and derivatives.
      !!
      !! # Description
      !! Combines the `aij` matrix with component mole amounts and, when
      !! `NTD` enables temperature dependence, evaluates temperature gradients.
      integer, intent(in) :: ntd !! Flag selecting temperature-dependent derivative calculations.
      integer, intent(in) :: nc !! Number of active components.

      real(pr), intent(in) :: T !! Temperature [K].
      real(pr), intent(in) :: rn(nc) !! Component mole amounts used in mixture rules.
      real(pr), intent(out) :: D !! Mixture attractive parameter.
      real(pr), intent(out) :: dDiT(nc) !! \(\partial^2 D/\partial n_i \partial T\).
      real(pr), intent(out) :: dDdT !! \(\partial D/\partial T\).
      real(pr), intent(out) :: dDdT2 !! \(\partial^2 D/\partial T^2\).
      real(pr), intent(out) :: dDi(nc) !! \(\partial D/\partial n_i\).
      real(pr), intent(out) :: dDij(nc, nc) !! \(\partial^2 D/\partial n_i \partial n_j\).

      real(pr) :: aij(nc, nc) !! Pair attractive-parameter matrix.
      real(pr) :: daijdT(nc, nc) !! Temperature derivative of the pair attractive parameters.
      real(pr) :: daijdT2(nc, nc) !! Second temperature derivative of the pair attractive parameters.
      real(pr) :: aux !! Intermediate accumulator.
      real(pr) :: aux2 !! Second intermediate accumulator.

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.

      call aijTder(NTD, nc, T, aij, daijdT, daijdT2)

      D = 0
      dDdT = 0
      dDdT2 = 0
      do i = 1, nc
         aux = 0
         aux2 = 0
         dDi(i) = 0
         dDiT(i) = 0

         do j = 1, nc
            dDi(i) = dDi(i) + 2*rn(j)*aij(i, j)
            if (NTD .eq. 1) then
               dDiT(i) = dDiT(i) + 2*rn(j)*daijdT(i, j)
               aux2 = aux2 + rn(j)*daijdT2(i, j)
            end if

            dDij(i, j) = 2*aij(i, j)
            aux = aux + rn(j)*aij(i, j)
         end do

         D = D + rn(i)*aux
         if (NTD .eq. 1) then
            dDdT = dDdT + rn(i)*dDiT(i)/2
            dDdT2 = dDdT2 + rn(i)*aux2
         end if
      end do
   end subroutine DandTnder
   ! ==========================================================================

   subroutine DELTAnder(nc, rn, D1m, dD1i, dD1ij)
      !! # DELTAnder
      !! Evaluate the mole-averaged cubic denominator parameter and derivatives.
      !!
      !! # Description
      !! Computes the mixture \(\delta_1\) parameter and its first and second
      !! derivatives with respect to component mole amounts.
      integer, intent(in) :: nc !! Number of active components.
      real(pr), intent(in) :: rn(nc) !! Component mole amounts used in mixture rules.
      real(pr), intent(out) :: D1m !! Mixture cubic denominator parameter.
      real(pr), intent(out) :: dD1i(nc) !! Mole-number gradient of the cubic denominator parameter.
      real(pr), intent(out) :: dD1ij(nc, nc) !! Mole-number Hessian of the cubic denominator parameter.

      real(pr) :: totn !! Total number of moles.

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.

      D1m = 0.0_pr
      do i = 1, nc
         D1m = D1m + rn(i) * del1(i)
      end do

      TOTN = sum(rn)
      D1m = D1m/totn

      do i = 1, nc
         dD1i(i) = (del1(i) - D1m)/totn
         do j = 1, nc
            dD1ij(i, j) = (2.0_pr*D1m - del1(i) - del1(j))/totn**2
         end do
      end do
   end subroutine DELTAnder

   ! ==========================================================================
   !  Repulsive parameter routines
   ! --------------------------------------------------------------------------
   subroutine Bnder(nc, rn, Bmix, dBi, dBij)
      !! # Bnder
      !! Evaluate the legacy quadratic-mixture covolume and derivatives.
      !!
      !! # Description
      !! Uses the global binary repulsive interaction data to construct the
      !! mixture covolume and its mole-number gradient and Hessian.
      integer, intent(in) :: nc !! Number of active components.
      real(pr), intent(in) :: rn(nc) !! Component mole amounts used in mixture rules.
      real(pr), intent(out) :: Bmix !! Mixture covolume.
      real(pr), intent(out) :: dBi(nc) !! Mole-number gradient of the mixture covolume.
      real(pr), intent(out) :: dBij(nc, nc) !! Mole-number Hessian of the mixture covolume.

      real(pr) :: totn !! Total number of moles.
      real(pr) :: aux(nc) !! Intermediate accumulator.

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.

      TOTN = sum(rn)
      Bmix = 0.0_pr
      aux = 0.0_pr

      do i = 1, nc
         do j = 1, nc
            bij(i, j) = (b(i) + b(j)) * 0.5_pr * (1.0_pr - lij(i, j))
            aux(i) = aux(i) + rn(j)*bij(i, j)
         end do
         Bmix = Bmix + rn(i)*aux(i)
      end do

      Bmix = Bmix/totn

      do i = 1, nc
         dBi(i) = (2*aux(i) - Bmix)/totn
         do j = 1, i
            dBij(i, j) = (2*bij(i, j) - dBi(i) - dBi(j))/totn
            dBij(j, i) = dBij(i, j)
         end do
      end do
   end subroutine Bnder
   ! ==========================================================================

   ! ==========================================================================
   !  Properties
   ! --------------------------------------------------------------------------
   function cubic_v0(z, p, t)
      !! # cubic_v0
      !! Return a cubic-EOS liquid-volume initial estimate for legacy solvers.
      !!
      !! # Description
      !! Calculates the covolume-based estimate using the globally configured
      !! EOS parameters; it is not a converged volume root.
      real(pr) :: z(:) !! Feed or overall composition.
      real(pr) :: p
      real(pr) :: t
      real(pr) :: cubic_v0

      real(pr) :: dbi(nc)
      real(pr) :: dbij(nc,nc)
      call bnder(nc, z, cubic_v0, dBi, dBij)
   end function
end module
module legacy_thermo_properties
   use yaeos__constants, only: R, pr
   use legacy_ar_models, only: ArVnder, vinit
   implicit none
contains
   
   subroutine TERMO(nc, MTYP, INDIC, T, P, rn, V, PHILOG, DLPHIP, DLPHIT, FUGN)
      !! # TERMO
      !! Evaluate legacy fugacity and pressure-derivative properties.
      !!
      !! # Description
      !! Dispatches the selected cubic model and mixing rule to calculate
      !! fugacity coefficients and requested pressure/temperature derivatives.
      !  MTYP      TYPE OF ROOT DESIRED (-1 vapor, 1 liquid, 0 lower Gibbs energy phase)
      !  rn        mixture mole numbers                        (input)
      !  t         temperature (k)                             (input)x, y
      !  p         pressure    (bar)                          (input)
      !  v         volume      (L)                            (output)
      !  PHILOG    vector of ln(phi(i)*P)                     (output)   INDIC < 5
      !  DLPHIT    t-derivative of ln(phi(i)) (const P, n)    (output)   INDIC = 2 or 4
      !  DLPHIP    P-derivative of ln(phi(i)) (const T, n)    (output)   INDIC < 5
      !  FUGN      comp-derivative of ln(phi(i)) (const t & P)(output)   INDIC > 2
      !  -------------------------------------------------------------------------
      integer, intent(in) :: nc !! Number of components
      integer, intent(in) :: indic !! Desired element, this should be setted with optionals
      integer, intent(in) :: mtyp !! Type of root desired (-1 vapor, 1 liquid, 0 lower Gr)
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), intent(in) :: p !! Pressure [bar]
      real(pr), intent(in) :: rn(nc) !! Mixture mole numbers

      real(pr), intent(out) :: v !! Volume [L]
      real(pr), intent(out) :: PHILOG(nc) !! ln(phi*p) vector
      real(pr), optional, intent(out) :: DLPHIT(nc) !! ln(phi) Temp derivative
      real(pr), optional, intent(out) :: DLPHIP(nc) !! ln(phi) Pressure derivative
      real(pr), optional, intent(out) :: FUGN(nc, nc) !! ln(phi) compositional derivative

      real(pr) :: ar
      real(pr) :: arv
      real(pr) :: artv
      real(pr) :: arv2
      real(pr) :: RT
      real(pr) :: Z
      real(pr) :: dpv
      real(pr) :: dpdt
      real(pr) :: Arn(nc)
      real(pr) :: ArVn(nc)
      real(pr) :: ArTn(nc)
      real(pr) :: Arn2(nc, nc)
      real(pr) :: DPDN(nc)
      real(pr) :: totn !! Total number of moles.
      integer :: ntemp
      integer :: igz
      integer :: nder
      integer :: i !! Loop/index variable.
      integer :: k


      !  The output PHILOG is actually the vector ln(phi(i)*P)
      NTEMP = 0
      IGZ = 0
      NDER = 1
      if (INDIC .gt. 2) NDER = 2
      if (INDIC .eq. 2 .or. INDIC .eq. 4) NTEMP = 1
      TOTN = sum(rn)
      ! if (P .le. 0.0d0) MTYP = 1
      call VCALC(MTYP, NC, NTEMP, rn, T, P, V)
      RT = R*T
      Z = V/(TOTN*RT) ! this is Z/P
      call ArVnder(nc, NDER, NTEMP, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
      DPV = -ArV2 - RT*TOTN/V**2
      DPDT = -ArTV + TOTN*R/V

      do I = 1, NC
         PHILOG(I) = -log(Z) + Arn(I)/RT
         DPDN(I) = RT/V - ArVn(I)
         if (present(dlphip)) DLPHIP(I) = -DPDN(I)/DPV/RT - 1.D0/P
         if (NTEMP .ne. 0) then
            if (present(dlphit)) then
               DLPHIT(I) = (ArTn(I) - Arn(I)/T)/RT + DPDN(I)*DPDT/DPV/RT + 1.D0/T
            end if
         end if
      end do

      if (present(fugn)) then
         do I = 1, NC
            do K = I, NC
               FUGN(I, K) = 1.D0/TOTN + (Arn2(I, K) + DPDN(I)*DPDN(K)/DPV)/RT
               FUGN(K, I) = FUGN(I, K)
            end do
         end do
      end if
   end subroutine TERMO

   subroutine zTVTERMO(nc, INDIC, T, rn, V, P, DPV, PHILOG, DLPHIP, DLPHIT, FUGN)
      !! # zTVTERMO
      !! Evaluate legacy thermodynamic properties at fixed composition, temperature, and volume.
      !!
      !! # Description
      !! Calculates pressure, its volume derivative, fugacity coefficients, and
      !! their pressure/temperature/composition derivatives from global model settings.
      !! Calculation of lnphi*P and derivatives
      !!  rn        mixture mole numbers                       (input)
      !!  t         temperature (k)                            (input)
      !!  v         volume      (L)                            (input)
      !!  p         pressure    (bar)                          (output)
      !!  PHILOG    vector of ln(phi(i)*P)                     (output)  0 < INDIC < 5
      !!  DLPHIT    t-derivative of ln(phi(i)) (const P, n)    (output)  0 < INDIC = 2 or 4
      !!  DLPHIP    P-derivative of ln(phi(i)) (const T, n)    (output)  0 < INDIC < 5
      !!  FUGN      comp-derivative of ln(phi(i)) (const t & P)(output)  2 < INDIC
      !!  -------------------------------------------------------------------------
      implicit none

      integer, intent(in) :: nc !! Number of active components.
      integer, intent(in) :: indic !! Selects the requested outputs (see the `INDIC` table above).
      real(pr), intent(in) :: t !! Temperature [K].
      real(pr), intent(in) :: rn(nc) !! Component mole amounts used in mixture rules.
      real(pr), intent(in) :: v !! Total volume [L].

      real(pr), intent(out) :: p !! Pressure [bar].
      real(pr), intent(out) :: dpv !! Pressure derivative with respect to volume.
      real(pr), intent(out) :: PHILOG(nc) !! \(\ln(\phi_i P)\) vector.
      real(pr), intent(out) :: DLPHIT(nc) !! Temperature derivatives of logarithmic fugacity coefficients.
      real(pr), intent(out) :: DLPHIP(nc) !! Pressure derivatives of logarithmic fugacity coefficients.
      real(pr), intent(out) :: FUGN(nc, nc) !! Composition derivatives of \(\ln \phi_i\) at constant \(T, P\).

      real(pr) :: Arn(nc)
      real(pr) :: ArVn(nc)
      real(pr) :: ArTn(nc)
      real(pr) :: Arn2(nc, nc)
      real(pr) :: DPDN(nc)
      real(pr) :: totn !! Total number of moles.
      real(pr) :: ar
      real(pr) :: arv
      real(pr) :: artv
      real(pr) :: arv2
      real(pr) :: RT
      real(pr) :: Z !! Feed or overall composition.
      real(pr) :: dpdt

      integer :: ntemp
      integer :: igz
      integer :: nder
      integer :: i !! Loop/index variable.
      integer :: k

      NTEMP = 0
      IGZ = 0
      NDER = 1

      if (INDIC .gt. 2) NDER = 2
      if (INDIC .eq. 2 .or. INDIC .eq. 4) NTEMP = 1

      TOTN = sum(rn)

      RT = R*T
      Z = V/(TOTN*RT) ! this is Z/P

      call ArVnder(nc, NDER, NTEMP, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
      P = TOTN*RT/V - ArV
      DPV = -ArV2 - RT*TOTN/V**2
      DPDT = -ArTV + TOTN*R/V

      if (INDIC > 0) then
         do I = 1, NC
            PHILOG(I) = -log(Z) + Arn(I)/RT
            DPDN(I) = RT/V - ArVn(I)
            DLPHIP(I) = -DPDN(I)/DPV/RT - 1.D0/P
            if (NTEMP .ne. 0) then
               DLPHIT(I) = (ArTn(I) - Arn(I)/T)/RT + DPDN(I)*DPDT/DPV/RT + 1.D0/T
            end if
         end do
      end if

      if (NDER .ge. 2) then
         do I = 1, NC
            do K = I, NC
               FUGN(I, K) = 1.D0/TOTN + (Arn2(I, K) + DPDN(I)*DPDN(K)/DPV)/RT
               FUGN(K, I) = FUGN(I, K)
            end do
         end do
      end if
   end subroutine zTVTERMO

   subroutine PUREFUG_CALC(nc, icomp, T, P, V, phi)
      !! # PUREFUG_CALC
      !! Calculate pure-component fugacity data for a selected legacy component.
      !!
      !! # Description
      !! Evaluates the configured EOS for component `icomp` at the supplied
      !! state and returns its fugacity coefficient.
      !! Fugacity of a pure component
      integer, intent(in) :: nc !! Number of active components.
      integer, intent(in) :: icomp !! One-based component index.
      real(pr), intent(in) :: T !! Temperature [K].
      real(pr), intent(in) :: P !! Pressure [bar].
      real(pr), intent(in) :: V !! Total volume [L].
      real(pr), intent(out) :: phi !! Pure-component fugacity coefficient.

      real(pr) :: rn(nc)
      real(pr) :: Ar !! Residual Helmholtz energy.
      real(pr) :: Arv
      real(pr) :: ArTV !! Mixed temperature/volume derivative of residual Helmholtz energy.
      real(pr) :: ArV2 !! Second volume derivative of residual Helmholtz energy.
      real(pr) :: Arn(nc) !! Mole-number gradient of residual Helmholtz energy.
      real(pr) :: ArVn(nc) !! Mixed volume/mole-number derivative of residual Helmholtz energy.
      real(pr) :: ArTn(nc) !! Mixed temperature/mole-number derivative of residual Helmholtz energy.
      real(pr) :: Arn2(nc, nc) !! Mole-number Hessian of residual Helmholtz energy.
      real(pr) :: RT
      real(pr) :: Z !! Feed or overall composition.
      real(pr) :: philog
      rn = 0.0
      rn(icomp) = 1.0
      RT = R*T
      Z = P*V/RT
      call ArVnder(nc, 0, 0, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
      PHILOG = -log(Z) + Arn(icomp)/RT
      phi = exp(PHILOG)
   end subroutine purefug_calc

   recursive subroutine VCALC(ITYP, nc, NTEMP, rn, T, P, V)
      !! # VCALC
      !!
      !! # Description
      !! ROUTINE FOR CALCULATION OF VOLUME, GIVEN PRESSURE
      integer, intent(in) :: ITYP !! TYPE OF ROOT DESIRED (-1 vapor, 1 liquid, 0 lower Gibbs energy phase)
      integer, intent(in) :: nc !! NO. OF COMPONENTS
      integer, intent(in) :: ntemp !! 1 if T-derivatives are required
      real(pr), intent(in) ::  rn(nc) !! FEED MOELS
      real(pr), intent(in) :: T !! TEMPERATURE
      real(pr), intent(in) :: P !! PRESURE
      real(pr), intent(out) :: V !! VOLUME

      real(pr) :: Ar !! Residual Helmholtz energy.
      real(pr) :: ArV !! Volume derivative of residual Helmholtz energy.
      real(pr) :: ArTV !! Mixed temperature/volume derivative of residual Helmholtz energy.
      real(pr) :: ArV2 !! Second volume derivative of residual Helmholtz energy.
      real(pr) :: Arn(nc) !! Mole-number gradient of residual Helmholtz energy.
      real(pr) :: ArVn(nc) !! Mixed volume/mole-number derivative of residual Helmholtz energy.
      real(pr) :: ArTn(nc) !! Mixed temperature/mole-number derivative of residual Helmholtz energy.
      real(pr) :: Arn2(nc, nc) !! Mole-number Hessian of residual Helmholtz energy.
      logical :: FIRST_RUN

      integer :: nder
      real(pr) :: totn !! Total number of moles.
      real(pr) :: B
      real(pr) :: CPV
      real(pr) :: S3R
      real(pr) :: ZETMIN
      real(pr) :: ZETA !! Reduced density \(\zeta = B/V\) used as the volume-solver variable.
      real(pr) :: ZETMAX
      real(pr) :: del
      real(pr) :: pcalc
      real(pr) :: der
      real(pr) :: AT
      real(pr) :: AVAP
      real(pr) :: VVAP

      integer :: iter !! Iteration count.

      NDER = 0
      FIRST_RUN = .true.
      TOTN = sum(rn)
      CPV = vinit(rn, p, t)
      B = CPV
      S3R = 1.D0/CPV
      ITER = 0

      ZETMIN = 0.D0
      !ZETMAX = 1.D0-0.01*T/5000        !.99D0  This is flexible for low T (V very close to B)
      ZETMAX = 1.D0 - 0.01*T/(10000*B) ! improvement for cases with heavy components
      if (ITYP .gt. 0) then
         ZETA = .5D0
      else
         ! IDEAL GAS ESTIMATE
         ZETA = min(.5D0, CPV*P/(TOTN*R*T))
      end if
  100 continue
      DEL = 1
      pcalc = 2*p

      do while(abs(DEL) > 1d-10 .and. iter < 100)
         V = CPV/ZETA
         ITER = ITER + 1
         call ArVnder(&
            nc, NDER, NTEMP, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2 &
         )
         PCALC = TOTN*R*T/V - ArV

         if (PCALC .gt. P) then
            ZETMAX = ZETA
         else
            ZETMIN = ZETA
         end if

         AT = (Ar + V*P)/(T*R) - TOTN*log(V)
         ! AT is something close to Gr(P,T)

         DER = (ArV2*V**2 + TOTN*R*T)*S3R ! this is dPdrho/B
         DEL = -(PCALC - P)/DER
         ZETA = ZETA + max(min(DEL, 0.1D0), -.1D0)

         if (ZETA .gt. ZETMAX .or. ZETA .lt. ZETMIN) &
            ZETA = .5D0*(ZETMAX + ZETMIN)
      end do

      if (ITYP .eq. 0) then
         ! FIRST RUN WAS VAPOUR; RERUN FOR LIQUID
         if (FIRST_RUN) then
            VVAP = V
            AVAP = AT
            FIRST_RUN = .false.
            ZETA = 0.5D0
            ZETMAX = 1.D0 - 0.01*T/500
            goto 100
         else
            if (AT .gt. AVAP) V = VVAP
         end if
      end if
   end subroutine vcalc
   ! ==========================================================================


end module
