module yaeos_c
   !! # Yaeos C Interface
   !! C interface intended to be used on external languanges. With an emphasis
   !! on using it on Python.
   !!
   !! # Description
   !! The interface holds two lists of models (one for `ArModels` and another
   !! for`GeModels`), and two lists of logicals that represent wich models are
   !! in those lists are allocated.
   !!
   !! When a model is instanciated/allocated, it is stored in the singleton
   !! object `x_model` and right after that the `extend_x_models_list` procedure
   !! is called. This procedure searches for the first `free_x_models` id and
   !! allocates the singleton model there, returning the `id` of the procedure.
   use iso_c_binding, only: c_double, c_int, c_int64_t
   use yaeos, only: ArModel, GeModel
   implicit none

   private

   ! CubicEoS
   public :: srk, pr76, pr78, rkpr, psrk
   ! Mixing rules
   public :: set_mhv, set_qmr, set_qmrtd, set_hv

   ! __del__
   public :: make_available_ar_models_list
   public :: make_available_ge_models_list

   ! GeMoels
   public :: nrtl
   public :: unifac_vle
   public :: uniquac
   public :: ln_gamma_ge
   public :: excess_gibbs_ge
   public :: excess_enthalpy_ge
   public :: excess_entropy_ge

   ! Thermoprops
   public :: lnphi_vt, lnphi_pt, pressure, volume, enthalpy_residual_vt
   public :: gibbs_residual_vt, entropy_residual_vt
   public :: Cv_residual_vt, Cp_residual_vt

   ! Phase equilibria
   public :: flash, flash_grid
   public :: saturation_pressure, saturation_temperature
   public :: pt2_phase_envelope
   public :: critical_point, critical_line
   public :: stability_zpt, tm

   type :: ArModelContainer
      !! Container type for ArModels
      class(ArModel), allocatable :: model
   end type ArModelContainer

   type :: GeModelContainer
      !! Container type for GeModels
      class(GeModel), allocatable :: model
   end type GeModelContainer

   class(ArModel), allocatable :: ar_model !! Singleton to hold temporal ArModels
   class(GeModel), allocatable :: ge_model !! Singleton to hold temporal GeModels

   ! Containers of models
   integer, parameter :: max_models  = 1000000
   logical :: free_ar_model(max_models) = .true.
   logical :: free_ge_model(max_models) = .true.

   class(ArModelContainer), allocatable :: ar_models(:)
   class(GeModelContainer), allocatable :: ge_models(:)

contains

   ! ==========================================================================
   !  Ge Models
   ! --------------------------------------------------------------------------
   ! NRTL
   subroutine nrtl(a, b, c, id)
      use yaeos, only: fNRTL => NRTL
      real(c_double), intent(in) :: a(:,:), b(:,:), c(:,:)
      integer(c_int), intent(out) :: id
      ge_model = fNRTL(a, b, c)
      call extend_ge_models_list(id)
   end subroutine nrtl

   ! UNIQUAC
   subroutine uniquac(id, qs, rs, aij, bij, cij, dij, eij)
      use yaeos, only: setup_uniquac
      integer(c_int), intent(out) :: id
      real(c_double), intent(in) :: qs(:)
      !! Molecule's relative areas \(Q_i\)
      real(c_double), intent(in) :: rs(size(qs))
      !! Molecule's relative volumes \(R_i\)
      real(c_double), intent(in) :: aij(size(qs),size(qs))
      !! Interaction parameters matrix \(a_{ij}\)
      real(c_double), intent(in) :: bij(size(qs),size(qs))
      !! Interaction parameters matrix \(b_{ij}\)
      real(c_double), intent(in) :: cij(size(qs),size(qs))
      !! Interaction parameters matrix \(c_{ij}\)
      real(c_double), intent(in) :: dij(size(qs),size(qs))
      !! Interaction parameters matrix \(d_{ij}\)
      real(c_double), intent(in) :: eij(size(qs),size(qs))
      !! Interaction parameters matrix \(e_{ij}\)

      ge_model = setup_uniquac(qs, rs, aij, bij, cij, dij, eij)
      call extend_ge_models_list(id)
   end subroutine uniquac

   ! UNIFAC
   subroutine unifac_vle(id, nc, ngs, g_ids, g_v)
      use yaeos, only: UNIFAC, setup_unifac, Groups
      integer(c_int), intent(out) :: id !! Saved model id
      integer(c_int), intent(in) :: nc !! Number of components
      integer(c_int), intent(in) :: ngs(nc) !! Number of groups at each molecule
      integer(c_int), intent(in) :: g_ids(:, :) !! Ids of groups for each molecule
      integer(c_int), intent(in) :: g_v(:, :) !! Number of groups for each molecule

      type(Groups) :: molecules(nc)

      call setup_groups(nc, ngs, g_ids, g_v, molecules)
      ge_model = setup_unifac(molecules)
      call extend_ge_models_list(id)
   end subroutine unifac_vle

   subroutine setup_groups(nc, ngs, g_ids, g_v, molecules)
      use yaeos, only: Groups
      integer(c_int), intent(in) :: nc !! Number of components
      integer(c_int), intent(in) :: ngs(nc) !! Number of groups at each molecule
      integer(c_int), intent(in) :: g_ids(:, :) !! Ids of groups for each molecule
      integer(c_int), intent(in) :: g_v(:, :) !! Number of groups for each molecule
      type(Groups), intent(out) :: molecules(nc)
      integer :: i

      do i=1,nc
         molecules(i)%groups_ids = g_ids(i, :ngs(i))
         molecules(i)%number_of_groups = g_v(i, :ngs(i))
      end do
   end subroutine setup_groups

   subroutine extend_ge_models_list(id)
      !! Find the first available model container and allocate the model
      !! there. Then return the found id.
      integer(c_int), intent(out) :: id
      integer :: i
      if (.not. allocated(ge_models)) allocate(ge_models(max_models))

      ! Find the first not allocated model
      do i=1,max_models
         if (free_ge_model(i)) then
            free_ge_model(i) = .false.
            id = i
            call move_alloc(ge_model, ge_models(i)%model)
            exit
         end if
      end do
      if (id == max_models) error stop 1
   end subroutine extend_ge_models_list

   subroutine make_available_ge_models_list(id)
      !! Make the geModel id available for allocation
      integer(c_int), intent(in) :: id
      free_ge_model(id) = .true.
   end subroutine make_available_ge_models_list

   ! Ge Thermoprops
   subroutine excess_gibbs_ge(id, n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: n(:)
      !! Moles vector
      real(c_double), intent(in) :: T
      !! Temperature [K]
      real(c_double), intent(out) :: Ge
      !! Excess gibbs energy
      real(c_double), optional, intent(inout) :: GeT
      !! \(\frac{dG^E}{dT}\)
      real(c_double), optional, intent(inout) :: GeT2
      !! \(\frac{d^2G^E}{dT^2}\)
      real(c_double), optional, intent(inout) :: Gen(size(n))
      !! \(\frac{dG^E}{dn_i}\)
      real(c_double), optional, intent(inout) :: GeTn(size(n))
      !! \(\frac{d^2G^E}{dTdn_i}\)
      real(c_double), optional, intent(inout) :: Gen2(size(n), size(n))
      !! \(\frac{d^2G^E}{dn_idn_j}\)

      call ge_models(id)%model%excess_gibbs(&
         n, T, Ge=Ge, GeT=GeT, GeT2=GeT2, Gen=Gen, GeTn=GeTn, Gen2=Gen2 &
         )
   end subroutine excess_gibbs_ge

   subroutine ln_gamma_ge(id, n, T, lngamma, dlngamma_dt, dlngamma_dn)
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: n(:)
      !! Moles vector
      real(c_double), intent(in) :: T
      !! Temperature [K]
      real(c_double), intent(out) :: lngamma(size(n))
      !! Natural logarithm of activity coefficients
      real(c_double), optional, intent(inout) :: dlngamma_dt(size(n))
      !! \(\frac{d\ln \gamma_i}{dT}\)
      real(c_double), optional, intent(inout) :: dlngamma_dn(size(n),size(n))
      !! \(\frac{d\ln \gamma_i}{dn_j}\)

      call ge_models(id)%model%ln_activity_coefficient(&
         n, T, lngamma=lngamma, dlngammadT=dlngamma_dt, dlngammadn=dlngamma_dn&
      )
   end subroutine ln_gamma_ge

   subroutine excess_enthalpy_ge(id, n, T, He, HeT, Hen)
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: n(:)
      !! Moles vector
      real(c_double), intent(in) :: T
      !! Temperature [K]
      real(c_double), intent(out) :: He
      !! Excess enthalpy
      real(c_double), optional, intent(inout) :: HeT
      !! \(\frac{dH^E}{dT}\)
      real(c_double), optional, intent(inout) :: Hen(size(n))
      !! \(\frac{dH^E}{dn}\)

      call ge_models(id)%model%excess_enthalpy(&
         n, T, He=He, HeT=HeT, Hen=Hen &
         )
   end subroutine excess_enthalpy_ge

   subroutine excess_entropy_ge(id, n, T, Se, SeT, Sen)
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: n(:)
      !! Moles vector
      real(c_double), intent(in) :: T
      !! Temperature [K]
      real(c_double), intent(out) :: Se
      !! Excess entropy
      real(c_double), optional, intent(inout) :: SeT
      !! \(\frac{dS^E}{dT}\)
      real(c_double), optional, intent(inout) :: Sen(size(n))
      !! \(\frac{dS^E}{dn}\)

      call ge_models(id)%model%excess_entropy(&
         n, T, Se=Se, SeT=SeT, Sen=Sen &
         )
   end subroutine excess_entropy_ge

   ! ==========================================================================
   !  Ar Models
   ! --------------------------------------------------------------------------
   subroutine extend_ar_models_list(id)
      !! Find the first available model container and allocate the model
      !! there. Then return the found id.
      integer(c_int), intent(out) :: id
      integer :: i
      if (.not. allocated(ar_models)) allocate(ar_models(max_models))

      ! Find the first not allocated model
      do i=1,max_models
         if (free_ar_model(i)) then
            free_ar_model(i) = .false.
            id = i
            call move_alloc(ar_model, ar_models(i)%model)
            exit
         end if
      end do
      if (id == max_models) error stop 1
   end subroutine extend_ar_models_list

   subroutine make_available_ar_models_list(id)
      !! Make the ArModel id available for allocation
      integer(c_int), intent(in) :: id
      free_ar_model(id) = .true.
   end subroutine make_available_ar_models_list

   ! ==========================================================================
   !  Cubic Mixing rules
   ! --------------------------------------------------------------------------
   subroutine set_qmrtd(ar_id, kij_0, kij_inf, t_star, lij)
      use yaeos, only: QMRTD, CubicEoS
      integer(c_int), intent(in) :: ar_id
      real(c_double), intent(in) :: kij_0(:, :)
      real(c_double), intent(in) :: kij_inf(:, :)
      real(c_double), intent(in) :: t_star(:, :)
      real(c_double), intent(in) :: lij(:, :)

      type(QMRTD) :: mixrule

      mixrule = QMRTD(k=kij_inf, k0=kij_0, Tref=t_star, l=lij)

      associate (ar_model => ar_models(ar_id)%model)
         select type(ar_model)
          class is(CubicEoS)
            deallocate(ar_model%mixrule)
            ar_model%mixrule = mixrule
         end select
      end associate
   end subroutine set_qmrtd

   subroutine set_mhv(ar_id, ge_id, q)
      !! Michelsen's Modified Huron-Vidal 1 with constant `q_1` parameter
      use yaeos, only: MHV, CubicEoS
      integer(c_int), intent(in) :: ar_id
      integer(c_int), intent(in) :: ge_id
      real(c_double), intent(in) :: q

      type(MHV) :: mixrule

      ar_model = ar_models(ar_id)%model
      ge_model = ge_models(ge_id)%model

      select type(ar_model)
       class is(CubicEoS)
         mixrule = MHV(ge=ge_model, b=ar_model%b, q=q)
         deallocate(ar_model%mixrule)
         ar_model%mixrule = mixrule
      end select

      call move_alloc(ar_model, ar_models(ar_id)%model)
   end subroutine set_mhv

   subroutine set_hv(ar_id, ge_id)
      !! Huron-Vidal Mixing rule
      use yaeos, only: HV, CubicEoS
      integer(c_int), intent(in) :: ar_id
      integer(c_int), intent(in) :: ge_id

      associate(&
         ar_model => ar_models(ar_id)%model,&
         ge_model => ge_models(ge_id)%model)
         select type(ar_model)
          class is(CubicEoS)
            deallocate(ar_model%mixrule)
            ar_model%mixrule = HV(ge=ge_models(ge_id)%model, bi=ar_model%b, del1=ar_model%del1)
         end select
      end associate
   end subroutine set_hv

   subroutine set_qmr(ar_id, kij, lij)
      use yaeos, only: QMR, CubicEoS
      integer(c_int), intent(in) :: ar_id
      real(c_double) :: kij(:, :)
      real(c_double) :: lij(:, :)

      type(QMR) :: mixrule

      ar_model = ar_models(ar_id)%model

      select type(ar_model)
       class is(CubicEoS)
         mixrule = QMR(k=kij, l=lij)
         deallocate(ar_model%mixrule)
         ar_model%mixrule = mixrule
      end select

      call move_alloc(ar_model, ar_models(ar_id)%model)
   end subroutine set_qmr

   ! ==========================================================================
   !  Cubic EoS implementations
   ! --------------------------------------------------------------------------
   subroutine pr76(tc, pc, w, id)
      use yaeos, only: PengRobinson76
      real(c_double), intent(in) :: tc(:), pc(:), w(:)
      integer(c_int), intent(out) :: id

      ar_model = PengRobinson76(tc, pc, w)
      call extend_ar_models_list(id)
   end subroutine pr76

   subroutine pr78(tc, pc, w, id)
      use yaeos, only: PengRobinson78
      real(c_double), intent(in) :: tc(:), pc(:), w(:)
      integer(c_int), intent(out) :: id

      ar_model = PengRobinson78(tc, pc, w)
      call extend_ar_models_list(id)
   end subroutine pr78

   subroutine srk(tc, pc, w, id)
      use yaeos, only: SoaveRedlichKwong
      real(c_double), intent(in) :: tc(:), pc(:), w(:)
      integer(c_int), intent(out) :: id
      ar_model = SoaveRedlichKwong(tc, pc, w)
      call extend_ar_models_list(id)
   end subroutine srk

   subroutine rkpr(tc, pc, w, zc, delta_1, k, id)
      use yaeos, only: fRKPR => RKPR
      real(c_double), intent(in) :: tc(:), pc(:), w(:), zc(:)
      real(c_double), optional, intent(in) :: delta_1(size(tc)), k(size(tc))
      integer(c_int), intent(out) :: id

      if (all(delta_1 == 0) .and. all(k == 0)) then
         ar_model = fRKPR(tc, pc, w, zc)
      else if (all(delta_1 == 0)) then
         ar_model = fRKPR(tc, pc, w, zc, k=k)
      else if (all(k == 0)) then
         ar_model = fRKPR(tc, pc, w, zc, delta_1=delta_1)
      else
         ar_model = fRKPR(tc, pc, w, zc, delta_1=delta_1, k=k)
      end if
      call extend_ar_models_list(id)
   end subroutine rkpr

   subroutine psrk(id, nc, tc, pc, w, c1, c2, c3, ngs, g_ids, g_v)
      use yaeos, only: Groups, fPSRK => PSRK
      integer(c_int), intent(out) :: id
      real(c_double), intent(in) :: tc(:), pc(:), w(:)
      real(c_double), intent(in) :: c1(:), c2(:), c3(:)
      integer, intent(in) :: nc
      integer, intent(in) :: ngs(nc)
      integer, intent(in) :: g_ids(:, :)
      integer, intent(in) :: g_v(:, :)

      type(Groups) :: molecules(nc)
      call setup_groups(nc, ngs, g_ids, g_v, molecules)
      ar_model = fPSRK(tc, pc, w, molecules, c1, c2, c3)
      call extend_ar_models_list(id)
   end subroutine psrk

   ! ==========================================================================
   !  Thermodynamic properties
   ! --------------------------------------------------------------------------
   subroutine residual_helmholtz(id, n, v, t, ar, ArT, ArV, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2)
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: n(:), v, t
      real(c_double), intent(out) :: ar
      real(c_double), optional, intent(out) :: &
         ArT, ArV, ArTV, ArV2, ArT2, Arn(size(n)), ArVn(size(n)), ArTn(size(n)), Arn2(size(n), size(n))

      call ar_models(id)%model%residual_helmholtz(&
         n=n, V=V, T=T, &
         Ar=Ar,  ArV=ArV, ArT=ArT, ArTV=ArTV, &
         ArV2=ARV2, ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2)
   end subroutine residual_helmholtz

   subroutine lnphi_vt(id, n, v, t, lnphi, dlnphidp, dlnphidt, dlnphidn)
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: n(:), v, t
      real(c_double), intent(out) :: lnphi(size(n))
      real(c_double) :: p

      real(c_double), optional, intent(in out) :: &
         dlnphidp(size(n)), dlnphidt(size(n)), dlnphidn(size(n), size(n))

      call ar_models(id)%model%lnphi_vt(&
         n, V, T, P, lnphi, dlnPhidP, dlnphidT, dlnPhidn &
         )
   end subroutine lnphi_vt

   subroutine lnphi_pt(id, n, p, t, root_type, lnphi, dlnphidp, dlnphidt, dlnphidn)
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: n(:), p, t
      character(len=15), intent(in) :: root_type
      real(c_double), intent(out) :: lnphi(size(n))

      real(c_double), optional, intent(in out) :: &
         dlnphidp(size(n)), dlnphidt(size(n)), dlnphidn(size(n), size(n))

      call ar_models(id)%model%lnphi_pt(&
         n, P=P, T=T, root_type=root_type, &
         lnphi=lnphi, dlnphidp=dlnPhidP, dlnphidt=dlnphidT, dlnphidn=dlnPhidn &
         )
   end subroutine lnphi_pt

   subroutine pressure(id, n, V, T, P, dPdV, dPdT, dPdn)
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: n(:), V, T
      real(c_double), intent(out) :: P
      real(c_double), optional, intent(in out) :: dPdV, dPdT, dPdn(size(n))

      call ar_models(id)%model%pressure(&
         n, V, T, P, dPdV, dPdT, dPdn &
         )
   end subroutine pressure

   subroutine volume(id, n, P, T, root_type, V)
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: n(:), P, T
      character(len=15), intent(in) :: root_type
      real(c_double), intent(out) :: V

      call ar_models(id)%model%volume(n=n, P=P, T=T, root_type=root_type, V=V)
   end subroutine volume

   subroutine enthalpy_residual_vt(id, n, V, T, Hr, HrT, HrV, Hrn)
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: n(:), V, T
      real(c_double), intent(out) :: Hr
      real(c_double), optional, intent(in out) :: HrT, HrV, Hrn(size(n))

      call ar_models(id)%model%enthalpy_residual_vt(&
         n, V, T, Hr, HrT, HrV, Hrn &
         )
   end subroutine enthalpy_residual_vt

   subroutine gibbs_residual_vt(id, n, V, T, Gr, GrT, GrV, Grn)
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: n(:), V, T
      real(c_double), intent(out) :: Gr
      real(c_double), optional, intent(in out) :: GrT, GrV, Grn(size(n))

      call ar_models(id)%model%gibbs_residual_vt(&
         n, V, T, Gr, GrT, GrV, Grn &
         )
   end subroutine gibbs_residual_vt

   subroutine entropy_residual_vt(id, n, V, T, Sr, SrT, SrV, Srn)
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: n(:), V, T
      real(c_double), intent(out) :: Sr
      real(c_double), optional, intent(in out) :: SrT, SrV, Srn(size(n))

      call ar_models(id)%model%entropy_residual_vt(&
         n, V, T, Sr, SrT, SrV, Srn &
         )
   end subroutine entropy_residual_vt

   subroutine Cv_residual_vt(id, n, V, T, Cv)
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: n(:), V, T
      real(c_double), intent(out) :: Cv

      call ar_models(id)%model%cv_residual_vt(&
         n, V, T, Cv &
         )
   end subroutine Cv_residual_vt

   subroutine Cp_residual_vt(id, n, V, T, Cp)
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: n(:), V, T
      real(c_double), intent(out) :: CP

      call ar_models(id)%model%cp_residual_vt(n, V, T, Cp)
   end subroutine Cp_residual_vt

   ! ==========================================================================
   ! Phase equilibria
   ! --------------------------------------------------------------------------
   subroutine critical_point(id, z0, zi, spec, max_iters, x, T, P, V)
      use yaeos, only: EquilibriumState, fcritical_point => critical_point
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z0(:)
      real(c_double), intent(in) :: zi(:)
      integer, intent(in) :: spec
      integer, intent(in) :: max_iters
      real(c_double), intent(out) :: x(size(z0))
      real(c_double), intent(out) :: T
      real(c_double), intent(out) :: P
      real(c_double), intent(out) :: V

      real(c_double) :: y(size(z0)), Vx, Vy, beta
      real(c_double) :: S

      type(EquilibriumState) :: crit

      S = 0
      crit = fcritical_point(&
         model=ar_models(id)%model, z0=z0, zi=zi, &
         S=S, spec=spec, max_iters=max_iters &
         )
      call equilibria_state_to_arrays(crit, x, y, P, T, V, Vy, beta)
   end subroutine critical_point

   subroutine critical_line(&
      id, a0, da0, &
      z0, zi, max_points, stop_pressure, &
      as, Vs, Ts, Ps)
      use yaeos, only: EquilibriumState, CriticalLine, &
         fcritical_line => critical_line, spec_CP
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z0(:)
      real(c_double), intent(in) :: zi(:)
      real(c_double), intent(in) :: a0
      real(c_double), intent(in) :: da0
      integer, intent(in) :: max_points
      real(c_double), intent(in) :: stop_pressure
      real(c_double), intent(out) :: as(max_points)
      real(c_double), intent(out) :: Ts(max_points)
      real(c_double), intent(out) :: Ps(max_points)
      real(c_double), intent(out) :: Vs(max_points)

      type(CriticalLine) :: cl

      integer :: i

      as = makenan()
      Ts = makenan()
      Ps = makenan()
      Vs = makenan()

      cl = fcritical_line(&
         model=ar_models(id)%model, a0=a0, &
         z0=z0, zi=zi, &
         ns=spec_CP%a, S=a0, ds0=da0, maxp=stop_pressure, max_points=max_points)

      do i=1,size(cl%a)
         as(i) = cl%a(i)
         Ts(i) = cl%T(i)
         Ps(i) = cl%P(i)
         Vs(i) = cl%V(i)
      end do
   end subroutine critical_line

   subroutine equilibria_state_to_arrays(eq_state, x, y, P, T, Vx, Vy, beta)
      use yaeos, only: EquilibriumState
      type(EquilibriumState) :: eq_state
      real(c_double), intent(out) :: x(:)
      real(c_double), intent(out) :: y(:)
      real(c_double), intent(out) :: P
      real(c_double), intent(out) :: T
      real(c_double), intent(out) :: Vx
      real(c_double), intent(out) :: Vy
      real(c_double), intent(out) :: Beta

      x = eq_state%x
      y = eq_state%y
      P = eq_state%p
      T = eq_state%T
      Vx = eq_state%Vx
      Vy = eq_state%Vy
      beta = eq_state%beta
   end subroutine equilibria_state_to_arrays

   subroutine flash(id, z, T, P, x, y, k0, Pout, Tout, Vx, Vy, beta)
      use yaeos, only: EquilibriumState, fflash => flash
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:)
      real(c_double), intent(in) :: T
      real(c_double), intent(in) :: P
      real(c_double), intent(in) :: k0(size(z))
      real(c_double), intent(out) :: x(size(z))
      real(c_double), intent(out) :: y(size(z))
      real(c_double), intent(out) :: Pout
      real(c_double), intent(out) :: Tout
      real(c_double), intent(out) :: Vx
      real(c_double), intent(out) :: Vy
      real(c_double), intent(out) :: beta

      type(EquilibriumState) :: result
      integer :: iters

      if (all(k0 == 0)) then
         result = fflash(ar_models(id)%model, z, t, p_spec=p, iters=iters)
      else
         result = fflash(ar_models(id)%model, z, t, p_spec=p, k0=k0, iters=iters)
      end if

      if (.not. allocated(result%x) .or. .not. allocated(result%y)) then
         Pout = P
         Tout = T
         x = z
         y = z
         beta = -1
         Vx = 1
         Vy = 1
         return
      end if

      call equilibria_state_to_arrays(result, x, y, Pout, Tout, Vx, Vy, beta)
   end subroutine flash

   subroutine saturation_pressure(id, z, T, kind, P0, P, x, y, Vx, Vy, beta)
      use yaeos, only: EquilibriumState, fsaturation_pressure => saturation_pressure
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:)
      real(c_double), intent(in) :: T
      character(len=15), intent(in) :: kind
      real(c_double), intent(in) :: P0

      real(c_double), intent(out) :: P
      real(c_double), intent(out) :: x(size(z))
      real(c_double), intent(out) :: y(size(z))
      real(c_double), intent(out) :: Vx, Vy, beta

      real(c_double) :: aux

      type(EquilibriumState) :: sat

      if (P0 == 0) then
         sat = fsaturation_pressure(ar_models(id)%model, z, T, kind)
      else
         sat = fsaturation_pressure(ar_models(id)%model, z, T, kind, P0=P0)
      end if
      call equilibria_state_to_arrays(sat, x, y, P, aux, Vx, Vy, beta)
   end subroutine saturation_pressure

   subroutine saturation_temperature(id, z, P, kind, T0, T, x, y, Vx, Vy, beta)
      use yaeos, only: EquilibriumState, fsaturation_temperature => saturation_temperature
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:)
      real(c_double), intent(in) :: P
      character(len=15), intent(in) :: kind
      real(c_double), intent(in) :: T0

      real(c_double), intent(out) :: T
      real(c_double), intent(out) :: x(size(z))
      real(c_double), intent(out) :: y(size(z))
      real(c_double), intent(out) :: Vx, Vy, beta

      real(c_double) :: aux

      type(EquilibriumState) :: sat

      if (T0 == 0) then
         sat = fsaturation_temperature(ar_models(id)%model, z, P, kind)
      else
         sat = fsaturation_temperature(ar_models(id)%model, z, P, kind, T0=T0)
      end if
      call equilibria_state_to_arrays(sat, x, y, aux, T, Vx, Vy, beta)
   end subroutine saturation_temperature

   subroutine pt2_phase_envelope(id, z, kind, max_points, Ts, Ps, tcs, pcs, T0, P0)
      use yaeos, only: &
         saturation_pressure, saturation_temperature, pt_envelope_2ph, &
         EquilibriumState, PTEnvel2, find_hpl
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:)
      integer, intent(in) :: max_points
      character(len=15), intent(in) :: kind
      real(c_double), intent(out) :: Ts(max_points)
      real(c_double), intent(out) :: Ps(max_points)
      real(c_double), intent(out) :: Tcs(5), Pcs(5)
      real(c_double), optional, intent(in) :: T0, P0

      real(8) :: nan
      type(EquilibriumState) :: sat
      type(PTEnvel2) :: env

      integer :: i, neval=0

      real(c_double) :: T, P

      neval = neval + 1
      nan = makenan()
      Ts = nan
      Ps = nan
      Tcs = nan
      Pcs = nan

      if (present(T0)) then
         T = T0
      else
         T = 150
      end if

      if (present(P0)) then
         P = P0
      else
         P = 1
      end if

      select case(kind)
       case("bubble")
         sat = saturation_pressure(ar_models(id)%model, z, T=T, kind=kind)
         env = pt_envelope_2ph(ar_models(id)%model, z, sat, points=max_points)
       case("dew")
         sat = saturation_temperature(ar_models(id)%model, z, P=P, kind=kind)
         env = pt_envelope_2ph(ar_models(id)%model, z, sat, points=max_points)
       case("liquid-liquid")
         ! sat = saturation_temperature(ar_models(id)%model, z, P=P, kind=kind)
         env = find_hpl(ar_models(id)%model, z, T, P)
      end select


      i = size(env%points)
      Ts(:i) = env%points%T
      Ps(:i) = env%points%P

      i = size(env%cps)
      Tcs(:i) = env%cps%T
      Pcs(:i) = env%cps%P
   end subroutine pt2_phase_envelope

   subroutine flash_grid(id, z, Ts, Ps, xs, ys, Vxs, Vys, betas, parallel)
      use yaeos, only: EquilibriumState, flash
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:)
      real(c_double), intent(in) :: Ts(:)
      real(c_double), intent(in) :: Ps(:)
      real(c_double), dimension(size(Ps), size(Ts), size(z)), intent(out) :: xs, ys
      real(c_double), dimension(size(Ps), size(Ts)), intent(out) :: Vxs, Vys, betas
      logical, intent(in) :: parallel

      class(ArModel), allocatable :: model
      type(EquilibriumState) :: flash_result

      real(8) :: T, P

      integer :: i, j, nt, np, iter

      model = ar_models(id)%model
      np = size(Ps)
      nt = size(Ts)

      if (parallel) then
         !$OMP PARALLEL DO PRIVATE(i, j, t, p, flash_result) SHARED(model, z, ts, ps, betas, Vxs, Vys, xs, ys)
         do i=1,np
            do j=1,nt
               T = Ts(j)
               P = Ps(i)
               flash_result = flash(model, z, T=T, P_spec=P, iters=iter)
               betas(i, j) = flash_result%beta

               Vxs(i, j) = flash_result%Vx
               Vys(i, j) = flash_result%Vy
               xs(i, j, :) = flash_result%x
               ys(i, j, :) = flash_result%y
            end do
         end do
         !$OMP END PARALLEL DO
      else
         do i=1,np
            do j=1,nt
               T = Ts(j)
               P = Ps(i)
               flash_result = flash(model, z, T=T, P_spec=P, iters=iter)
               betas(i, j) = flash_result%beta
               print *, i, j, flash_result%iters, flash_result%beta

               Vxs(i, j) = flash_result%Vx
               Vys(i, j) = flash_result%Vy
               xs(i, j, :) = flash_result%x
               ys(i, j, :) = flash_result%y
            end do
         end do
      end if
   end subroutine flash_grid

   subroutine stability_zpt(id, z, P, T, w_min, tm_val, all_mins)
      use yaeos, only: min_tpd
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:), P, T
      real(c_double), intent(out) :: w_min(size(z))
      real(c_double), intent(out) :: tm_val
      real(c_double), intent(out) :: all_mins(size(z), size(z))

      call min_tpd(&
         ar_models(id)%model, z=z, P=P, T=T, &
         mintpd=tm_val, w=w_min, all_minima=all_mins)
   end subroutine stability_zpt

   subroutine tm(id, z, w, P, T, tm_value)
      use yaeos, only: ftm => tm
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:), w(size(z)), P, T
      real(c_double), intent(out) :: tm_value

      tm_value = ftm(model=ar_models(id)%model, z=z, w=w, P=P, T=T)
   end subroutine tm

   ! ==========================================================================
   ! Auxiliar
   ! --------------------------------------------------------------------------
   function makenan()
      real(c_double) :: makenan
      makenan = 0
      makenan = makenan/makenan
   end function makenan
end module yaeos_c
