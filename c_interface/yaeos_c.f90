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
   public :: srk, pr76, pr78, rkpr, psrk, get_ac_b_del1_del2
   ! Mixing rules
   public :: set_mhv, set_qmr, set_qmrtd, set_hv, set_hvnrtl
   ! Multifluid equations
   public :: multifluid_gerg2008

   ! __del__
   public :: make_available_ar_models_list
   public :: make_available_ge_models_list

   ! GeModels
   public :: nrtl
   public :: unifac_vle
   public :: unifac_psrk
   public :: uniquac
   public :: unifac_dortmund
   public :: ln_gamma_ge
   public :: excess_gibbs_ge
   public :: excess_enthalpy_ge
   public :: excess_entropy_ge

   ! Thermoprops
   public :: lnphi_vt, lnphi_pt, pressure, volume, enthalpy_residual_vt
   public :: gibbs_residual_vt, entropy_residual_vt
   public :: Cv_residual_vt, Cp_residual_vt

   ! Phase equilibria
   public :: flash, flash_vt, flash_grid, solve_mp_flash
   public :: flash_ge
   public :: saturation_pressure, saturation_temperature
   public :: pure_saturation_line
   public :: pt_mp_phase_envelope, px_mp_phase_envelope, tx_mp_phase_envelope
   public :: generalized_isopleth
   public :: critical_point, critical_line, find_llcl
   public :: stability_zpt, tm
   public :: stability_zt_ge

   ! Helpers
   public :: find_self_intersections

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

   ! Dortmund
   subroutine unifac_dortmund(id, nc, ngs, g_ids, g_v)
      use yaeos, only: UNIFAC, setup_dortmund, Groups
      integer(c_int), intent(out) :: id !! Saved model id
      integer(c_int), intent(in) :: nc !! Number of components
      integer(c_int), intent(in) :: ngs(nc) !! Number of groups at each molecule
      integer(c_int), intent(in) :: g_ids(:, :) !! Ids of groups for each molecule
      integer(c_int), intent(in) :: g_v(:, :) !! Number of groups for each molecule

      type(Groups) :: molecules(nc)

      call setup_groups(nc, ngs, g_ids, g_v, molecules)
      ge_model = setup_dortmund(molecules)
      call extend_ge_models_list(id)
   end subroutine unifac_dortmund

   ! PSRK
   subroutine unifac_psrk(id, nc, ngs, g_ids, g_v)
      use yaeos, only: UNIFAC, setup_psrk, Groups
      integer(c_int), intent(out) :: id !! Saved model id
      integer(c_int), intent(in) :: nc !! Number of components
      integer(c_int), intent(in) :: ngs(nc) !! Number of groups at each molecule
      integer(c_int), intent(in) :: g_ids(:, :) !! Ids of groups for each molecule
      integer(c_int), intent(in) :: g_v(:, :) !! Number of groups for each molecule

      type(Groups) :: molecules(nc)

      call setup_groups(nc, ngs, g_ids, g_v, molecules)
      ge_model = setup_psrk(molecules)
      call extend_ge_models_list(id)
   end subroutine unifac_psrk

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

   subroutine set_mhv(ar_id, ge_id, q, lij)
      !! Michelsen's Modified Huron-Vidal 1 with constant `q_1` parameter
      use yaeos, only: MHV, CubicEoS
      integer(c_int), intent(in) :: ar_id
      integer(c_int), intent(in) :: ge_id
      real(c_double), intent(in) :: q
      real(c_double), intent(in) :: lij(:, :)

      type(MHV) :: mixrule

      ar_model = ar_models(ar_id)%model
      ge_model = ge_models(ge_id)%model

      select type(ar_model)
       class is(CubicEoS)
         mixrule = MHV(ge=ge_model, b=ar_model%b, q=q, lij=lij)
         deallocate(ar_model%mixrule)
         ar_model%mixrule = mixrule
      end select

      call move_alloc(ar_model, ar_models(ar_id)%model)
   end subroutine set_mhv

   subroutine set_hvnrtl(ar_id, alpha, gji, use_kij, kij)
      !! Huron-Vidal NRTL mixing rule
      use yaeos, only: fHV_NRTL => HV_NRTL, init_hvnrtl, CubicEoS
      integer(c_int), intent(in) :: ar_id
      real(c_double), intent(in) :: alpha(:, :)
      real(c_double), intent(in) :: gji(:, :)
      logical, intent(in) :: use_kij(:, :)
      real(c_double), intent(in) :: kij(:, :)

      type(fHV_NRTL) :: mixrule

      ar_model = ar_models(ar_id)%model

      select type(ar_model)
       class is(CubicEoS)
         mixrule = init_hvnrtl(b=ar_model%b, &
            del1=ar_model%del1, alpha=alpha, gji=gji, use_kij=use_kij, &
            kij=kij)
         deallocate(ar_model%mixrule)
         ar_model%mixrule = mixrule
      end select

      call move_alloc(ar_model, ar_models(ar_id)%model)
   end subroutine set_hvnrtl

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

   subroutine get_ac_b_del1_del2(id, ac, b, del1, del2, nc)
      use yaeos, only: CubicEoS, size
      integer(c_int), intent(in) :: id
      integer, intent(in) :: nc
      real(c_double), dimension(nc), intent(out) :: &
         ac, b, del1, del2


      associate(model => ar_models(id)%model)
         select type(model)
          class is(CubicEoS)
            ac(:nc) = model%ac
            b(:nc) = model%b
            del1(:nc) = model%del1
            del2(:nc) = model%del2
         end select
      end associate
   end subroutine get_ac_b_del1_del2

   ! ==========================================================================
   !  Multifluid equations
   ! --------------------------------------------------------------------------
   subroutine multifluid_gerg2008(ids, id)
      use yaeos, only: gerg_2008, GERG2008
      integer, intent(in) :: ids(:)
      integer(c_int), intent(out) :: id

      integer :: i
      ar_model = gerg_2008(ids)
      call extend_ar_models_list(id)
   end subroutine multifluid_gerg2008


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
   subroutine critical_point(id, z0, zi, spec, S, max_iters, x, T, P, V)
      use yaeos, only: EquilibriumState, fcritical_point => critical_point
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z0(:)
      real(c_double), intent(in) :: zi(:)
      integer, intent(in) :: spec
      real(c_double), intent(in) :: S
      integer, intent(in) :: max_iters
      real(c_double), intent(out) :: x(size(z0))
      real(c_double), intent(out) :: T
      real(c_double), intent(out) :: P
      real(c_double), intent(out) :: V

      real(c_double) :: y(size(z0)), Vx, Vy, beta

      type(EquilibriumState) :: crit

      crit = fcritical_point(&
         model=ar_models(id)%model, z0=z0, zi=zi, &
         S=S, spec=spec, max_iters=max_iters &
         )
      call equilibria_state_to_arrays(crit, x, y, P, T, V, Vy, beta)
   end subroutine critical_point

   subroutine critical_line(&
      id, ns, S, ds0, &
      a0, v0, t0, p0, z0, zi, stability_analysis, max_points, stop_pressure, &
      as, Vs, Ts, Ps, CEP_x, CEP_y, CEP_P, CEP_Vx, CEP_Vy, CEP_T)
      use yaeos, only: EquilibriumState, CriticalLine, &
         fcritical_line => critical_line, spec_CP
      integer(c_int), intent(in) :: id
      integer(c_int), intent(in) :: ns
      real(c_double), intent(in) :: S
      real(c_double), intent(in) :: dS0
      real(c_double), intent(in) :: a0
      real(c_double), intent(in) :: v0
      real(c_double), intent(in) :: t0
      real(c_double), intent(in) :: p0
      real(c_double), intent(in) :: z0(:)
      real(c_double), intent(in) :: zi(:)
      logical, intent(in) :: stability_analysis
      integer, intent(in) :: max_points
      real(c_double), intent(in) :: stop_pressure
      real(c_double), intent(out) :: as(max_points)
      real(c_double), intent(out) :: Ts(max_points)
      real(c_double), intent(out) :: Ps(max_points)
      real(c_double), intent(out) :: Vs(max_points)
      real(c_double), intent(out) :: CEP_x(size(z0))
      real(c_double), intent(out) :: CEP_y(size(z0))
      real(c_double), intent(out) :: CEP_P
      real(c_double), intent(out) :: CEP_Vx
      real(c_double), intent(out) :: CEP_Vy
      real(c_double), intent(out) :: CEP_T

      type(CriticalLine) :: cl

      integer :: i

      as = makenan()
      Ts = makenan()
      Ps = makenan()
      Vs = makenan()

      if (v0 == 0 .and. t0 == 0 .and. p0 == 0) then
         cl = fcritical_line(&
            model=ar_models(id)%model, a0=a0, &
            z0=z0, zi=zi, &
            ns0=ns, S0=S, ds0=ds0, &
            stability_analysis=stability_analysis, maxp=stop_pressure, max_points=max_points)
      elseif (t0 == 0 .and. p0 == 0) then
         cl = fcritical_line(&
            model=ar_models(id)%model, a0=a0, &
            z0=z0, zi=zi, &
            ns0=ns, S0=S, ds0=ds0, &
            v0=v0, &
            stability_analysis=stability_analysis, maxp=stop_pressure, max_points=max_points)
      elseif (v0==0 .and. p0 == 0) then
         cl = fcritical_line(&
            model=ar_models(id)%model, a0=a0, &
            z0=z0, zi=zi, &
            ns0=ns, S0=S, ds0=ds0, &
            t0=t0, &
            stability_analysis=stability_analysis, maxp=stop_pressure, max_points=max_points)
      else if (t0 == 0 .and. v0==0)  then
         cl = fcritical_line(&
            model=ar_models(id)%model, a0=a0, &
            z0=z0, zi=zi, &
            ns0=ns, S0=S, ds0=ds0, &
            p0=p0, &
            stability_analysis=stability_analysis, maxp=stop_pressure, max_points=max_points)
      else if (t0 == 0) then
         cl = fcritical_line(&
            model=ar_models(id)%model, a0=a0, &
            z0=z0, zi=zi, &
            ns0=ns, S0=S, ds0=ds0, &
            v0=v0, p0=p0, &
            stability_analysis=stability_analysis, maxp=stop_pressure, max_points=max_points)
      else if (v0 == 0) then
         cl = fcritical_line(&
            model=ar_models(id)%model, a0=a0, &
            z0=z0, zi=zi, &
            ns0=ns, S0=S, ds0=ds0, &
            t0=t0, p0=p0, &
            stability_analysis=stability_analysis, maxp=stop_pressure, max_points=max_points)
      else if (p0 == 0) then
         cl = fcritical_line(&
            model=ar_models(id)%model, a0=a0, &
            z0=z0, zi=zi, &
            ns0=ns, S0=S, ds0=ds0, &
            v0=v0, t0=t0, &
            stability_analysis=stability_analysis, maxp=stop_pressure, max_points=max_points)
      else
         cl = fcritical_line(&
            model=ar_models(id)%model, a0=a0, &
            z0=z0, zi=zi, &
            ns0=ns, S0=S, ds0=ds0, &
            v0=v0, t0=t0, p0=p0, &
            stability_analysis=stability_analysis, maxp=stop_pressure, max_points=max_points)

      end if

      do i=1,size(cl%a)
         as(i) = cl%a(i)
         Ts(i) = cl%T(i)
         Ps(i) = cl%P(i)
         Vs(i) = cl%V(i)
      end do

      if (allocated(cl%CEP%y)) then
         CEP_x = cl%CEP%x
         CEP_y = cl%CEP%y
         CEP_P = cl%CEP%P
         CEP_Vx = cl%CEP%Vx
         CEP_Vy = cl%CEP%Vy
         CEP_T = cl%CEP%T
      else
         ! If the critical line is not defined, set the CEP values to NaN
         CEP_x = makenan()
         CEP_y = makenan()
         CEP_P = makenan()
         CEP_Vx = makenan()
         CEP_Vy = makenan()
         CEP_T = makenan()
      end if
   end subroutine critical_line

   subroutine find_llcl(id, z0, zi, P, Tstart, a, T, V)
      use yaeos__equilibria, only: ffind_llcl => find_llcl
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z0(:)
      real(c_double), intent(in) :: zi(:)
      real(c_double), intent(in) :: Tstart
      real(c_double), intent(in) :: P
      real(c_double), intent(out) :: a
      real(c_double), intent(out) :: T
      real(c_double), intent(out) :: V

      T = Tstart
      call ffind_llcl(&
         model=ar_models(id)%model, z0=z0, zi=zi, P=P, &
         a=a, T=T, V=V &
         )
   end subroutine find_llcl

   subroutine stability_zpt(id, z, P, T, w_min, min_tm, all_mins)
      use yaeos, only: min_tpd, tm
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:), P, T
      real(c_double), intent(out) :: w_min(size(z))
      real(c_double), intent(out) :: min_tm
      real(c_double), intent(out) :: all_mins(size(z), size(z)+1)

      real(c_double) :: d_i(size(z))

      integer :: i

      call min_tpd(&
         ar_models(id)%model, z=z, P=P, T=T, &
         mintpd=min_tm, w=w_min, all_minima=all_mins &
         )
   end subroutine stability_zpt

   subroutine stability_zt_ge(id, z, T, w_min, min_tm, all_mins)
      use yaeos, only: min_tpd, tm, pr
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:), T
      real(c_double), intent(out) :: w_min(size(z))
      real(c_double), intent(out) :: min_tm
      real(c_double), intent(out) :: all_mins(size(z), size(z)+1)

      real(c_double) :: d_i(size(z))

      integer :: i

      call min_tpd(&
         ge_models(id)%model, z=z, P=1._pr, T=T, &
         mintpd=min_tm, w=w_min, all_minima=all_mins &
         )
   end subroutine stability_zt_ge

   subroutine tm(id, z, w, P, T, tm_value)
      use yaeos, only: ftm => tm
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:), w(size(z)), P, T
      real(c_double), intent(out) :: tm_value

      tm_value = ftm(model=ar_models(id)%model, z=z, w=w, P=P, T=T)
   end subroutine tm

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
         beta = -1.0
         Vx = 1.0
         Vy = 1.0
         return
      end if

      call equilibria_state_to_arrays(result, x, y, Pout, Tout, Vx, Vy, beta)
   end subroutine flash

   subroutine solve_mp_flash(np, id, z, P, T, kinds_x, kind_w, max_iters, x_l0, w0, betas0, x_l, w, betas, iters, F)
      use yaeos, only: solve_mp_flash_point
      integer(c_int), intent(in) :: np
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:)
      real(c_double), intent(in) :: P
      real(c_double), intent(in) :: T
      integer(c_int), intent(in) :: kinds_x(np)
      integer(c_int), intent(in) :: kind_w
      integer(c_int), intent(in) :: max_iters
      real(c_double), intent(in) :: x_l0(np, size(z))
      real(c_double), intent(in) :: w0(size(z))
      real(c_double), intent(in) :: betas0(np+1)
      real(c_double), intent(out) :: x_l(np, size(z))
      real(c_double), intent(out) :: w(size(z))
      real(c_double), intent(out) :: betas(np+1)
      integer(c_int), intent(out) :: iters
      real(c_double), intent(out) :: F(size(z)*np+np+1+2)

      character(len=14) :: x_kinds(np), w_kind

      real(c_double) :: X(size(z)*np+np+1+2)
      integer :: nc
      integer :: ns1, ns2
      real(c_double) :: S1, S2

      real(c_double) :: dF(size(F), size(F))
      real(c_double) :: K(np, size(z))

      logical :: less_phases
      integer :: beta_0_index
      integer :: i
      real(c_double) :: betas_in(np+1)

      nc = size(z)

      call convert_kind(kinds_x, x_kinds)
      call convert_kind(kind_w, w_kind)

      ! where (betas0 == 0.0)
      !    betas_in = 1e-20
      ! elsewhere
      !    betas_in = betas0
      ! end where
      betas_in = betas0

      X = [(log(x_l0(i, :)/w0), i=1,np), betas_in, log(P), log(T)]

      ns1 = nc*np+np+1+1
      ns2 = nc*np+np+1+2
      S1 = X(ns1)
      S2 = X(ns2)

      call solve_mp_flash_point(&
         ar_models(id)%model, z=z, np=np, kinds_x=x_kinds, kind_w=w_kind, &
         X=X, ns1=ns1, S1=S1, ns2=ns2, S2=S2, max_iters=max_iters, F=F, &
         less_phases=less_phases, beta_0_index=beta_0_index, iters=iters &
         )

      do i=1,np
         K(i, :) = exp(X((i-1)*nc+1 : i*nc))
      end do

      betas = X(np*nc+1 : np*nc+np+1)

      w = z/(matmul(betas(:np), K(:np, :)) + betas(np+1))

      do i=1,np
         x_l(i, :) = K(i, :) * w
      end do

   end subroutine solve_mp_flash

   subroutine flash_vt(id, z, T, V, x, y, k0, Pout, Tout, Vx, Vy, beta)
      use yaeos, only: EquilibriumState, fflash => flash
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:)
      real(c_double), intent(in) :: T
      real(c_double), intent(in) :: V
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
         result = fflash(ar_models(id)%model, z, T, v_spec=V, iters=iters)
      else
         result = fflash(ar_models(id)%model, z, T, v_spec=V, k0=k0, iters=iters)
      end if

      if (.not. allocated(result%x) .or. .not. allocated(result%y)) then
         Pout = -1.0
         Tout = T
         x = z
         y = z
         beta = -1.0
         Vx = 1.0
         Vy = 1.0
         return
      end if

      call equilibria_state_to_arrays(result, x, y, Pout, Tout, Vx, Vy, beta)
   end subroutine flash_vt

   subroutine flash_ge(id, z, T, x, y, k0, Pout, Tout, Vx, Vy, beta)
      use yaeos, only: EquilibriumState, fflash => flash
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:)
      real(c_double), intent(in) :: T
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
         result = fflash(ge_models(id)%model, z, t, iters=iters)
      else
         result = fflash(ge_models(id)%model, z, t, k0=k0, iters=iters)
      end if

      if (.not. allocated(result%x) .or. .not. allocated(result%y)) then
         Tout = T
         x = z
         y = z
         beta = -1
         Vx = 1
         Vy = 1
         return
      end if

      call equilibria_state_to_arrays(result, x, y, Pout, Tout, Vx, Vy, beta)
   end subroutine flash_ge

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

   subroutine saturation_pressure(id, z, T, kind, P0, y0, P, x, y, Vx, Vy, beta)
      use yaeos, only: EquilibriumState, fsaturation_pressure => saturation_pressure
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:)
      real(c_double), intent(in) :: T
      character(len=15), intent(in) :: kind
      real(c_double), intent(in) :: P0
      real(c_double), intent(in) :: y0(size(z))

      real(c_double), intent(out) :: P
      real(c_double), intent(out) :: x(size(z))
      real(c_double), intent(out) :: y(size(z))
      real(c_double), intent(out) :: Vx, Vy, beta

      real(c_double) :: aux

      type(EquilibriumState) :: sat

      if (P0 /= 0 .and. all(y0 /= 0)) then
         sat = fsaturation_pressure(ar_models(id)%model, z, T, kind, P0=P0, y0=y0)
      elseif (P0 /= 0) then
         sat = fsaturation_pressure(ar_models(id)%model, z, T, kind, P0=P0)
      elseif (all(y0 /= 0)) then
         sat = fsaturation_pressure(ar_models(id)%model, z, T, kind, y0=y0)
      else
         sat = fsaturation_pressure(ar_models(id)%model, z, T, kind)
      end if
      call equilibria_state_to_arrays(sat, x, y, P, aux, Vx, Vy, beta)
      if (sat%kind == "failed") then
         ! If the saturation pressure calculation failed, set P to NaN
         P = makenan()
      end if
   end subroutine saturation_pressure

   subroutine saturation_temperature(id, z, P, kind, T0, y0, T, x, y, Vx, Vy, beta)
      use yaeos, only: &
         EquilibriumState, &
         fsaturation_temperature => saturation_temperature, &
         k_wilson
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:)
      real(c_double), intent(in) :: P
      character(len=15), intent(in) :: kind
      real(c_double), intent(in) :: T0
      real(c_double), intent(in) :: y0(size(z))

      real(c_double), intent(out) :: T
      real(c_double), intent(out) :: x(size(z))
      real(c_double), intent(out) :: y(size(z))
      real(c_double), intent(out) :: Vx, Vy, beta

      real(c_double) :: aux

      type(EquilibriumState) :: sat

      if (T0 /= 0 .and. all(y0 /= 0)) then
         sat = fsaturation_temperature(ar_models(id)%model, z, P, kind, y0=y0, T0=T0)
      else if (all(y0 /= 0)) then
         sat = fsaturation_temperature(ar_models(id)%model, z, P, kind, y0=y0)
      else if (T0 /= 0) then
         sat = fsaturation_temperature(ar_models(id)%model, z, P, kind, T0=T0)
      else
         sat = fsaturation_temperature(ar_models(id)%model, z, P, kind)
      end if

      call equilibria_state_to_arrays(sat, x, y, aux, T, Vx, Vy, beta)
      if (sat%kind == "failed") then
         ! If the saturation temperature calculation failed, set T to NaN
         T = makenan()
      end if
   end subroutine saturation_temperature

   subroutine pure_saturation_line(id, comp_id, stop_P, stop_T, P, T, Vx, Vy)
      use yaeos, only: fsat => pure_saturation_line, PurePsat, pr
      integer(c_int), intent(in) :: id
      integer(c_int), intent(in) :: comp_id
      real(c_double), intent(in) :: stop_P
      real(c_double), intent(in) :: stop_T
      real(c_double), intent(out) :: P(10000)
      real(c_double), intent(out) :: T(10000)
      real(c_double), intent(out) :: Vx(10000)
      real(c_double), intent(out) :: Vy(10000)

      integer :: npoints
      type(PurePsat) :: sat

      real(8) :: nan

      nan = 0
      nan = nan/nan

      T = nan
      P = nan
      Vx = nan
      Vy = nan

      sat = fsat(ar_models(id)%model, component=comp_id, minP=stop_P, minT=stop_T)

      npoints = minval([size(sat%T), size(T)])

      T(:npoints) = sat%T(:npoints)
      P(:npoints) = sat%P(:npoints)
      Vx(:npoints) = sat%Vx(:npoints)
      Vy(:npoints) = sat%Vy(:npoints)
   end subroutine pure_saturation_line

   ! ==========================================================================
   ! Multi-phase envelopes
   ! --------------------------------------------------------------------------
   subroutine pt_mp_phase_envelope(&
      id, z, np, x_l0, w0, betas0, P0, T0, ns0, ds0, &
      beta_w, kinds_x, kind_w, max_points, stop_pressure, &
      x_ls, ws, betas, Ps, Ts, iters, ns, main_kinds, ref_kinds, &
      Pcs, Tcs &
      )
      use yaeos, only: PTEnvelMP, pt_envelope
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:)
      integer(c_int), intent(in) :: np
      real(c_double), intent(in) :: x_l0(np, size(z))
      real(c_double), intent(in) :: w0(size(z))
      real(c_double), intent(in) :: betas0(np)
      real(c_double), intent(in) :: P0
      real(c_double), intent(in) :: T0
      integer(c_int), intent(in) :: kinds_x(np)
      integer(c_int), intent(in) :: kind_w
      real(c_double), intent(in) :: beta_w
      real(c_double), intent(in) :: stop_pressure

      integer(c_int), intent(in) :: ns0
      real(c_double), intent(in) :: ds0
      integer(c_int), intent(in) :: max_points

      real(c_double), intent(out) :: x_ls(max_points, np, size(z))
      real(c_double), intent(out) :: ws(max_points, size(z))
      real(c_double), intent(out) :: betas(max_points, np)
      real(c_double), intent(out) :: Ps(max_points)
      real(c_double), intent(out) :: Ts(max_points)

      integer(c_int), intent(out) :: iters(max_points)
      integer(c_int), intent(out) :: ns(max_points)

      character(len=14), intent(out) :: main_kinds(max_points, np)
      character(len=14), intent(out) :: ref_kinds(max_points)
      real(c_double), intent(out) :: Tcs(max_points)
      real(c_double), intent(out) :: Pcs(max_points)


      integer :: i, j

      type(PTEnvelMP) :: pt_mp
      character(len=14) :: x_kinds(np), w_kind

      x_ls = makenan()
      ws = makenan()
      betas = makenan()
      Ps = makenan()
      Ts = makenan()
      Tcs = makenan()
      Pcs = makenan()

      call convert_kind(kinds_x, x_kinds)
      call convert_kind(kind_w, w_kind)

      pt_mp = pt_envelope(&
         model=ar_models(id)%model, np=np, z=z, kinds_x=x_kinds, kind_w=w_kind,&
         x_l0=x_l0, w0=w0, betas0=betas0, &
         P0=P0, T0=T0, ns0=ns0, ds0=ds0, &
         beta_w=beta_w, points=max_points, max_pressure=stop_pressure &
         )

      do i=1,size(pt_mp%points)
         do j=1,np
            x_ls(i, j, :) = pt_mp%points(i)%x_l(j, :)
         end do
         ws(i, :) = pt_mp%points(i)%w
         betas(i, :) = pt_mp%points(i)%betas
         Ps(i) = pt_mp%points(i)%P
         Ts(i) = pt_mp%points(i)%T
         iters(i) = pt_mp%points(i)%iters
         ns(i) = pt_mp%points(i)%ns
         main_kinds(i, :) = pt_mp%points(i)%kinds_x
         ref_kinds(i) = pt_mp%points(i)%kind_w
      end do

      do i=1,size(pt_mp%Tc)
         Tcs(i) = pt_mp%Tc(i)
         Pcs(i) = pt_mp%Pc(i)
      end do
   end subroutine pt_mp_phase_envelope

   subroutine px_mp_phase_envelope(&
      id, z0, zi, np, T, x_l0, w0, betas0, P0, alpha0, ns0, ds0, &
      beta_w, kinds_x, kind_w, max_points, &
      x_ls, ws, betas, Ps, alphas, iters, ns, main_kinds, ref_kinds, Pcs, acs &
      )
      use yaeos, only: PXEnvelMP, px_envelope
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z0(:)
      real(c_double), intent(in) :: zi(:)
      integer(c_int), intent(in) :: np
      real(c_double), intent(in) :: T
      real(c_double), intent(in) :: x_l0(np, size(z0))
      real(c_double), intent(in) :: w0(size(z0))
      real(c_double), intent(in) :: betas0(np)
      real(c_double), intent(in) :: P0
      real(c_double), intent(in) :: alpha0
      real(c_double), intent(in) :: beta_w

      integer(c_int), intent(in) :: ns0
      real(c_double), intent(in) :: ds0
      integer(c_int), intent(in) :: kinds_x(np)
      integer(c_int), intent(in) :: kind_w
      integer(c_int), intent(in) :: max_points

      real(c_double), intent(out) :: x_ls(max_points, np, size(z0))
      real(c_double), intent(out) :: ws(max_points, size(z0))
      real(c_double), intent(out) :: betas(max_points, np)
      real(c_double), intent(out) :: Ps(max_points)
      real(c_double), intent(out) :: alphas(max_points)

      integer(c_int), intent(out) :: iters(max_points)
      integer(c_int), intent(out) :: ns(max_points)

      character(len=14), intent(out) :: main_kinds(max_points, np)
      character(len=14), intent(out) :: ref_kinds(max_points)
      real(c_double), intent(out) :: Pcs(max_points)
      real(c_double), intent(out) :: acs(max_points)

      character(len=14) :: x_kinds(np), w_kind

      integer :: i, j

      type(PXEnvelMP) :: px_mp

      x_ls = makenan()
      ws = makenan()
      betas = makenan()
      Ps = makenan()
      alphas = makenan()
      Pcs = makenan()
      acs = makenan()

      call convert_kind(kinds_x, x_kinds)
      call convert_kind(kind_w, w_kind)

      px_mp = px_envelope(&
         model=ar_models(id)%model, np=np, z0=z0, zi=zi, T=T, x_l0=x_l0, &
         w0=w0, betas0=betas0, P0=P0, alpha0=alpha0, ns0=ns0, ds0=ds0, &
         beta_w=beta_w, kinds_x=x_kinds, kind_w=w_kind, points=max_points &
         )

      do i=1,size(px_mp%points)
         do j=1,np
            x_ls(i, j, :) = px_mp%points(i)%x_l(j, :)
         end do
         ws(i, :) = px_mp%points(i)%w
         betas(i, :) = px_mp%points(i)%betas
         Ps(i) = px_mp%points(i)%P
         alphas(i) = px_mp%alpha(i)
         iters = px_mp%points(i)%iters
         ns = px_mp%points(i)%ns
         main_kinds(i, :) = px_mp%points(i)%kinds_x
         ref_kinds(i) = px_mp%points(i)%kind_w
      end do

      do i=1,size(px_mp%Pc)
         Pcs(i) = px_mp%Pc(i)
         acs(i) = px_mp%ac(i)
      end do
   end subroutine px_mp_phase_envelope

   subroutine tx_mp_phase_envelope(&
      id, z0, zi, np, P, beta_w, kinds_x, kind_w, x_l0, w0, betas0, T0, alpha0, ns0, ds0, max_points, &
      x_ls, ws, betas, Ts, alphas, iters, ns, main_kinds, ref_kinds, Tcs, acs &
      )
      use yaeos, only: TXEnvelMP, tx_envelope
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z0(:)
      real(c_double), intent(in) :: zi(:)
      integer(c_int), intent(in) :: np
      real(c_double), intent(in) :: P
      real(c_double), intent(in) :: beta_w
      integer(c_int), intent(in) :: kinds_x(np)
      integer(c_int), intent(in) :: kind_w
      real(c_double), intent(in) :: x_l0(np, size(z0))
      real(c_double), intent(in) :: w0(size(z0))
      real(c_double), intent(in) :: betas0(np)
      real(c_double), intent(in) :: T0
      real(c_double), intent(in) :: alpha0

      integer(c_int), intent(in) :: ns0
      real(c_double), intent(in) :: ds0
      integer(c_int), intent(in) :: max_points

      real(c_double), intent(out) :: x_ls(max_points, np, size(z0))
      real(c_double), intent(out) :: ws(max_points, size(z0))
      real(c_double), intent(out) :: betas(max_points, np)
      real(c_double), intent(out) :: Ts(max_points)
      real(c_double), intent(out) :: alphas(max_points)

      integer(c_int), intent(out) :: iters(max_points)
      integer(c_int), intent(out) :: ns(max_points)
      character(len=14), intent(out) :: main_kinds(max_points, np)
      character(len=14), intent(out) :: ref_kinds(max_points)
      real(c_double), intent(out) :: Tcs(max_points)
      real(c_double), intent(out) :: acs(max_points)

      integer :: i, j

      type(TXEnvelMP) :: tx_mp
      character(len=14) :: x_kinds(np), w_kind

      x_ls = makenan()
      ws = makenan()
      betas = makenan()
      Ts = makenan()
      alphas = makenan()
      Tcs = makenan()
      acs = makenan()

      call convert_kind(kinds_x, x_kinds)
      call convert_kind(kind_w, w_kind)

      tx_mp = tx_envelope(&
         model=ar_models(id)%model, kinds_x=x_kinds, kind_w=w_kind, &
         np=np, z0=z0, zi=zi, P=P, x_l0=x_l0, &
         w0=w0, betas0=betas0, T0=T0, alpha0=alpha0, ns0=ns0, ds0=ds0, &
         beta_w=beta_w, points=max_points &
         )

      do i=1,size(tx_mp%points)
         do j=1,np
            x_ls(i, j, :) = tx_mp%points(i)%x_l(j, :)
         end do
         ws(i, :) = tx_mp%points(i)%w
         betas(i, :) = tx_mp%points(i)%betas
         Ts(i) = tx_mp%points(i)%T
         alphas(i) = tx_mp%alpha(i)
         iters = tx_mp%points(i)%iters
         ns = tx_mp%points(i)%ns
         main_kinds(i, :) = tx_mp%points(i)%kinds_x
         ref_kinds(i) = tx_mp%points(i)%kind_w
      end do

      do i=1,size(tx_mp%Tc)
         Tcs(i) = tx_mp%Tc(i)
         acs(i) = tx_mp%ac(i)
      end do
   end subroutine tx_mp_phase_envelope

   ! ==========================================================================
   ! Generalized isolines
   ! --------------------------------------------------------------------------
   subroutine generalized_isopleth(&
      id, nc, np, nstab, kinds_x, kind_w, z, x_l0, w0, betas0, P0, T0, &
      spec_variable, spec_variable_value, ns0, S0, dS0, &
      ws_stab, &
      x_ls, ws, betas, Ts, Ps, w_more_stable, found_unstability, max_points &
      )
      !! Calculates a generalized isopleth line for a given model and conditions.
      !! This subroutine computes the generalized isopleth line for a given model and conditions,
      use yaeos, only: GeneralizedIsoZLine, create_generalized_isoz_line
      integer(c_int), intent(in) :: id
      integer(c_int), intent(in) :: nc
      integer(c_int), intent(in) :: np
      integer(c_int), intent(in) :: nstab
      real(c_double), intent(in) :: z(nc)
      real(c_double), intent(in) :: x_l0(np, nc), w0(nc), betas0(np+1), P0, T0
      integer(c_int), intent(in) :: kinds_x(np), kind_w
      integer, intent(in) :: spec_variable
      real(c_double), intent(in) :: spec_variable_value
      integer, intent(in) :: ns0
      real(c_double), intent(in) :: S0, dS0
      real(c_double), intent(in) :: ws_stab(nstab, nc)
      integer(c_int), intent(in) :: max_points

      real(c_double), intent(out) :: x_ls(max_points, np, nc)
      real(c_double), intent(out) :: ws(max_points, nc)
      real(c_double), intent(out) :: betas(max_points, np+1)
      real(c_double), intent(out) :: Ts(max_points)
      real(c_double), intent(out) :: Ps(max_points)
      real(c_double), intent(out) :: w_more_stable(nc)
      logical, intent(out) :: found_unstability

      character(len=14) :: kx(np), kw

      type(GeneralizedIsoZLine) :: line
      integer :: i, l
      
      call convert_kind(kinds_x, kx)
      call convert_kind(kind_w, kw)

      line = create_generalized_isoz_line(&
         ar_models(id)%model, nc=nc, np=np, nstab=nstab, kinds_x=kx, kind_w=kw, z=z, &
         x_l0=x_l0, w0=w0, betas0=betas0, P0=P0, T0=T0, &
         spec_variable=spec_variable, spec_variable_value=spec_variable_value, ns0=ns0, S0=S0, dS0=dS0, &
         ws_stab=ws_stab, max_points=max_points &
         )

      found_unstability = line%found_unstability


      x_ls = makenan()
      ws = makenan()
      Ts = makenan()
      Ps = makenan()
      betas = makenan()
      w_more_stable = makenan()

      do i=1,min(size(line%points), max_points)
         x_ls(i, :, :) = line%points(i)%x_l(:, :)
         ws(i, :) = line%points(i)%w
         Ts(i) = line%points(i)%T
         Ps(i) = line%points(i)%P
         betas(i, :) = line%points(i)%betas
      end do

      if (found_unstability) then
         w_more_stable = line%w_more_stable
      end if

   end subroutine generalized_isopleth



   ! ==========================================================================
   ! Auxiliar
   ! --------------------------------------------------------------------------
   function makenan()
      real(c_double) :: makenan
      makenan = 0
      makenan = makenan/makenan
   end function makenan

   elemental subroutine convert_kind(numeric, string)
      use yaeos, only: root_kinds
      integer(c_int), intent(in) :: numeric
      character(len=14), intent(out) :: string

      select case(numeric)
       case (root_kinds%liquid)
         string = "liquid"
       case (root_kinds%vapor)
         string = "vapor"
       case (root_kinds%stable)
         string = "stable"
       case default
         string = "unknown"
      end select
   end subroutine convert_kind

   subroutine convert_gerg_components(string, numeric)
      !! Converts a string representation of a component to its numeric identifier
      use yaeos, only: G2008Components
      use iso_fortran_env, only: error_unit

      character(len=*), intent(in) :: string !! String representation of the component
      integer(c_int), intent(out) :: numeric !! Numeric identifier for the component

      select case (string)
       case ("methane", "C1", "CH4", "nC1")
         numeric = G2008Components%methane
       case ("nitrogen", "N2")
         numeric = G2008Components%nitrogen
       case ("carbon_dioxide", "CO2", "carbon dioxide")
         numeric = G2008Components%carbon_dioxide
       case ("ethane", "C2", "nC2")
         numeric = G2008Components%ethane
       case ("propane", "C3", "nC3")
         numeric = G2008Components%propane
       case ("nbutane", "C4", "nC4")
         numeric = G2008Components%nbutane
       case ("isobutane", "iC4")
         numeric = G2008Components%isobutane
       case ("npentane", "C5", "nC5")
         numeric = G2008Components%npentane
       case ("isopentane", "iC5")
         numeric = G2008Components%isopentane
       case ("nhexane", "C6", "nC6")
         numeric = G2008Components%nhexane
       case ("nheptane", "C7", "nC7")
         numeric = G2008Components%nheptane
       case ("noctane", "C8", "nC8")
         numeric = G2008Components%noctane
       case ("nonane", "C9", "nC9")
         numeric = G2008Components%nonane
       case ("decane", "C10", "nC10")
         numeric = G2008Components%decane
       case ("hydrogen", "H2")
         numeric = G2008Components%hydrogen
       case ("oxygen", "O2")
         numeric = G2008Components%oxygen
       case ("carbon_monoxide", "CO")
         numeric = G2008Components%carbon_monoxide
       case ("water", "H2O")
         numeric = G2008Components%water
       case ("hydrogen_sulfide", "H2S")
         numeric = G2008Components%hydrogen_sulfide
       case ("helium", "He")
         numeric = G2008Components%helium
       case ("argon", "Ar")
         numeric = G2008Components%argon
       case default
         write(error_unit,*) "Unknown component: ", trim(string)
         numeric = -1
      end select
   end subroutine convert_gerg_components

   subroutine find_self_intersections(x, y, x_inter, y_inter, i, j)
      use yaeos__math, only: intersect_one_line, Point
      integer, parameter :: max_intersections = 10
      real(c_double), intent(in) :: x(:), y(:)
      real(c_double), intent(out) :: x_inter(max_intersections), y_inter(max_intersections)
      integer(c_int), intent(out) :: i(max_intersections), j(max_intersections)

      integer :: k

      type(Point), allocatable :: intersections(:)

      intersections = intersect_one_line(x, y)

      x_inter = makenan()
      y_inter = makenan()

      ! Populate the output arrays
      do concurrent (k=1:min(max_intersections, size(intersections)))
         x_inter(k) = intersections(k)%x
         y_inter(k) = intersections(k)%y
         i(k) = intersections(k)%i
         j(k) = intersections(k)%j
      end do
   end subroutine find_self_intersections

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
end module yaeos_c
