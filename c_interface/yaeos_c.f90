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
   public :: srk, pr76, pr78, rkpr
   ! Mixing rules
   public :: set_mhv, set_qmr

   ! __del__
   public :: make_available_ar_models_list
   public :: make_available_ge_models_list
   ! GeMoels
   public :: nrtl
   public :: unifac_vle
   public :: uniquac
   public :: ln_gamma

   ! Thermoprops
   public :: lnphi_vt, lnphi_pt, pressure, volume
   
   ! Phase equilibria
   public :: flash, flash_grid
   public :: saturation_pressure
   public :: pt2_phase_envelope

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
   subroutine nrtl(a, b, c, id)
      use yaeos, only: fNRTL => NRTL
      real(c_double), intent(in) :: a(:,:), b(:,:), c(:,:)
      integer(c_int), intent(out) :: id
      ge_model = fNRTL(a, b, c)
      call extend_ge_models_list(id)
   end subroutine nrtl

   subroutine unifac_vle(id, nc, ngs, g_ids, g_v)
      use yaeos, only: UNIFAC, setup_unifac, Groups
      integer(c_int), intent(out) :: id
      integer(c_int), intent(in) :: nc !! Number of components
      integer(c_int), intent(in) :: ngs(nc) !! Number of groups at each molecule
      integer(c_int), intent(in) :: g_ids(:, :) !! Ids of groups for each molecule
      integer(c_int), intent(in) :: g_v(:, :) !! Number of groups for each molecule

      type(Groups) :: molecules(nc)
      integer :: i

      do i=1,nc
         molecules(i)%groups_ids = g_ids(i, :ngs(i))
         molecules(i)%number_of_groups = g_v(i, :ngs(i))
      end do

      ge_model = setup_unifac(molecules)
      call extend_ge_models_list(id)
   end subroutine

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
   end subroutine

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

   subroutine ln_gamma(id, n, T, lngamma)
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: n(:)
      real(c_double), intent(in) :: T
      real(c_double), intent(out) :: lngamma(size(n))
      call ge_models(id)%model%ln_activity_coefficient(n, T, lngamma)
   end subroutine ln_gamma

   ! =============================================================================
   !  Ar Models
   ! -----------------------------------------------------------------------------
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

      if (present(delta_1) .and. present(k)) then
         ar_model = fRKPR(tc, pc, w, zc, delta_1=delta_1, k=k)
      elseif (present(delta_1))  then
         ar_model = fRKPR(tc, pc, w, zc, delta_1=delta_1)
      elseif (present(k))  then
         ar_model = fRKPR(tc, pc, w, zc, k=k)
      else
         ar_model = fRKPR(tc, pc, w, zc)
      end if
      call extend_ar_models_list(id)
   end subroutine rkpr

   ! ==========================================================================
   !  Thermodynamic properties
   ! --------------------------------------------------------------------------
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

   ! ==========================================================================
   ! Phase equilibria
   ! --------------------------------------------------------------------------
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

   subroutine flash(id, z, T, P, x, y, Pout, Tout, Vx, Vy, beta)
      use yaeos, only: EquilibriumState, fflash => flash
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:)
      real(c_double), intent(in) :: T
      real(c_double), intent(in) :: P
      real(c_double), intent(out) :: x(size(z))
      real(c_double), intent(out) :: y(size(z))
      real(c_double), intent(out) :: Pout
      real(c_double), intent(out) :: Tout
      real(c_double), intent(out) :: Vx
      real(c_double), intent(out) :: Vy
      real(c_double), intent(out) :: beta

      type(EquilibriumState) :: result
      integer :: iters

      result = fflash(ar_models(id)%model, z, t, p_spec=p, iters=iters)
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

   subroutine saturation_pressure(id, z, T, kind, P, x, y, Vx, Vy, beta)
      use yaeos, only: EquilibriumState, fsaturation_pressure => saturation_pressure
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:)
      real(c_double), intent(in) :: T
      character(len=15), intent(in) :: kind

      real(c_double), intent(out) :: P
      real(c_double), intent(out) :: x(size(z))
      real(c_double), intent(out) :: y(size(z))
      real(c_double), intent(out) :: Vx, Vy, beta

      real(c_double) :: aux

      type(EquilibriumState) :: sat

      sat = fsaturation_pressure(ar_models(id)%model, z, T, kind)
      call equilibria_state_to_arrays(sat, x, y, P, aux, Vx, Vy, beta)
   end subroutine saturation_pressure

   subroutine pt2_phase_envelope(id, z, kind, max_points, Ts, Ps, tcs, pcs, T0, P0)
      use yaeos, only: &
         saturation_pressure, saturation_temperature, pt_envelope_2ph, &
         EquilibriumState, PTEnvel2
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:)
      integer, intent(in) :: max_points
      character(len=15), intent(in) :: kind
      real(c_double), intent(out) :: Ts(max_points)
      real(c_double), intent(out) :: Ps(max_points)
      real(c_double), intent(out) :: Tcs(5), Pcs(5)
      real(c_double), optional, intent(in) :: T0, P0

      real(8) :: makenan, nan
      type(EquilibriumState) :: sat
      type(PTEnvel2) :: env

      integer :: i, neval=0

      real(c_double) :: T, P

      makenan=0

      neval = neval + 1
      nan = makenan/makenan
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
       case("dew")
         sat = saturation_temperature(ar_models(id)%model, z, P=P, kind=kind)
       case("liquid-liquid")
         sat = saturation_temperature(ar_models(id)%model, z, P=P, kind=kind)
      end select


      env = pt_envelope_2ph(ar_models(id)%model, z, sat, points=max_points)
      i = size(env%points)
      Ts(:i) = env%points%T
      Ps(:i) = env%points%P

      i = size(env%cps)
      Tcs(:i) = env%cps%T
      Pcs(:i) = env%cps%P
   end subroutine pt2_phase_envelope

   subroutine flash_grid(id, z, Ts, Ps, xs, ys, Vxs, Vys, betas)
      use yaeos, only: EquilibriumState, flash
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:)
      real(c_double), intent(in) :: Ts(:)
      real(c_double), intent(in) :: Ps(:)
      real(c_double), dimension(size(Ps), size(Ts), size(z)), intent(out) :: xs, ys
      real(c_double), dimension(size(Ps), size(Ts)), intent(out) :: Vxs, Vys, betas

      class(ArModel), allocatable :: model
      type(EquilibriumState) :: flash_result

      real(8) :: T, P

      integer :: i, j, nt, np, iter

      model = ar_models(id)%model
      np = size(Ps)
      nt = size(Ts)

      !$OMP PARALLEL DO PRIVATE(i, j, t, p, flash_result) shared(model, z, ts, ps, betas, Vxs, Vys, xs, ys)
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
   end subroutine
end module yaeos_c
