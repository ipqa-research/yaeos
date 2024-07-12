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
   public :: srk, pr76, pr78
   ! Mixing rules
   public :: set_mhv, set_qmr
   
   ! __del__
   public :: make_available_ar_models_list
   public :: make_available_ge_models_list
   ! GeMoels
   public :: nrtl
   public :: ln_gamma
   
   ! Thermoprops
   public :: fug_vt
   ! Phase equilibria
   public :: flash
   public :: saturation_pressure

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
   end subroutine

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

   subroutine srk(tc, pc, w, kij, lij, id)
      use yaeos, only: SoaveRedlichKwong
      real(c_double), intent(in) :: tc(:), pc(:), w(:), kij(:, :), lij(:, :)
      integer(c_int), intent(out) :: id
      ar_model = SoaveRedlichKwong(tc, pc, w, kij, lij)
      call extend_ar_models_list(id)
   end subroutine srk

   ! ==========================================================================
   !  Thermodynamic properties
   ! --------------------------------------------------------------------------
   subroutine fug_vt(id, n, v, t, lnfug, dlnphidp, dlnphidt, dlnphidn)
      use yaeos, only: fugacity_vt
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: n(:), v, t
      real(c_double), intent(out) :: lnfug(size(n))
      real(c_double) :: p

      real(c_double), optional, intent(in out) :: &
         dlnphidp(size(n)), dlnphidt(size(n)), dlnphidn(size(n), size(n))

      call fugacity_vt(&
         ar_models(id)%model, &
         n, V, T, P, lnfug, dlnPhidP, dlnphidT, dlnPhidn &
         )
   end subroutine fug_vt

   ! ==========================================================================
   ! Phase equilibria
   ! --------------------------------------------------------------------------
   subroutine equilibria_state_to_arrays(eq_state, x, y, P, T, Vx, Vy, beta)
      use yaeos, only: EquilibriaState
      type(EquilibriaState) :: eq_state
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
   end subroutine

   subroutine flash(id, z, T, P, x, y, k0, Pout, Tout, Vx, Vy, beta)
      use yaeos, only: EquilibriaState, fflash => flash
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:)
      real(c_double), intent(in) :: T
      real(c_double), intent(in) :: P
      real(c_double), optional, intent(in out) :: k0(size(z))
      real(c_double), intent(out) :: x(size(z))
      real(c_double), intent(out) :: y(size(z))
      real(c_double), intent(out) :: Pout
      real(c_double), intent(out) :: Tout
      real(c_double), intent(out) :: Vx
      real(c_double), intent(out) :: Vy
      real(c_double), intent(out) :: beta
      
      type(EquilibriaState) :: result
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
   end subroutine

   subroutine saturation_pressure(id, z, T, kind, P, x, y, Vx, Vy, beta)
      use yaeos, only: EquilibriaState, fsaturation_pressure => saturation_pressure
      integer(c_int), intent(in) :: id
      real(c_double), intent(in) :: z(:)
      real(c_double), intent(in) :: T
      character(len=15), intent(in) :: kind

      real(c_double), intent(out) :: P
      real(c_double), intent(out) :: x(size(z))
      real(c_double), intent(out) :: y(size(z))
      real(c_double), intent(out) :: Vx, Vy, beta

      real(c_double) :: aux
      
      type(EquilibriaState) :: sat

      sat = fsaturation_pressure(ar_models(id)%model, z, T, kind)
      call equilibria_state_to_arrays(sat, x, y, P, aux, Vx, Vy, beta)
   end subroutine
end module yaeos_c
