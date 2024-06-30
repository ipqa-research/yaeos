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
   use yaeos, only: &
   ! Generic Models
   & ArModel, GeModel, &

   ! Cubic Models
   & CubicEoS, &
   & SoaveRedlichKwong, PengRobinson76, PengRobinson78, &
   & CubicMixRule, QMR, MHV, &

   ! Ge Models
   & fNRTL => NRTL, &

   ! Thermodynamic properties
   & fugacity_vt
   implicit none

   private

   public :: pr76
   public :: srk
   public :: fug_vt
   public :: make_available_ar_models_list
   public :: make_available_ge_models_list
   public :: set_mhv, set_qmr
   public :: nrtl

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
      use yaeos, only: QMR
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
   !  Cubic EoS
   ! --------------------------------------------------------------------------
   subroutine pr76(tc, pc, w, id) bind(C, name="PR76")
      real(c_double), intent(in) :: tc(:), pc(:), w(:)
      integer(c_int), intent(out) :: id

      ar_model = PengRobinson76(tc, pc, w)
      call extend_ar_models_list(id)
   end subroutine pr76

   subroutine srk(tc, pc, w, kij, lij, id)
      real(c_double), intent(in) :: tc(:), pc(:), w(:), kij(:, :), lij(:, :)
      integer(c_int), intent(out) :: id
      ar_model = SoaveRedlichKwong(tc, pc, w, kij, lij)
      call extend_ar_models_list(id)
   end subroutine srk

   ! ==========================================================================
   !  Thermodynamic properties
   ! --------------------------------------------------------------------------
   subroutine fug_vt(id, n, v, t, lnfug, dlnphidp, dlnphidt, dlnphidn)
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
end module yaeos_c
