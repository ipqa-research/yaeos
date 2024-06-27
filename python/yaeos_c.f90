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

   public :: pr76, srk, fug_vt

   type :: ArModelContainer
      class(ArModel), allocatable :: model
   end type
   
   type :: GeModelContainer
      class(GeModel), allocatable :: model
   end type

   ! Singleton to hold temporal ArModels
   class(ArModel), allocatable :: model

   ! Containers of models
   integer, parameter :: max_models  = 100000
   logical :: free_ar_model(max_models) = .true.
   logical :: free_ge_model(max_models) = .true.

   class(ArModelContainer), allocatable :: ar_models(:)
   class(GeModelContainer), allocatable :: ge_models(:)


contains

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
            call move_alloc(model, ar_models(i)%model)
            exit
         end if
      end do
      if (id == max_models) error stop 1
   end subroutine
   
   subroutine make_available_ar_models_list(id)
      !! Make the ArModel id available for allocation
      integer(c_int), intent(out) :: id
      integer :: i
      free_ar_model(i) = .true.
   end subroutine

   subroutine pr76(tc, pc, w, kij, lij, id) bind(C, name="PR76")
      real(c_double), intent(in) :: tc(:), pc(:), w(:), kij(:, :), lij(:, :)
      integer(c_int), intent(out) :: id
      
      model = PengRobinson76(tc, pc, w, kij, lij)
      call extend_ar_models_list(id)
   end subroutine
   
   subroutine srk(tc, pc, w, kij, lij, id)
      real(c_double), intent(in) :: tc(:), pc(:), w(:), kij(:, :), lij(:, :)
      integer(c_int), intent(out) :: id
      model = SoaveRedlichKwong(tc, pc, w, kij, lij)
      call extend_ar_models_list(id)
   end subroutine

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
   end subroutine
end module
