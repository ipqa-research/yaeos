module ar_models
   !! Simple derived types that hold the important parameters
   use constants, only: pr
   use hyperdual_mod

   implicit none

   type :: ArModel
      !! Residual Helmholtz model, this type contains the basics attributes that
      !! this kind of model must have. 
      integer :: size
      character(len=30), allocatable :: names(:) !! Components names
      real(pr), allocatable :: z(:) !! Global composition
      character(len=:), allocatable :: thermo_model !! Thermodynamic model name
      character(len=:), allocatable :: mixing_rule !! Mixing rule to use name
      procedure(dual_property), nopass, pointer :: residual_helmholtz
         !! Ar(model, z, v, t): Residual Helmholtz function
   end type

   interface size
      module procedure :: size_model
   end interface

   interface alloc
      module procedure :: alloc_ArModel
   end interface
   
   ! Dual Property
   abstract interface
     pure subroutine dual_property(model, z, v, t, property)
        import hyperdual
        import ArModel
        class(ArModel), intent(in) :: model
        type(hyperdual), intent(in) :: z(model%size), v, t
        type(hyperdual), intent(out) :: property
     end subroutine dual_property
   end interface

contains
   ! ===========================================================================
   ! Constructor
   ! ---------------------------------------------------------------------------
   subroutine init_ArModel(model, ares)
      class(ArModel) :: model
      procedure(dual_property) :: ares

      model%residual_helmholtz => ares
   end subroutine init_ArModel

   subroutine alloc_ArModel(model, n)
      type(ArModel) :: model
      integer :: n, stat

      model%size = n

      allocate(model%names(n), stat=stat)
      allocate(model%z(n), stat=stat)
   end subroutine
   ! ===========================================================================

   ! ===========================================================================
   ! Thermo properties
   ! ---------------------------------------------------------------------------
   subroutine residual_helmholtz(model, z, v, t, ar, dar, dar2)
      class(ArModel) :: model
      real(pr), intent(in) :: z(model%size)
      real(pr), intent(in) :: v, t

      real(pr), intent(out) :: ar
      real(pr), intent(out) :: dar(model%size + 2)
      real(pr), intent(out) :: dar2(model%size + 2, model%size + 2)

      ar = 0
      dar = 0
      dar2 = 0

      call dualderiv( &
         model, z, v, t, &
         model%residual_helmholtz, ar, dar, dar2 &
      )
   end subroutine residual_helmholtz
   ! ===========================================================================
   
   pure function size_model(model) result(smodel)
      class(ArModel), intent(in) :: model
      integer :: smodel

      smodel = size(model%names)
   end function

   pure subroutine dualderiv( &
      model, z, v, t, &
      f_in, &
      f, df, df2)
      class(ArModel), intent(in) :: model

      real(pr), intent(in) :: v, t, z(model%size)
      procedure(dual_property) :: f_in

      real(pr), intent(out) :: f
      real(pr), intent(out) :: df(model%size + 2)
      real(pr), intent(out) :: df2(model%size + 2, model%size + 2)

      type(hyperdual) :: X(model%size + 2)
      type(hyperdual) :: y

      integer :: i, j, n

      n = model%size + 2
      df = 0
      df2 = 0

      X = [z, v, t]
      do i = 1, n
         X = [z, v, t]
         X(i)%f1 = 1
         X(i)%f2 = 1

         call f_in(model, X(:n - 2), X(n - 1), X(n), y)
         df(i) = y%f1
         df2(i, i) = y%f12

         do j = i, n
            X = [z, v, t]
            X(i)%f1 = 1
            X(j)%f2 = 1

            call f_in(model, X(:n - 2), X(n - 1), X(n), y)
            df2(i, j) = y%f12
            df2(j, i) = df2(i, j)
         end do
      end do

      f = y%f0
   end subroutine
end module ar_models
