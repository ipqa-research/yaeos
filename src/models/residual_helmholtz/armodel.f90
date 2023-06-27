module ar_models
   !! Simple derived types that hold the important parameters
   use constants, only: pr
   use hyperdual_mod
   use yaeos_interfaces, only: dual_property

   implicit none

   private
   public :: set_ar_function, residual_helmholtz

   procedure(dual_property), pointer :: ar_fun

contains
   ! ===========================================================================
   ! Constructor
   ! ---------------------------------------------------------------------------
   subroutine set_ar_function(ar)
      ! Setup the general ar_function
      procedure(dual_property) :: ar
      ar_fun => ar
   end subroutine
   ! ===========================================================================

   ! ===========================================================================
   ! Get all Residual Helmholtz derivatives
   ! ---------------------------------------------------------------------------
   subroutine residual_helmholtz(z, v, t, ar, dar, dar2)
      real(pr), intent(in) :: z(:)
      real(pr), intent(in) :: v, t

      real(pr), intent(out) :: ar
      real(pr), intent(out) :: dar(size(z) + 2)
      real(pr), intent(out) :: dar2(size(z) + 2, size(z) + 2)

      ar = 0
      dar = 0
      dar2 = 0

      call dualderiv( &
         z, v, t, &
         ar_fun, ar, dar, dar2 &
      )
   end subroutine residual_helmholtz
   ! ===========================================================================

   pure subroutine dualderiv( &
      z, v, t, &
      f_in, &
      f, df, df2)

      real(pr), intent(in) :: v, t, z(:)
      procedure(dual_property) :: f_in

      real(pr), intent(out) :: f
      real(pr), intent(out) :: df (size(z) + 2)
      real(pr), intent(out) :: df2(size(z) + 2, size(z) + 2)

      type(hyperdual) :: X(size(z) + 2)
      type(hyperdual) :: y

      integer :: i, j, n

      n = size(z) + 2
      df = 0
      df2 = 0

      X = [z, v, t]
      do i = 1, n
         X = [z, v, t]
         X(i)%f1 = 1
         X(i)%f2 = 1

         call f_in(X(:n - 2), X(n - 1), X(n), y)
         df(i) = y%f1
         df2(i, i) = y%f12
         do j = i, 0
            X = [z, v, t]
            X(i)%f1 = 1
            X(j)%f2 = 1

            call f_in(X(:n - 2), X(n - 1), X(n), y)
            df2(i, j) = y%f12
            df2(j, i) = df2(i, j)
         end do
      end do

      f = y%f0
   end subroutine
end module ar_models
