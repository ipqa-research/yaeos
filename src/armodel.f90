module yaeos_ar_models
   !-| 
   ! # Residual free Helmholtz energy models.
   !
   !   This module holds the general logic for residual Helmholtz 
   !   free energy models. This library assumes that always the reduced
   !   \(\alpha^r = \frac{A^r}{RT}\) energy is used.
   !
   !   A Helmholtz free energy model is considered to be a just single
   !   subroutine that respondes to the function \(Ar = f(z, v, t)\).
   !   Using hyper-dual numbers as variables and output.
   !
   !   This function is stored in a single variable `ar_fun` that will later
   !   be called to obtain specific derivatives and, later on, use it inside
   !   other routines to calculate bulk properties.
   !
   !   In principle, just defining a subroutine that responds to the interface
   !   provided by `dual_property` and setting that function as the desired
   !   Ar function to use is enough to provide all the library provided 
   !   functionality.
   !
   ! ```fortran
   ! subroutine Ar(z, v, t, ar)
   !    use yaeos_adiff
   !    type(hyperdual), intent(in) :: z(:)
   !    type(hyperdual), intent(in) :: v
   !    type(hyperdual), intent(in) :: t
   !    type(hyperdual), intent(out) :: ar
   !
   !    ! A very complicated residual helmholtz function of a mixture
   !    ar = sum(z) * v * t
   ! end subroutine
   !
   ! ...
   !
   !program
   ! use yaeos, only: set_ar_function, pressure
   ! call set_ar_function(Ar)
   ! call pressure(z, v, t, p, dp, dp2) ! This calculates the pressure and
   !                                    ! derivatives.
   !end program
   !```
   !

   use yaeos_constants, only: pr
   use yaeos_autodiff
   use yaeos_interfaces, only: dual_property

   implicit none

   private
   public :: set_ar_function, residual_helmholtz, ar_fun

   procedure(dual_property), pointer :: ar_fun

contains
   ! ===========================================================================
   ! Constructor
   ! ---------------------------------------------------------------------------
   subroutine set_ar_function(ar)
      !-| Select which Ar subroutine to use.
      !   this subroutine sets up the pointer procedure to the desired
      !   Helmholtz subroutine
      procedure(dual_property) :: ar !| Ar subroutine to use
      ar_fun => ar
   end subroutine
   ! ===========================================================================

   ! ===========================================================================
   ! Get all Residual Helmholtz derivatives
   ! ---------------------------------------------------------------------------
   subroutine residual_helmholtz(z, v, t, ar, dar, dar2)
      !-| General residual Helmholtz subroutine.
      !
      !   This subroutine handles all the necesary calls to calculate the
      !   residual Helmholtz free energy with just using the pointed Ar
      !   subroutine.
      real(pr), intent(in) :: z(:) !| Composition
      real(pr), intent(in) :: v !| Volume
      real(pr), intent(in) :: t !| Temperature

      real(pr), intent(out) :: ar !| Residual Helmholtz energy
      real(pr), intent(out) :: dar(size(z) + 2) !| First derivatives
      real(pr), intent(out) :: dar2(size(z) + 2, size(z) + 2) !| Second derivatives

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
      !-| Subroutine that takes care of he hyperdual numbers logic to obtain
      !   the function and it's derivatives.
      real(pr), intent(in) :: v, t, z(:)
      procedure(dual_property) :: f_in !| Subroutine to use

      real(pr), intent(out) :: f !| Function value
      real(pr), intent(out) :: df(size(z) + 2) !| Gradient
      real(pr), intent(out) :: df2(size(z) + 2, size(z) + 2) !| Hessian

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
         do j = i, n-1
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
end module
