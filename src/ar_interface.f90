module ar_interface
   !! # Legacy residual Helmholtz callbacks
   !! Procedure-pointer interfaces used by the deprecated global-state EOS code.
   !!
   !! # Description
   !! `ar_fun` points to a residual Helmholtz routine and `vinit` points to a
   !! liquid-volume initializer. New applications should use concrete `ArModel`
   !! types instead of assigning these global pointers.
   !!
   !! # Examples
   !! ```fortran
   !! ar_fun => my_residual_helmholtz
   !! vinit => my_volume_initializer
   !! ```
   use yaeos__constants, only: pr, R
   use iso_fortran_env, only: error_unit

   implicit none

   procedure(Ares), pointer :: ar_fun !! Legacy residual Helmholtz callback pointer.
   procedure(initial_volume), pointer :: vinit !! Legacy initial-volume callback pointer.

   abstract interface
      subroutine Ares(z, v, t, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
         !! # Ares
         !! Callback interface for legacy residual Helmholtz models.
         !!
         !! # Description
         !! Calculates the residual Helmholtz energy and derivative outputs
         !! consumed by the legacy property routines.
         import pr
         real(pr), intent(in) :: z(:) !! Component mole amounts.
         real(pr), intent(in) :: v !! Total volume [L].
         real(pr), intent(in) :: t !! Temperature [K].
         real(pr), intent(out) :: Ar !! Residual Helmholtz energy.
         real(pr), intent(out) :: ArV !! Volume derivative of residual Helmholtz energy.
         real(pr), intent(out) :: ArTV !! Mixed temperature/volume derivative of residual Helmholtz energy.
         real(pr), intent(out) :: ArV2 !! Second volume derivative of residual Helmholtz energy.
         real(pr), dimension(size(z)), intent(out) :: Arn !! Mole-number gradient of residual Helmholtz energy.
         real(pr), dimension(size(z)), intent(out) :: ArVn !! Mixed volume/mole-number derivative of residual Helmholtz energy.
         real(pr), dimension(size(z)), intent(out) :: ArTn !! Mixed temperature/mole-number derivative of residual Helmholtz energy.
         real(pr), intent(out) :: Arn2(size(z), size(z)) !! Composition Hessian.
      end subroutine

      function initial_volume(z, p, t)
         !! # initial_volume
         !! Return an initial volume estimate for a legacy pressure solve.
         import pr
         real(pr) :: z(:) !! Component mole amounts.
         real(pr) :: p !! Pressure [bar].
         real(pr) :: t !! Temperature [K].
         real(pr) :: initial_volume !! Initial total volume [L].
      end function
   end interface
end module
