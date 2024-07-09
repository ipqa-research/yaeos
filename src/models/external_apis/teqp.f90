module yaeos__models_external_apis_feos
   !! # `teqp` API
   !! API to work with the `teqp` rust library models.
   !!
   !! # Description
   !! Work in progress
   !!
   !! # Examples
   !!
   !! ```fortran
   !!  A basic code example
   !! ```
   !!
   !! # References
   !! [teqp](https://github.com/ianhbell/teqp)
   use yaeos__models_ar, only: ArModel
   use yaeos__constants, only: pr
   implicit none

   type, abstract, extends(ArModel) :: ArModelteqp
      real(pr), allocatable :: parameters(:) 
         !! Relevant model parameteres, should be specific from model
   contains
   end type

contains

   subroutine model_setter(self, setup_string)
      !! Setup the model from a string
      class(ArModelteqp), intent(out) :: self
      character(len=*), intent(in) :: setup_string
      
      ! Setup logic from some input string? 
      ! I'd prefer a more specific implementation but I'm not sure how easy
      ! to extend it could be
   end subroutine

   subroutine residual_helmholtz(&
            self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2 &
         )
         !! Residual Helmholtz model generic interface.
         !!
         !! This interface represents how an Ar model should be implemented.
         !! By our standard, a Resiudal Helmholtz model takes as input:
         !!
         !! - The mixture's number of moles vector.
         !! - Volume, by default in liters.
         !! - Temperature, by default in Kelvin.
         !!
         !! All the output arguments are optional. While this keeps a long
         !! signature for the implementation, this is done this way to take
         !! advantage of any inner optimizations to calculate derivatives
         !! inside the procedure.
         !!
         !! Once the model is implemented, the signature can be short like
         !! `model%residual_helmholtz(n, v, t, ArT2=dArdT2)`
         class(ArModelteqp), intent(in) :: self !! ArModel
         real(pr), intent(in) :: n(:) !! Moles vector
         real(pr), intent(in) :: v !! Volume [L]
         real(pr), intent(in) :: t !! Temperature [K]
         real(pr), optional, intent(out) :: Ar !! Residual Helmoltz energy
         real(pr), optional, intent(out) :: ArV !! \(\frac{dAr}{dV}\)
         real(pr), optional, intent(out) :: ArT !! \(\frac{dAr}{dT}\)
         real(pr), optional, intent(out) :: ArT2 !! \(\frac{d^2Ar}{dT^2}\)
         real(pr), optional, intent(out) :: ArTV !! \(\frac{d^2Ar}{dTV}\)
         real(pr), optional, intent(out) :: ArV2 !! \(\frac{d^2Ar}{dV^2}\)
         real(pr), optional, intent(out) :: Arn(size(n)) !! \(\frac{dAr}{dn_i}\)
         real(pr), optional, intent(out) :: ArVn(size(n)) !! \(\frac{d^2Ar}{dVn_i}\)
         real(pr), optional, intent(out) :: ArTn(size(n)) !! \(\frac{d^2Ar}{dTn_i}\)
         real(pr), optional, intent(out) :: Arn2(size(n), size(n))!! \(\frac{d^2Ar}{dn_{ij}}\)

         ! Here the function that calculates residual helmholtz energy
         ! from teqp would be called.
   end subroutine
end module