module yaeos_models_ar
    !! Module that defines the basics of a residual Helmholtz energy.
    use yaeos_constants, only: pr
    implicit none

    type, abstract :: ArModel
        !! Abstract residual Helmholtz model.
        integer :: id
        character(len=:), allocatable :: name
    contains
        procedure(abs_residual_helmholtz), deferred :: residual_helmholtz
            !! Method to calculate residual helmholtz energy and derivatives.
        procedure(abs_volume_initializer), deferred :: get_v0
            !! Volume initializer
    end type

    abstract interface
        subroutine abs_residual_helmholtz(&
            self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2 &
        )
         !! Residual Helmholtz model generic interface
         import ArModel, pr
         class(ArModel), intent(in) :: self
         real(pr), intent(in) :: n(:)
         real(pr), intent(in) :: v, t
         real(pr), optional, intent(out) :: Ar, ArV, ArT, ArT2, ArTV, ArV2
         real(pr), optional, dimension(size(n)), intent(out) :: Arn, ArVn, ArTn
         real(pr), optional, intent(out) :: Arn2(size(n), size(n))
      end subroutine

      function abs_volume_initializer(self, n, p, t)
         !! Initialization of volume
         import ArModel, pr
         class(ArModel), intent(in) :: self
         real(pr), intent(in) :: n(:)
         real(pr), intent(in) :: p
         real(pr), intent(in) :: t
         real(pr) :: abs_volume_initializer
      end function
   end interface

end module