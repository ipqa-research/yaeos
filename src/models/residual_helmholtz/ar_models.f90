module yaeos_models_ar
    !! Module that defines the basics of a residual Helmholtz energy.
    use yaeos_constants, only: pr
    implicit none

    type, abstract :: ArModel
        !! Abstract residual Helmholtz model.
        integer :: id
        character(len=:), allocatable :: name
        procedure(abs_residual_helmholtz), pointer :: residual_helmholtz => NULL()
            !! Method to calculate residual helmholtz energy and derivatives.
        procedure(abs_volume_initializer), pointer :: get_v0 => NULL()
            !! Volume initializer
    end type

    abstract interface
       pure subroutine abs_residual_helmholtz( &
             self, n, v, t, &
             ar, dardn, dardv, dardt, &
             dardn2, darnv, dardnt, dardvt, dardv2, dardt2 &
          )
          import pr
          import ArModel
          class(ArModel), intent(in) :: self
          real(pr), intent(in) :: n(:), v, t
          real(pr), optional, intent(out) :: ar, dardn(size(n)), dardv, dardt
          real(pr), optional, intent(out) :: dardn2(size(n), size(n)), darnv(size(n)), dardnt(size(n))
          real(pr), optional, intent(out) :: dardvt, dardv2, dardt2
       end subroutine

       pure function abs_volume_initializer(self, n, p, t) result(v)
           import pr
           import ArModel
           class(ArModel), intent(in) :: self
           real(pr), intent(in) :: n(:)
           real(pr), intent(in) :: p, t
           real(pr) :: v
       end function
    end interface
end module