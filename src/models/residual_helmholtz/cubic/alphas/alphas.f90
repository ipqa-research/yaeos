module yaeos_models_ar_cubic_alphas
    use yaeos_constants, only: pr
    use yaeos_substance, only: substances
    use yaeos_models_ar_genericcubic, only: CubicEoS, AlphaFunction
    implicit none

    type, extends(AlphaFunction) :: AlphaSoave
        real(pr), allocatable :: k(:)
    contains
        procedure :: alpha
    end type

contains

    subroutine alpha(self, Tr, a, dadt, dadt2)
        class(AlphaSoave), intent(in) :: self
        real(pr), intent(in) :: Tr(:)
        real(pr), intent(out) :: a(:), dadt(:), dadt2(:)

        associate(k => self%k)
            a = (1 + k*(1 - sqrt(Tr)))**2
            dadT = k*(k*(sqrt(Tr) - 1) - 1)/sqrt(Tr)
            dadT2 = (1.0d0/2.0d0)*k*(k + 1)/Tr**(1.5_pr)
        end associate

    end subroutine
end module