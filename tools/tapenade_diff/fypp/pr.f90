! =============================================================================
! -> Setup or select your model here


! =============================================================================
module model
    use yaeos___constants, only: pr
    use yaeos___models_ar, only: ArModel
    implicit none

    type :: PR
    real(8), allocatable :: tc(:), pc(:), w(:)
        real(8), allocatable :: ac(:), b(:)
        real(8), allocatable :: del1(:), del2(:)
        real(8), allocatable :: k(:)
        real(8), allocatable :: kij(:, :), lij(:, :)
    end type

contains

    subroutine ar(self, n, v, t, arval)
        type(PR) :: self
        real(8), intent(in) :: n(:)
        real(8), intent(in) :: v, t
        real(8), intent(out) :: arval
        ! Generic Cubic Variables
        real(8) :: b_v
        ! Alpha function variables
        real(8) :: a(size(n))
        ! QMR Variables
        real(8) :: amix, bmix, d1mix, d2mix
        integer :: i, j
        
        a = self%ac * (1 + self%k*(1 - sqrt(t/self%tc)))**2
        do i=1, size(n)
            do j=1, size(n)
                amix = amix + sqrt(a(i)*a(j)) * (1 - self%kij(i, j))
                bmix = bmix + (self%b(i) + self%b(j))/2 * (1 - self%lij(i, j))
                d1mix = self%delta_1(1)
                d2mix = self%delta_2(1)
            end do
        end do
        b_v = bmix/v
        arval = (&
                - sum(n) * log(1.0 - b_v) &
                - amix / (R*t*bmix)*1.0 / (d1mix - d2mix) &
                * log((1.0 + d1mix * b_v) / (1.0 + d2mix * b_v)) &
        ) * (R * t)
    end subroutine

    pure function volume_initalizer(n, p, t) result(v0)
        real(8), intent(in) :: n(:)
        real(8), intent(in) :: p
        real(8), intent(in) :: t
        real(8) :: v0
        v0 = sum(n*b)/sum(b)
    end function
end module