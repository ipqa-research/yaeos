module tapenade_pr
    use yaeos_tapenade_ar_api, only: ArModelTapenade
    implicit none

    real(8), allocatable :: kij(:, :), lij(:, :)
    real(8), allocatable :: ac(:), b(:), k(:)
    real(8), allocatable :: tc(:), pc(:), w(:)
    real(8), parameter :: r=0.08314472
    real(8), parameter :: del1 = 1. + sqrt(2.)
    real(8), parameter :: del2 = 1. - sqrt(2.)

    type(ArModelTapenade) :: model
contains
    subroutine setup_model(tc_in, pc_in, w_in, kij_in, lij_in)
        real(8) :: tc_in(:)
        real(8) :: pc_in(:)
        real(8) :: w_in(:)
        real(8) :: kij_in(:, :)
        real(8) :: lij_in(:, :)

        tc = tc_in
        pc = pc_in
        w = w_in

        ac = 0.45723553 * R**2 * tc**2 / pc
        b = 0.07779607 * R * tc/pc
        k = 0.37464 + 1.54226 * w - 0.26993 * w**2

        kij = kij_in
        lij = lij_in

        model%ar => ar
        model%ar_d => ar_d
        model%ar_b => ar_b
        model%ar_d_b => ar_d_b
        model%ar_d_d => ar_d_d
    end subroutine

    subroutine ar(n, v, t, arval)
        real(8), intent(in) :: n(:), v, t
        real(8), intent(out) :: arval
        real(8) :: amix, a(size(n)), ai(size(n)), z2(size(n)), nij
        real(8) :: bmix
        real(8) :: b_v
        real(8) :: aij(size(n), size(n)), bij(size(n), size(n))
        integer :: i, j

        a = sqrt(ac * (1.0 + k * (1.0 - sqrt(t/tc)))**2)
        
        amix = 0.0
        bmix = 0.0

        do i=1,size(n)-1
            do j=i+1,size(n)
                nij = n(i) * n(j)
                amix = amix + 2 * nij * (a(i) * a(j)) * (1 - kij(i, j))
                bmix = bmix + nij * (b(i) + b(j)) * (1 - lij(i, j))
             end do
        end do

        amix = amix + sum(n**2*a**2)
        bmix = bmix + sum(n**2 * b)

        bmix = bmix/sum(n)

        b_v = bmix/v
        arval = (&
              - sum(n) * log(1.0 - b_v) &
              - amix / (R*t*bmix)*1.0 / (del1 - del2) &
              * log((1.0 + del1 * b_v) / (1.0 + del2 * b_v)) &
        ) * (r * t)
    end subroutine

    pure function volume_initalizer(n, p, t) result(v0)
         real(8), intent(in) :: n(:)
         real(8), intent(in) :: p
         real(8), intent(in) :: t
         real(8) :: v0
         v0 = sum(n*b)/sum(b)
    end function
end module
