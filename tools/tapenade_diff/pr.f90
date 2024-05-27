module tapenade_pr
    use yaeos_constants, only: pr, R
    use yaeos_tapenade_ar_api, only: ArModelTapenade
    implicit none

    !type, extends(ArModelTapenade) :: TPR76
    type :: TPR76
        real(8), allocatable :: kij(:, :), lij(:, :)
        real(8), allocatable :: ac(:), b(:), k(:)
        real(8), allocatable :: tc(:), pc(:), w(:)
        real(8) :: del1 = 1. + sqrt(2.)
        real(8) :: del2 = 1. - sqrt(2.)
    ! contains
    !  procedure :: ar
    !  procedure :: ar_d
    !  procedure :: ar_b
    !  procedure :: ar_d_b
    !  procedure :: ar_d_d
    !  procedure :: v0
    end type
contains
    type(TPR76) function setup_model(tc, pc, w, kij, lij)
        real(8), intent(in) :: tc(:)
        real(8), intent(in) :: pc(:)
        real(8), intent(in) :: w(:)
        real(8), intent(in) :: kij(:, :)
        real(8), intent(in) :: lij(:, :)

        setup_model%components%tc = tc
        setup_model%components%pc = pc
        setup_model%components%w = w

        setup_model%ac = 0.45723553 * R**2 * tc**2 / pc
        setup_model%b = 0.07779607 * R * tc/pc
        setup_model%k = 0.37464 + 1.54226 * w - 0.26993 * w**2

        setup_model%kij = kij
        setup_model%lij = lij
    end function

    subroutine ar(model, n, v, t, arval)
        type(TPR76), intent(in) :: model
        real(8), intent(in) :: n(:), v, t
        real(8), intent(out) :: arval
        real(8) :: amix, a(size(n)), ai(size(n)), z2(size(n)), nij
        real(8) :: bmix
        real(8) :: b_v
        real(8) :: aij(size(n), size(n)), bij(size(n), size(n))
        integer :: i, j

        real(8) :: tc(size(n)), pc(size(n)), w(size(n))
        real(8) :: ac(size(n)), b(size(n)), k(size(n))
        real(8) :: del1, del2
        real(8) :: kij(size(n), size(n)), lij(size(n), size(n))

        tc = model%components%tc
        pc = model%components%pc
        w = model%components%w
        ac = model%ac
        b = model%b
        k = model%k

        kij = model%kij
        lij = model%lij

        del1 = model%del1
        del2 = model%del2

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
        bmix = bmix + sum(n**2 * model%b)

        bmix = bmix/sum(n)

        b_v = bmix/v
        arval = (&
              - sum(n) * log(1.0 - b_v) &
              - amix / (R*t*bmix)*1.0 / (del1 - del2) &
              * log((1.0 + del1 * b_v) / (1.0 + del2 * b_v)) &
        ) * (r * t)
    end subroutine

    pure function volume_initalizer(model, n, p, t) result(v0)
         type(TPR76), intent(in) :: model
         real(8), intent(in) :: n(:)
         real(8), intent(in) :: p
         real(8), intent(in) :: t
         real(8) :: v0
         v0 = sum(n*model%b)/sum(model%b)
    end function
end module
