module autodiff_hyperdual_pr76
    use yaeos__constants, only: pr, R
    use yaeos__adiff_hyperdual_ar_api
    implicit none

    type, extends(ArModelAdiff) :: hdPR76
        real(pr), allocatable :: kij(:, :), lij(:, :)
        real(pr), allocatable :: ac(:), b(:), k(:)
    contains
        procedure :: Ar => arfun
        procedure :: get_v0 => v0
    end type

    real(pr), private, parameter :: del1 = 1._pr + sqrt(2._pr)
    real(pr), private, parameter :: del2 = 1._pr - sqrt(2._pr)

contains

    type(hdPR76) function setup(tc, pc, w, kij, lij) result(self)
        !! Seup an Autodiff_PR76 model
        real(pr) :: tc(:)
        real(pr) :: pc(:)
        real(pr) :: w(:)
        real(pr) :: kij(:, :)
        real(pr) :: lij(:, :)

        self%components%tc = tc
        self%components%pc = pc
        self%components%w = w

        self%ac = 0.45723553_pr * R**2 * tc**2 / pc
        self%b = 0.07779607_pr * R * tc/pc
        self%k = 0.37464_pr + 1.54226_pr * w - 0.26993_pr * w**2

        self%kij = kij
        self%lij = lij
    end function

    function arfun(self, n, v, t) result(ar)
        class(hdPR76) :: self
        type(hyperdual), intent(in) :: n(:), v, t
        type(hyperdual) :: ar
    
        type(hyperdual) :: amix, a(size(n)), ai(size(n))
        type(hyperdual) :: bmix
        type(hyperdual) :: b_v, nij

        integer :: i, j

        associate(pc => self%components%pc, ac => self%ac, b => self%b, k => self%k, &
                  kij => self%kij, lij => self%lij, tc => self%components%tc &
        )
            a = 1.0_pr + k * (1.0_pr - sqrt(t/tc))
            a = ac * a ** 2
            ai = sqrt(a)

        
            amix = 0.0_pr
            bmix = 0.0_pr

            do i=1,size(n)-1
                do j=i+1,size(n)
                    nij = n(i) * n(j)
                    amix = amix + 2 * nij * (ai(i) * ai(j)) * (1 - kij(i, j))
                    bmix = bmix + nij * (b(i) + b(j)) * (1 - lij(i, j))
                end do
            end do

            amix = amix + sum(n**2*a)
            bmix = bmix + sum(n**2 * b)

            bmix = bmix/sum(n)

            b_v = bmix/v
            ar = (&
                - sum(n) * log(1.0_pr - b_v) &
                - amix / (R*t*bmix)*1.0_pr / (del1 - del2) &
                * log((1.0_pr + del1 * b_v) / (1.0_pr + del2 * b_v)) &
            ) * (r * t)
            end associate
    end function

    function v0(self, n, p, t)
        !! Initialization of volume
        class(hdPR76), intent(in) :: self
        real(pr), intent(in) :: n(:)
        real(pr), intent(in) :: p
        real(pr), intent(in) :: t
        real(pr) :: v0

        v0 = sum(n * self%b) / sum(n)
    end function
end module
