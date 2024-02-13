module hyperdual_pr76
    use yaeos_constants, only: pr, R
    use yaeos_ar_models_hyperdual
    use yaeos_substance, only: Substances
    implicit none

    type, extends(ArModelAdiff) :: PR76
        type(Substances) :: composition
        real(pr), allocatable :: kij(:, :), lij(:, :)
        real(pr), allocatable :: ac(:), b(:), k(:)
        real(pr), allocatable :: tc(:), pc(:), w(:)
    contains
        procedure :: Ar => arfun
        procedure :: get_v0 => v0
    end type

    real(pr), parameter :: del1 = 1._pr + sqrt(2._pr)
    real(pr), parameter :: del2 = 1._pr - sqrt(2._pr)

contains

    type(PR76) function setup(tc_in, pc_in, w_in, kij_in, lij_in) result(self)
        !-| Setup the enviroment to use the PengRobinson 76 Equation of State
        !   It uses the Cubic Van der Waals mixing rules
        real(pr) :: tc_in(:)
        real(pr) :: pc_in(:)
        real(pr) :: w_in(:)
        real(pr) :: kij_in(:, :)
        real(pr) :: lij_in(:, :)

        self%tc = tc_in
        self%pc = pc_in
        self%w = w_in

        self%ac = 0.45723553_pr * R**2 * self%tc**2 / self%pc
        self%b = 0.07779607_pr * R * self%tc/self%pc
        self%k = 0.37464_pr + 1.54226_pr * self%w - 0.26993_pr * self%w**2

        self%kij = kij_in
        self%lij = lij_in
    end function

    function arfun(self, n, v, t) result(ar)
        class(PR76) :: self
        type(hyperdual), intent(in) :: n(:), v, t
        type(hyperdual) :: ar
    
        type(hyperdual) :: amix, a(size(n)), ai(size(n)), z2(size(n)), zij
        type(hyperdual) :: bmix
        type(hyperdual) :: b_v, nij

        integer :: i, j

        associate(pc => self%pc, ac => self%ac, b => self%b, k => self%k, &
                  kij => self%kij, lij => self%kij, tc => self%tc &
        )
            ! a = ac * (1.0_pr + k * (1.0_pr - sqrt(t/self%tc)))**2
            
            ! amix = 0.0_pr
            ! bmix = 0.0_pr
            ! do i=1,size(n)
            !     do j=1,size(n)
            !         amix = amix + n(i) * n(j) * sqrt(a(i) * a(j)) * (1._pr - kij(i, j))
            !         bmix = bmix + n(i) * n(j) * 0.5_pr * (b(i) + b(j)) * (1._pr - lij(i, j))
            !     end do
            ! end do

            ! bmix = bmix/sum(n)

            ! b_v = bmix/v
            
            ! ar = (&
            !     - sum(n) * log(1.0_pr - b_v) &
            !     - amix / (R*t*bmix)*1.0_pr / (del1 - del2) &
            !     * log((1.0_pr + del1 * b_v) / (1.0_pr + del2 * b_v)) &
            ! ) * R*t
            a = sqrt(ac * (1.0_pr + k * (1.0_pr - sqrt(t/tc)))**2)
        
            amix = 0.0_pr
            bmix = 0.0_pr

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
            ar = (&
                - sum(n) * log(1.0_pr - b_v) &
                - amix / (R*t*bmix)*1.0_pr / (del1 - del2) &
                * log((1.0_pr + del1 * b_v) / (1.0_pr + del2 * b_v)) &
            ) * (r * t)
            end associate
    end function

    function v0(self, n, p, t)
        !! Initialization of volume
        class(PR76), intent(in) :: self
        real(pr), intent(in) :: n(:)
        real(pr), intent(in) :: p
        real(pr), intent(in) :: t
        real(pr) :: v0

        v0 = sum(n * self%b) / sum(n)
    end function

    subroutine adiff_pr76 
        class(ArModel), allocatable :: eos

        integer, parameter :: n=2
        real(pr) :: tc(n), pc(n), w(n), kij(n, n), lij(n, n)
        real(pr) :: z(n), lnfug(n), dlnphidp(n), dlnphidt(n), dlnphidn(n, n), v, t, p
        real(pr) :: ar
        real(pr) :: art, arv, arv2, art2, artv
        real(pr) :: arn(n), arvn(n), artn(n), arn2(n,n) 

        real(8) :: time, std, mean
        real(8) :: et, st
        integer :: i, nevals=1e3

        z = [0.3_pr, 0.7_pr]
        tc = [190._pr, 310._pr]
        pc = [14._pr, 30._pr]
        w = [0.001_pr, 0.03_pr]

        kij = reshape([0._pr, 0.1_pr, 0.1_pr, 0._pr], [n,n]) 
        lij = kij / 2._pr

        eos = setup(tc, pc, w, kij, lij)
        v = 1
        t = 150

        call cpu_time(st)
        do i=1,nevals
            call eos%residual_helmholtz(&
                    z, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
                    ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
            )
        end do
        call cpu_time(et)
        print *, (et-st)/nevals
        print *, "Ar: ", ar

        print *, "ArV: ", arV
        print *, "ArT: ", arT

        print *, "ArT2: ", arT2
        print *, "ArV2: ", ArV2
        
        print *, "ArTV: ", ArTV

        print *, "ArVn: ", ArVn
        print *, "ArTn: ", ArTn

        print *, "Arn2: ", Arn2
    end subroutine
end module