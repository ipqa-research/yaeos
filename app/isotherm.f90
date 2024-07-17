program main
    use yaeos, only: pr, ArModel, PengRobinson76, SoaveRedlichKwong
    use yaeos__models_ar, only: nvolume => volume
    use yaeos__models_solvers, only: volume_michelsen
    use forsus, only: Substance, forsus_dir
    use fortime, only: timer
    implicit none
    class(ArModel), allocatable :: eos

    integer, parameter :: nc = 2
    real(pr), dimension(nc) :: n, tc, pc, w
    real(pr), dimension(nc, nc) :: kij, lij

    type(timer) :: tim

    type(Substance) :: sus(2)

    real(pr) :: v, T, P
    real(pr) :: lnphip(nc)

    integer :: i, nevals=100000

    forsus_dir = "build/dependencies/forsus/data/json"
    sus(1) = Substance("propane", only=["critical"])
    sus(2) = Substance("n-butane", only=["critical"])

    n = [0.99, 0.01]
    tc = sus%critical%critical_temperature%value
    pc = sus%critical%critical_pressure%value/1e5
    w = sus%critical%acentric_factor%value

    eos = SoaveRedlichKwong(tc, pc, w)
    T = 350
    P = 500

    do i=1,10000
        P = real(i, pr)/10
        call eos%volume(n, P, T, V, "stable")
        print *, V, P
    end do

    print "(/)"

    do i=1,1000
        P = real(i, pr)/10
        call volume_michelsen(eos, n, P, T, V, "stable")
        print *, V, P
    end do
end program main
