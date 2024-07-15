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

    type(Substance) :: sus(2)

    real(pr) :: v, T, P
    real(pr) :: lnphip(nc)

    forsus_dir = "build/dependencies/forsus/data/json"
    sus(1) = Substance("propane", only=["critical"])
    sus(2) = Substance("n-butane", only=["critical"])

    n = [0.7, 0.3]
    tc = sus%critical%critical_temperature%value
    pc = sus%critical%critical_pressure%value/1e5
    w = sus%critical%acentric_factor%value

    eos = SoaveRedlichKwong(tc, pc, w)
    T = 350
    P = 10

    call eos%volume(n, P=P,T=T, root_type="stable", V=V)
    print *, V, P
    call volume_michelsen(eos, n, P=P,T=T, root_type="stable", V=V)
    print *, V, P


    call eos%pressure(n, V, T, P)
    print *, T, V, P
    call eos%lnphi_tp(n, T, P, V=V, root_type="vapor", lnphip=lnphip)
    print *, T, P, V, lnphip
end program main
