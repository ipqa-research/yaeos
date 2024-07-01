program main
    use yaeos, only: pr, ArModel, PengRobinson76, volume, pressure, SoaveRedlichKwong
    use forsus, only: Substance, forsus_dir
    implicit none
    class(ArModel), allocatable :: eos

    integer, parameter :: nc = 2
    real(pr), dimension(nc) :: n, tc, pc, w
    real(pr), dimension(nc, nc) :: kij, lij

    type(Substance) :: sus(2)

    real(pr) :: v
    real(pr) :: p0, pf, dp, p
    real(pr) :: t0, tf, dt, t
    integer :: i, j, n_p_points=500, n_t_points=5

    forsus_dir = "build/dependencies/forsus/data/json"
    sus(1) = Substance("propane", only=["critical"])
    sus(2) = Substance("n-butane", only=["critical"])

    n = [0.7, 0.3]
    tc = sus%critical%critical_temperature%value
    pc = sus%critical%critical_pressure%value/1e5
    w = sus%critical%acentric_factor%value
    kij = reshape([0., 0.0, 0.0, 0.], [nc,nc])
    lij = kij / 2

    eos = SoaveRedlichKwong(tc, pc, w, kij, lij)

    p0 = 100
    pf = 0.1
    dp = (pf-p0)/n_p_points

    t0 = 150
    tf = 350
    dt = (tf-t0)/n_t_points

    T = 300
    do i=1,1000
        V = real(i, pr)/1000._pr
        call pressure(eos, n, V, T, P=P)
        print *, V, P
    end do



    ! do j=0, n_t_points - 1
    !     t = t0 + dt * j
    !     print *, "# ", t
    !     do i=0,n_p_points-1
    !         p = p0 + dp * i
    !         call volume(eos, n, p, t, v, root_type="stable")
    !         print *, v, p
    !     end do
    !     print *, ""
    !     print *, ""
    ! end do
end program main
