module models
    use yaeos, only: pr, R, CubicEoS

contains

    type(CubicEoS) function binary_PR76() result(eos)
        use yaeos, only: PengRobinson76
        integer, parameter :: n=2
        real(pr) :: tc(n), pc(n), w(n)
        real(pr) :: kij(n, n), lij(n, n)
        tc = [190._pr, 310._pr]
        pc = [14._pr, 30._pr]
        w = [0.001_pr, 0.03_pr]

        kij = reshape([0., 0.0, 0.0, 0.], [n,n]) 
        lij = kij / 2 

        eos = PengRobinson76(tc, pc, w, kij, lij)
    end function

end module