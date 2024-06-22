module fixtures_models
    use yaeos, only: pr, R, CubicEoS, NRTL
    use yaeos__tapenade_ar_api, only: ArModelTapenade

contains

    type(CubicEoS) function binary_PR76() result(eos)
        use yaeos, only: PengRobinson76
        integer, parameter :: n=2
        real(pr) :: tc(n), pc(n), w(n)
        real(pr) :: kij(n, n), lij(n, n)
        tc = [190._pr, 310._pr]
        pc = [14._pr, 30._pr]
        w = [0.001_pr, 0.03_pr]

        kij = reshape([0., 0.1, 0.1, 0.], [n,n]) 
        lij = kij / 2 
        eos = PengRobinson76(tc, pc, w, kij, lij)
    end function
    
    type(CubicEoS) function binary_PR78() result(eos)
        use yaeos, only: PengRobinson78
        integer, parameter :: n=2
        real(pr) :: tc(n), pc(n), w(n)
        real(pr) :: kij(n, n), lij(n, n)
        tc = [190._pr, 310._pr]
        pc = [14._pr, 30._pr]
        w = [0.001_pr, 0.03_pr]

        kij = reshape([0., 0.1, 0.1, 0.], [n,n]) 
        lij = kij / 2 

        eos = PengRobinson78(tc, pc, w, kij, lij)
    end function
    
    type(CubicEoS) function binary_SRK() result(eos)
        use yaeos, only: SoaveRedlichKwong
        integer, parameter :: n=2
        real(pr) :: tc(n), pc(n), w(n)
        real(pr) :: kij(n, n), lij(n, n)
        tc = [190._pr, 310._pr]
        pc = [14._pr, 30._pr]
        w = [0.001_pr, 0.03_pr]

        kij = reshape([0., 0.1, 0.1, 0.], [n,n]) 
        lij = kij / 2 

        eos = SoaveRedlichKwong(tc, pc, w, kij, lij)
    end function
    
    type(CubicEoS) function binary_RKPR() result(eos)
        use yaeos, only: RKPR
        integer, parameter :: n=2
        real(pr) :: tc(n), pc(n), w(n), zc(n)
        real(pr) :: kij(n, n), lij(n, n)
        tc = [190._pr, 310._pr]
        pc = [14._pr, 30._pr]
        w = [0.001_pr, 0.03_pr]
        zc = [0.23, 0.26]

        kij = reshape([0., 0.1, 0.1, 0.], [n,n]) 
        lij = kij / 2 

        eos = RKPR(tc, pc, w, zc, kij, lij)
    end function

    type(hdPR76) function binary_PR76_hd() result(eos)
        use autodiff_hyperdual_pr76, only: setup, hdPR76
        integer, parameter :: n=2
        real(pr) :: tc(n), pc(n), w(n)
        real(pr) :: kij(n, n), lij(n, n)
        tc = [190._pr, 310._pr]
        pc = [14._pr, 30._pr]
        w = [0.001_pr, 0.03_pr]

        kij = reshape([0., 0.1, 0.1, 0.], [n,n]) 
        lij = kij / 2

        eos = setup(tc, pc, w, kij, lij)
    end function
    
    type(TPR76) function binary_PR76_tape() result(eos)
        use autodiff_tapenade_pr76, only: setup_model, TPR76
        integer, parameter :: n=2
        real(pr) :: tc(n), pc(n), w(n)
        real(pr) :: kij(n, n), lij(n, n)
        tc = [190._pr, 310._pr]
        pc = [14._pr, 30._pr]
        w = [0.001_pr, 0.03_pr]

        kij = reshape([0., 0.1, 0.1, 0.], [n,n]) 
        lij = kij / 2

        eos = setup_model(tc, pc, w, kij, lij)
    end function
end module