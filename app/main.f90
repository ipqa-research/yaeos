program main
    use yaeos, only: pr, R, Substances, AlphaSoave, CubicEoS, GenericCubic_Ar
    use yaeos, only: ArModel, PengRobinson78
    implicit none

    type(CubicEoS), target :: eos

    integer, parameter :: n=2
    real(pr) :: z(n)
    real(pr) :: v=1.0, t=150.0, p

    real(pr) :: ar
    real(pr) :: art, arv, arv2, art2, artv
    real(pr) :: arn(n), arvn(n), artn(n), arn2(n,n) 
    real(pr) :: lnfug(n), dlnphidp(n), dlnphidt(n), dlnphidn(n, n)


    real(pr) :: tc(n), pc(n), w(n), kij(n, n), lij(n, n)

    z = [0.3, 0.7]
    Tc = [304.2, 425.1]
    Pc = [73.8, 38.0]
    w = [0.2236, 0.200164]


    eos = PengRobinson78(Tc, Pc, w)

    V = 2.6
    T = 250
    
    call eos%lnphi_vt(&
         z, V, T, P, lnfug, dlnPhidP, dlnphidT, dlnPhidn &
    )

    print *, lnfug

    print *, dlnphidn(1, :)
    print *, dlnphidn(2, :)



end program
