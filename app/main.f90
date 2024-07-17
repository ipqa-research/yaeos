program main
    use yaeos, only: pr, R, Substances, AlphaSoave, CubicEoS, GenericCubic_Ar
    use yaeos, only: ArModel, PengRobinson76
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

    z = [0.3_pr, 0.7_pr]
    tc = [190._pr, 310._pr]
    pc = [14._pr, 30._pr]
    w = [0.001_pr, 0.03_pr]

    kij = reshape([0., 0.1, 0.1, 0.], [n,n]) 
    lij = kij / 2 

    eos = PengRobinson76(tc, pc, w, kij, lij)
    
    v = 1
    t = 150

    call eos%lnphi_vt(&
         z, V, T, P, lnfug, dlnPhidP, dlnphidT, dlnPhidn &
    )

    call eos%residual_helmholtz(&
            z, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
            ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
    )

    print *, "Ar: ", ar

    print *, "ArV: ", arV
    print *, "ArT: ", arT

    print *, "ArT2: ", arT2
    print *, "ArV2: ", ArV2
    
    print *, "ArTV: ", ArTV
    
    print *, "Arn: ", Arn

    print *, "ArVn: ", ArVn
    print *, "ArTn: ", ArTn

    print *, "Arn2: ", Arn2
end program
