program main
    use yaeos, only: pr, R, Substances, AlphaSoave, CubicEoS, GenericCubic_Ar, fugacity_vt, fugacity_tp, vcalc, QMR
    use yaeos, only: ArModel, PengRobinson76
    implicit none

    type(Substances) :: compos
    type(CubicEoS), target :: eos, eos2
    type(AlphaSoave) :: alpha
    type(QMR) :: mixrule

    class(ArModel), pointer :: this

    integer, parameter :: n=2
    real(pr) :: z(n), b, dbi(n), dbij(n,n)
    real(pr) :: v=1.0, t=150.0, p

    real(pr) :: ar
    real(pr) :: art, arv, arv2, art2, artv
    real(pr) :: arn(n), arvn(n), artn(n), arn2(n,n) 
    real(pr) :: lnfug(n), dlnphidp(n), dlnphidt(n), dlnphidn(n, n)

    class(ArModel), allocatable :: models(:)

    real(pr) :: tc(n), pc(n), w(n), kij(n, n), lij(n, n)

    integer :: i

    z = [0.3_pr, 0.7_pr]
    tc = [190._pr, 310._pr]
    pc = [14._pr, 30._pr]
    w = [0.001_pr, 0.03_pr]

    kij = reshape([0., 0.1, 0.1, 0.], [n,n]) 
    lij = kij / 2 

    eos = PengRobinson76(tc, pc, w, kij, lij)
    
    v = 1
    t = 150

    call fugacity_vt(eos, &
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

    print *, "ArVn: ", ArVn
    print *, "ArTn: ", ArTn

    print *, "Arn2: ", Arn2
end program