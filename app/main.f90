

program main
    use yaeos_constants, only: pr, R
    use yaeos_substance, only: substances
    use yaeos_models_ar_genericcubic, only: CubicEoS, GenericCubic_Ar

    use mixing, only: QMR
    use alpha_soave, only: AlphaSoave
    implicit none

    type(Substances) :: compos
    type(CubicEoS) :: eos
    type(AlphaSoave) :: alpha
    type(QMR) :: mixrule

    integer, parameter :: n=2
    real(pr) :: z(n), b, dbi(n), dbij(n,n)
    real(pr) :: v=1.0, t=150.0

    real(pr) :: ar
    real(pr) :: art, arv, arv2, art2, artv
    real(pr) :: arn(n), arvn(n), artn(n), arn2(n,n) 

    integer :: i

    z = [0.3_pr, 0.7_pr]
    compos%tc = [190._pr, 310._pr]
    compos%pc = [14._pr, 30._pr]
    compos%w = [0.001_pr, 0.03_pr]

    alpha%k = 0.37464_pr + 1.54226_pr * compos%w - 0.26993_pr * compos%w**2

    mixrule%k = reshape([0., 0.1, 0.1, 0.], [n,n])
    mixrule%l = mixrule%k / 2

    eos = set_eos(compos, alpha, mixrule)

    do i=1,1e6
        call GenericCubic_Ar(&
        eos, z, v, t, ar, arv, art, artv, arv2, art2, arn, arvn, artn, arn2)
    end do

    print *, ar
    print *, arv
    print *, art
    print *, artv
    print *, arv2
    print *, art2
    print *, arn
    print *, arvn
    print *, artn
    print *, arn2
contains
    type(CubicEoS) function set_eos(composition, alphafun, mixrule)
        use yaeos_substance, only: Substances
        use yaeos_models_ar_genericcubic, only: CubicEoS, AlphaFunction, CubicMixRule
        
        class(Substances) :: composition
        class(AlphaFunction) :: alphafun
        class(CubicMixRule) :: mixrule

        set_eos%ac = 0.45723553_pr * R**2 * composition%tc**2 / composition%pc
        set_eos%b = 0.07779607_pr * R * composition%tc/composition%pc

        set_eos%del1 = [1 + sqrt(2.0_pr)]
        set_eos%del2 = [1 - sqrt(2.0_pr)]

        set_eos%components = composition
        set_eos%alpha = alphafun
        set_eos%mixrule = mixrule
    end function
end program