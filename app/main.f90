program main
    use yaeos, only: pr, R, Substances, AlphaSoave, CubicEoS, GenericCubic_Ar, fugacity_vt, fugacity_tp, vcalc
    use yaeos, only: ArModel
    use mixing, only: QMR
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
    real(pr) :: lnfug(n), dlnphidp(n), dlnphidt(n), dlnphidn(n)

    class(ArModel), allocatable :: models(:)

    integer :: i

    z = [0.3_pr, 0.7_pr]
    compos%tc = [190._pr, 310._pr]
    compos%pc = [14._pr, 30._pr]
    compos%w = [0.001_pr, 0.03_pr]

    alpha%k = 0.37464_pr + 1.54226_pr * compos%w - 0.26993_pr * compos%w**2

    mixrule%k = reshape([0., 0.1, 0.1, 0.], [n,n]) * 0
    mixrule%l = mixrule%k / 2 * 0

    eos = set_eos(compos, alpha, mixrule)
    eos2 = set_eos(compos, alpha, mixrule)

    models = [eos, eos2]
    this => eos
    ! do i=1,2
    !     print *, loc(models(i))
    !     call models(i)%residual_helmholtz(&
    !         z, v, t, ar, arv, art, artv, arv2, art2, arn, arvn, artn, arn2&
    !     )
    ! end do

    call fugacity_vt(eos, &
         z, V, T, P, lnfug, dlnPhidP, dlnphidT, dlnPhidn &
    )

    p = 1.0
    call fugacity_tp(eos, &
         z, T, P, V, 1, lnfug, dlnPhidP, dlnphidT, dlnPhidn &
    )

contains
    type(CubicEoS) function set_eos(composition, alphafun, mixrule)
        use yaeos_substance, only: Substances
        use yaeos_models_ar_genericcubic, only: CubicEoS, AlphaFunction, CubicMixRule
        
        class(Substances) :: composition
        
        class(AlphaFunction), optional :: alphafun
        class(CubicMixRule), optional :: mixrule

        type(AlphaSoave) :: default_alphafun
        type(QMR) :: default_mixrule

        real(pr) :: k(size(composition%w))

        integer :: i, nc

        set_eos%ac = 0.45723553_pr * R**2 * composition%tc**2 / composition%pc
        set_eos%b = 0.07779607_pr * R * composition%tc/composition%pc
        set_eos%del1 = [1 + sqrt(2.0_pr)]
        set_eos%del2 = [1 - sqrt(2.0_pr)]
        
        set_eos%components = composition

        if (present(alphafun)) then
            set_eos%alpha = alphafun
        else
            default_alphafun%k = 0.37464_pr + 1.54226_pr * compos%w - 0.26993_pr * compos%w**2
            set_eos%alpha = default_alphafun
        end if

        if (present(mixrule)) then
            set_eos%mixrule = mixrule
        else
            default_mixrule%k=reshape([(0._pr, i=1,n**2)],[nc,nc])
            default_mixrule%l=reshape([(0._pr, i=1,n**2)],[nc,nc])
            set_eos%mixrule = default_mixrule
        end if
    end function
end program