module yaeos_models_ar_cubic_srk
    use yaeos_constants, only: R, pr
    use yaeos_substance, only: Substances
    use yaeos_models_ar_genericcubic, only: CubicEoS, AlphaFunction, CubicMixRule
    use yaeos_models_ar_genericcubic_quadratic_mixing, only: QMR
    use yaeos_models_ar_cubic_alphas, only: AlphaSoave
    implicit none

    private

    public :: SoaveRedlichKwong

contains

    type(CubicEoS) function SoaveRedlichKwong(tc, pc, w, kij, lij) result(model)
        real(pr), intent(in) :: tc(:), pc(:), w(:)
        real(pr), optional, intent(in) :: kij(:, :), lij(:, :)

        type(Substances) :: composition
        type(QMR) :: mixrule
        type(AlphaSoave) :: alpha
        integer :: nc
        integer :: i

        nc = size(tc)
        
        composition%tc = tc
        composition%pc = pc
        composition%w = w

        alpha%k = 0.48_pr + 1.574_pr * composition%w - 0.175_pr * composition%w**2

        if (present(kij)) then
            mixrule%k = kij
        else
            mixrule%k = reshape([(i, i=1,nc**2)], [nc, nc])
        endif
        
        if (present(lij)) then
            mixrule%l = lij
        else
            mixrule%l = reshape([(i, i=1,nc**2)], [nc, nc])
        endif

        model%components = composition
        model%ac = 0.427480_pr * R**2 * composition%tc**2/composition%pc
        model%b = 0.086640_pr * R * composition%tc/composition%pc
        model%del1 = [1]
        model%del2 = [0]
        model%alpha = alpha
        model%mixrule = mixrule
    end function
end module
