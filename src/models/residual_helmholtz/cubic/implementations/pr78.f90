module yaeos_models_ar_cubic_pengrobinson78
    use yaeos_constants, only: R, pr
    use yaeos_substance, only: Substances
    use yaeos_models_ar_genericcubic, only: CubicEoS, AlphaFunction, CubicMixRule
    use yaeos_models_ar_genericcubic_quadratic_mixing, only: QMR
    use yaeos_models_ar_cubic_alphas, only: AlphaSoave
    implicit none

    private

    public :: PengRobinson78
contains
    type(CubicEoS) function PengRobinson78(tc, pc, w, kij, lij) result(model)
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

        alpha%k = 0.37464_pr &
                 + 1.54226_pr * composition%w &
                 - 0.26993_pr * composition%w**2

        where (composition%w <=0.491)
              alpha%k = 0.37464 + 1.54226 * composition%w(i) - 0.26992 * composition%w**2
        elsewhere
              alpha%k = 0.379642 + 1.48503 * composition%w - 0.164423 * composition%w**2 + 0.016666 * composition%w**3
        end where

        if (present(kij)) then
            mixrule%k = kij
        else
            mixrule%k = reshape([(0, i=1,nc**2)], [nc, nc])
        endif
        
        if (present(lij)) then
            mixrule%l = lij
        else
            mixrule%l = reshape([(0, i=1,nc**2)], [nc, nc])
        endif

        model%components = composition
        model%ac = 0.45723553_pr * R**2 * composition%tc**2 / composition%pc
        model%b = 0.07779607_pr * R * composition%tc/composition%pc
        model%del1 = [1 + sqrt(2.0_pr)]
        model%del2 = [1 - sqrt(2.0_pr)]
        model%alpha = alpha
        model%mixrule = mixrule
    end function
end module
