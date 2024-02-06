module yaeos_models_ar_cubic_srk
    !! SoaveRedlichKwong implementation.
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
        !! SoaveRedlichKwong.
        !!
        !! Using the critical constants setup the parameters to use the 
        !! SoaveRedlichKwong Equation of State
        !!
        !! - \[\alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2\]
        !! - \[k = 0.48 + 1.574 \omega - 0.175 \omega^2 \]
        !! - \[a_c = 0.427480  R^2 * T_c^2/P_c\]
        !! - \[b = 0.086640  R T_c/P_c\]
        !! - \[\delta_1 = 1\]
        !! - \[\delta_2 = 0\]
        !!
        !! There is also the optional posibility to include the k_{ij} and l_{ij}
        !! matrices. Using by default Classic Van der Waals mixing rules.
        !!
        !! After setting up the model, it is possible to redefine either the
        !! mixing rule or the alpha function using a different derived type
        !! defined outside the function.
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
            mixrule%k = reshape([(0, i=1,nc**2)], [nc, nc])
        endif
        
        if (present(lij)) then
            mixrule%l = lij
        else
            mixrule%l = reshape([(0, i=1,nc**2)], [nc, nc])
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
