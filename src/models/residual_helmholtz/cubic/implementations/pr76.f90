module yaeos_models_ar_cubic_pengrobinson76
    !! PengRobinson 76 implementation.
    use yaeos_constants, only: R, pr
    use yaeos_substance, only: Substances
    use yaeos_models_ar_genericcubic, only: CubicEoS, AlphaFunction, CubicMixRule
    use yaeos_models_ar_genericcubic_quadratic_mixing, only: QMR
    use yaeos_models_ar_cubic_alphas, only: AlphaSoave
    implicit none

    private

    public :: PengRobinson76

contains

    type(CubicEoS) function PengRobinson76(tc, pc, w, kij, lij) result(model)
        !! PengRobinson76.
        !!
        !! Using the critical constants setup the parameters to use the 
        !! PengRobinson Equation of State
        !!
        !! - \[\alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2\]
        !! - \[k = 0.37464 + 1.54226 * \omega - 0.26993 \omega^2 \]
        !! - \[a_c = 0.45723553  R^2 T_c^2 / P_c\]
        !! - \[b = 0.07779607r  R T_c/P_c\]
        !! - \[\delta_1 = 1 + \sqrt{2}\]
        !! - \[\delta_2 = 1 - \sqrt{2}\]
        !!
        !! There is also the optional posibility to include the \(k_{ij}\) and
        !! \(l_{ij}\) matrices. Using by default Classic Van der Waals mixing
        !! rules.
        !!
        !! After setting up the model, it is possible to redefine either the
        !! mixing rule or the alpha function using a different derived type
        !! defined outside the function.
        real(pr), intent(in) :: tc(:) !! Critical Temperatures [K]
        real(pr), intent(in) :: pc(:) !! Critical Pressures [bar]
        real(pr), intent(in) :: w(:) !! Acentric Factors
        real(pr), optional, intent(in) :: kij(:, :) !! \(k_{ij}\) matrix
        real(pr), optional, intent(in) :: lij(:, :) !! \(l_{ij}\) matrix

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