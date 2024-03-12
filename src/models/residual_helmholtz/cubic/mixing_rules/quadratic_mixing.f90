module yaeos_models_ar_cubic_quadratic_mixing
    !! Quadratic Mixing Rules for Cubic EoS.
    use yaeos_constants, only: pr
    use yaeos_substance, only: substances
    use yaeos_models_ar_genericcubic, only: CubicMixRule
    implicit none

    type, extends(CubicMixRule) :: QMR
        !! Quadratic Mixing Rule (QMR) derived type. Classic Van der Waals mixing
        !! rules.
        !!
        !! QMR depends on binary interaction parameters, on a Cubic EoS
        !! the mixture is obtained by the combination of an attractive and
        !! repulsive parameter matrices.
        !!
        !! By default the attractive parameter matrix is calculated with:
        !! \[a_{ij} = \sqrt{a_i a_j}(1 - k_{ij})\]
        !! generating the \(a_{ij}\) matrix, but this procedure can be overriden
        !! replacing the `aij` pointer procedure.
        real(pr), allocatable :: k(:, :) !! Attractive Binary Interatction parameter matrix
        real(pr), allocatable :: l(:, :) !! Repulsive Binary Interatction parameter matrix
        procedure(get_aij), pointer :: aij => kij_constant
        !! Procedure to calculate \(a_{ij}\) matrix. Can be overloaded
        !! by any method that respets the interface [[get_aij(interface)]].
    contains
        procedure :: Dmix !! Attractive parameter mixing rule
        procedure :: Bmix !! Repulsive parameter mixing rule
        procedure :: D1mix => D1mix_constant
    end type QMR

    type, extends(QMR) :: QMR_RKPR
    contains
        procedure :: D1mix => RKPR_D1Mix
    end type QMR_RKPR

    abstract interface
        subroutine get_aij(self, ai, daidt, daidt2, aij, daijdt, daijdt2)
            !! Combining rule for the attractive parameter.
            !!
            !! From previously calculated attractive parameters calculate the
            !! \(a_{ij}\) matrix and it's corresponding derivatives.
            import pr, QMR
            class(QMR), intent(in) :: self
            real(pr), intent(in) :: ai(:), daidt(:), daidt2(:)
            real(pr), intent(out):: aij(:, :), daijdt(:, :), daijdt2(:, :)
        end subroutine get_aij
    end interface

contains

    subroutine Dmix(self, n, T, &
        ai, daidt, daidt2, &
        D, dDdT, dDdT2, dDi, dDidT, dDij)
        !! Attractive parameter mixing rule with quadratic mix.
        !!
        !! Takes the all the pure components attractive parameters and their
        !! derivatives with respect to temperature and mix them with the
        !! Van der Waals quadratic mixing rule:
        !!
        !! \[
        !!   D = \sum_i \sum_j n_i n_j a_{ij} = n^2 a_{mix}
        !! \]
        !!
        !! Inside the routine the \(a_{ij}\) matrix is calculated using the
        !! procedure contained in the `QMR` object, this procedures defaults
        !! to the common combining rule: \(a_{ij} = \sqrt{a_i a_j} (1 - k_{ij}) \)
        !!
        !! The procedure can be overloaded by a common one that respects the
        !! interface [[get_aij(interface)]]
        !!
        !! ```fortran
        !! type(QMR) :: my_mixing_rule
        !! my_mixing_rule%aij => new_aij_procedure
        !! ```
        class(QMR), intent(in) :: self !! Mixing rule object.
        real(pr), intent(in) :: T !! Temperature [K]
        real(pr), intent(in) :: n(:) !! Moles vector [mol]
        real(pr), intent(in) :: ai(:) !! Pure components attractive parameters \(a_i\)
        real(pr), intent(in) :: daidt(:) !! \(\frac{da_i}{dT}\)
        real(pr), intent(in) :: daidt2(:) !! \(\frac{d^2a_i}{dT^2}\)

        real(pr), intent(out) :: D !! Mixture attractive parameter \(n^2a_{mix}\)
        real(pr), intent(out) :: dDdT !! \(\frac{dD}{dT}\)
        real(pr), intent(out) :: dDdT2 !! \(\frac{d^2D}{dT^2}\)
        real(pr), intent(out) :: dDi(:) !! \(\frac{dD}{dn_i}\)
        real(pr), intent(out) :: dDidT(:) !! \(\frac{d^2D}{dTn_i}\)
        real(pr), intent(out) :: dDij(:, :)!! \(\frac{d^2D}{dn_{ij}}\)

        integer :: i, j, nc
        real(pr) :: aux, aux2
        real(pr) :: aij(size(ai), size(ai))
        real(pr) :: daijdt(size(ai), size(ai))
        real(pr) :: daijdt2(size(ai), size(ai))

        nc = size(ai)

        if (associated(self%aij)) then
            call self%aij(ai, daidt, daidt2, aij, daijdt, daijdt2)
        else
            write(*, *) "ERROR: aij matrix calculation not defined"
            call exit(1)
        end if

        D = 0
        dDdT = 0
        dDdT2 = 0
        do i = 1, nc
            aux = 0
            aux2 = 0
            dDi(i) = 0
            dDidT(i) = 0

            do j = 1, nc
                dDi(i) = dDi(i) + 2*n(j)*aij(i, j)

                dDidT(i) = dDidT(i) + 2*n(j)*daijdT(i, j)
                aux2 = aux2 + n(j)*daijdT2(i, j)

                dDij(i, j) = 2*aij(i, j)
                aux = aux + n(j)*aij(i, j)
            end do

            D = D + n(i)*aux

            dDdT = dDdT + n(i)*dDidT(i) * 0.5_pr
            dDdT2 = dDdT2 + n(i)*aux2
        end do
    end subroutine Dmix

    subroutine Bmix(self, n, bi, B, dBi, dBij)
        !! Mixture repulsive parameter.
        !!
        !! Calculate the mixture's repulsive parameter and it's derivatives
        !! with respect to composition:
        !!
        !! \[
        !!    B = \sum_i \sum_j n_i n_j \frac{b_i + b_j}{2} (1 - l_{ij})
        !! \]
        !!
        class(QMR), intent(in) :: self !! Mixing rule object.
        real(pr), intent(in) :: n(:) !! Moles vector.
        real(pr), intent(in) :: bi(:) !! Pure components repulsive parameters.
        real(pr), intent(out) :: B !! Mixture repulsive parameter.
        real(pr), intent(out) :: dBi(:) !! \(\frac{dB}{dn_i}\)
        real(pr), intent(out) :: dBij(:, :) !!\(\frac{d^2B}{dn_{ij}}\)

        real(pr) :: bij(size(n), size(n))

        real(pr) :: totn, aux(size(n))

        integer :: i, j, nc

        nc = size(n)
        TOTN = sum(n)
        B = 0.0_pr
        aux = 0.0_pr

        do i = 1, nc
            do j = 1, nc
                bij(i, j) = 0.5_pr * (bi(i) + bi(j)) * (1.0_pr - self%l(i,j))
                aux(i) = aux(i) + n(j) * bij(i, j)
            end do
            B = B + n(i)*aux(i)
        end do

        B = B/totn

        do i = 1, nc
            dBi(i) = (2*aux(i) - B)/totn
            do j = 1, i
                dBij(i, j) = (2*bij(i, j) - dBi(i) - dBi(j))/totn
                dBij(j, i) = dBij(i, j)
            end do
        end do
    end subroutine Bmix

    subroutine D1mix_constant(self, n, d1i, D1, dD1i, dD1ij)
        class(QMR), intent(in) :: self
        real(pr), intent(in) :: n(:)
        real(pr), intent(in) :: d1i(:)
        real(pr), intent(out) :: D1
        real(pr), intent(out) :: dD1i(:)
        real(pr), intent(out) :: dD1ij(:, :)

        D1 = d1i(1)
        dD1i = 0
        dD1ij = 0
    end subroutine D1mix_constant

    subroutine RKPR_D1mix(self, n, d1i, D1, dD1i, dD1ij)
        !! RKPR \(\delta_1\) parameter mixing rule.
        class(QMR_RKPR), intent(in) :: self
        real(pr), intent(in) :: n(:)
        real(pr), intent(in) :: d1i(:)
        real(pr), intent(out) :: D1
        real(pr), intent(out) :: dD1i(:)
        real(pr), intent(out) :: dD1ij(:, :)

        integer :: i, j, nc
        real(pr) :: totn

        nc = size(n)
        totn = sum(n)

        D1 = 0.0_pr
        do i = 1, nc
            D1 = D1 + n(i) * d1i(i)
        end do

        D1 = D1/totn

        do i = 1, nc
            dD1i(i) = (d1i(i) - D1)/totn
            do j = 1, nc
                dD1ij(i, j) = (2.0_pr*D1 - d1i(i) - d1i(j))/totn**2
            end do
        end do
    end subroutine RKPR_D1mix

    subroutine kij_constant(&
        self, a, dadt, dadt2, &
        aij, daijdt, daijdt2 &
        )
        !! Combining rule that uses constant \(k_{ij}\) values.
        class(QMR), intent(in) :: self
        real(pr), intent(in) :: a(:) !! Pure components attractive parameters (\a\)
        real(pr), intent(in) :: dadt(:) !! \(\frac{da}{dT}\)
        real(pr), intent(in) :: dadt2(:) !! \(\frac{d^2a}{dT^2}\)
        real(pr), intent(out) :: aij(:, :) !! \(a_{ij}\) Matrix
        real(pr), intent(out) :: daijdt(:, :) !! \(\frac{da_{ij}}{dT}\)
        real(pr), intent(out) :: daijdt2(:, :)!! \(\frac{d^2a_{ij}}{dT^2}\)

        integer :: i, j
        real(pr) :: sqrt_aii_ajj
        real(pr) :: ai2(size(a)), inner_sum, inner_sum_2
        real(pr) :: aij_daidt

        ai2 = a * a

        do i=1, size(a)
            aij(i, i) = a(i)
            daijdt(i, i) = dadt(i)
            daijdt2(i, i) = dadt2(i)

            do j=i+1, size(a)
                sqrt_aii_ajj = sqrt(a(i) * a(j))

                aij(i, j)    = sqrt_aii_ajj * (1 - self%k(i, j))

                inner_sum = a(i) * dadt(j) + a(j) * dadt(i)
                aij_daidt = aij(i, j) * (0.5_pr * inner_sum)

                daijdt(i, j) = 0.5_pr * aij(i, j) * (inner_sum) / (a(i)*a(j))
                daijdT2(i, j) = &
                    (1 - self%k(i, j))*(dadT(j)*dadT(i)/sqrt(a(i)*a(j)) &
                    + sqrt(a(i)/a(j))*(dadT2(j) - dadT(j)**2/(2*a(j))) &
                    + sqrt(a(j)/a(i))*(dadT2(i) - dadT(i)**2/(2*a(i))))/2

                aij(j, i) = aij(i, j)
                daijdt(j, i) = daijdt(i, j)
                daijdt2(j, i) = daijdt2(i, j)
            end do
        end do
    end subroutine kij_constant
end module yaeos_models_ar_cubic_quadratic_mixing
