module yaeos__models_ar_cubic_mixing_base
    !! # Mixing rules core math
    !! Procedures of the core calculations of CubicEoS mixing rules.
    !!
    !! # Description
    !! This module holds all the basic math to use mixing rules in other codes.
    !! Keeping it simple and accesible.
    !!
    !! # Examples
    !!
    !! ```fortran
    !! bi = [0.2, 0.3]
    !! lij = reshape([0.0, 0.2, 0.2, 0], [2,2])
    !!
    !! ! Calculate B parameter with Quadratric Mixing Rules.
    !! call bmix_qmr(n, bi, lij, b, dbi, dbij)
    !! 
    !! ```
    !!
    !! # References
    use yaeos__constants, only: pr
    implicit none
contains

    pure subroutine bmix_linear(n, bi, b, dbi, dbij)
        real(pr), intent(in) :: n(:)
        real(pr), intent(in) :: bi(:)
        real(pr), intent(out) :: b, dbi(:), dbij(:, :)

        b = sum(n*bi)
        dbi = bi
        dbij = 0
    end subroutine

    pure subroutine bmix_qmr(n, bi, lij, b, dbi, dbij)
        real(pr), intent(in) :: n(:)
        real(pr), intent(in) :: bi(:)
        real(pr), intent(in) :: lij(:, :)
        real(pr), intent(out) :: b, dbi(:), dbij(:, :)
        
        real(pr) :: bij(size(n), size(n))

        real(pr) :: totn, aux(size(n))

        integer :: i, j, nc

        nc = size(n)
        TOTN = sum(n)
        B = 0
        dBi = 0
        dBij = 0
        aux = 0

        do i = 1, nc
            do j = 1, nc
                bij(i, j) = 0.5_pr * (bi(i) + bi(j)) * (1.0_pr - lij(i,j))
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
    end subroutine

    pure subroutine d1mix_rkpr(n, d1i, d1, dd1i, dd1ij)
        !! RKPR \(\delta_1\) parameter mixing rule.
        !!
        !! The RKPR EoS doesn't have a constant \(\delta_1\) value for each 
        !! component, so a proper mixing rule should be provided. A linear
        !! combination is used.
        !!
        !! \[
        !!     \Delta_1 = \sum_i^N n_i \delta_{1i}
        !! \]
        !!
        real(pr), intent(in) :: n(:)
        real(pr), intent(in) :: d1i(:)
        real(pr), intent(out) :: D1
        real(pr), intent(out) :: dD1i(:)
        real(pr), intent(out) :: dD1ij(:, :)

        integer :: i, j, nc
        real(pr) :: totn

        nc = size(n)
        totn = sum(n)

        D1 = sum(n * d1i)/totn

        do i = 1, nc
            dD1i(i) = (d1i(i) - D1)/totn
            do j = 1, nc
                dD1ij(i, j) = (2 * D1 - d1i(i) - d1i(j))/totn**2
            end do
        end do
    end subroutine

    subroutine lamdba_hv(d1, dd1i, dd1ij, L, dLi, dLij)
        !! Infinite pressure limit parameter \(\Lambda\)
        !!
        !! \[
        !! \Lambda = \frac{1}{\delta_1 + \delta_2} \ln \frac{1 + \delta_1}{1 + \delta_2}
        !! \]
        real(pr), intent(in) :: d1
        real(pr), intent(in) :: dd1i(:)
        real(pr), intent(in) :: dd1ij(:, :)
        real(pr), intent(out) :: L
        real(pr), intent(out) :: dLi(:)
        real(pr), intent(out) :: dLij(:, :)

        real(pr) :: f, g, h
        real(pr), dimension(size(dd1i)) :: df, dg, dh
        real(pr), dimension(size(dd1i), size(dd1i)) :: d2f, d2g, d2h

        integer :: i, j, nc

        nc = size(dd1i)

        f = d1 + 1
        g = (d1 + 1)*d1 + d1 - 1
        h = log((d1+1)**2 / 2)

        L = f/g * h

        df = dd1i
        dg = 2*(d1 + 1)*dd1i
        dh = 2 * dd1i/(d1 + 1)

        dLi = f/g * dh - f*h*dg/g**2 + h * df/g
        
        do concurrent (i=1:nc, j=1:nc)
            d2f(i, j) = dd1ij(i, j)
            d2g(i, j) = 2*dd1ij(i, j)*(d1 + 1) + 2*dd1i(i)*dd1i(j)
            d2h(i, j) = 2*(dd1ij(i, j)/(d1 + 1) - dd1i(i)*dd1i(j)/(d1 + 1)**2)
        end do

        dLij =  0
        do concurrent (i=1:nc, j=1:nc)
            dLij(i, j) = f/g * d2h(i, j) - f * dh(i) * dg(j) / g**2 + dh(i)*df(j)/g
            dLij(i, j) = dLij(i, j) - f * h * d2g(i, j)/g**2 + f*dg(i)*dh(j)/g**2 &
            - 2*f*h*dg(i)*dg(j)/g**3 + h*df(j)*dg(i)/g**2
            dLij(i, j) = dLij(i, j) + h * d2f(i,j)/g + df(i)*dh(j)/g - h*df(i)*dg(j)/g**2
        end do
    end subroutine
end module