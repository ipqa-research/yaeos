module yaeos_models_ar_genericcubic_quadratic_mixing
    use yaeos_constants, only: pr
    use yaeos_substance, only: substances
    use yaeos_models_ar_genericcubic, only: CubicMixRule
    implicit none
   
    type, extends(CubicMixRule) :: QMR
        real(pr), allocatable :: k(:, :), l(:, :)
        procedure(get_aij), pointer :: aij => kij_constant
    contains
        procedure :: Dmix
        procedure :: Bmix
    end type
    
    abstract interface
        subroutine get_aij(&
            self, ai, daidt, daidt2, &
            aij, daijdt, daijdt2)
            import pr, QMR
            class(QMR) :: self
            real(pr), intent(in) :: ai(:), daidt(:), daidt2(:)
            real(pr), intent(out):: aij(:, :), daijdt(:, :), daijdt2(:, :)
        end subroutine
    end interface

contains

    subroutine Dmix(self, n, T, &
        ai, daidt, daidt2, &
        D, dDdT, dDdT2, dDi, dDidT, dDij)
        class(QMR), intent(in) :: self
        real(pr), intent(in) :: T, n(:)
        real(pr), intent(in) :: ai(:), daidt(:), daidt2(:)
        real(pr), intent(out) :: D, dDdT, dDdT2, dDi(:), dDidT(:), dDij(:, :)

        integer :: i, j, nc
        real(pr) :: aux, aux2
        real(pr) :: aij(size(ai), size(ai))
        real(pr) :: daijdt(size(ai), size(ai))
        real(pr) :: daijdt2(size(ai), size(ai))

        nc = size(ai)

        if (associated(self%aij)) &
             call self%aij(ai, daidt, daidt2, aij, daijdt, daijdt2)

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
    end subroutine
    
    subroutine Bmix(self, n, bi, B, dBi, dBij)
      !! Mixture repulsive parameter.
      class(QMR), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: bi(:)
      real(pr), intent(out) ::  B, dBi(:), dBij(:, :)

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
    end subroutine

    subroutine kij_constant(&
        self, a, dadt, dadt2, &
        aij, daijdt, daijdt2 &
                            )
      class(QMR) :: self
      real(pr), intent(in) :: a(:), dadt(:), dadt2(:)
      real(pr), intent(out) :: aij(:, :), daijdt(:, :), daijdt2(:, :)

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
            daijdt2(i, j) = &
                  aij_daidt / (ai2(i) * ai2(j)) &
                - aij_daidt * dadt(j) / (a(i) * ai2(j)) &
                - aij_daidt * dadt(i) / (a(j) * ai2(i)) &
                + aij(i, j) * (&
                      0.5_pr * (a(i) * dadt2(j) + a(j) * dadt2(i)) &
                      + dadt(i) * dadt(j)&
                      ) / (a(i) * a(j))
            
            aij(j, i) = aij(i, j)
            daijdt(j, i) = daijdt(i, j)
            daijdt2(j, i) = daijdt2(i, j)
         end do
      end do
    end subroutine
end module