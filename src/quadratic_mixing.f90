module yaeos__models_ar_cubic_quadratic_mixing
   !! Quadratic Mixing Rules for Cubic EoS.
   use yaeos__constants, only: pr, solving_volume
   use yaeos__substance, only: substances
   use yaeos__models_ar_genericcubic, only: CubicMixRule
   use yaeos__models_ar_cubic_mixing_base, only: bmix_qmr
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
   contains
      procedure :: aij => kij_constant !! Default attractive parameter combining rule
      procedure :: Dmix !! Attractive parameter mixing rule
      procedure :: Bmix !! Repulsive parameter mixing rule
      procedure :: D1mix => RKPR_D1mix
   end type QMR

   type, extends(QMR) :: QMRTD
      real(pr), allocatable :: k0(:, :)
      real(pr), allocatable :: Tref(:, :)
   contains
      procedure :: aij => kij_exp_tdep
   end type QMRTD

   abstract interface
      subroutine get_aij(self, T, ai, daidt, daidt2, aij, daijdt, daijdt2)
         !! Combining rule for the attractive parameter.
         !!
         !! From previously calculated attractive parameters calculate the
         !! \(a_{ij}\) matrix and it's corresponding derivatives.
         import pr, QMR
         class(QMR), intent(in) :: self
         real(pr), intent(in) :: T
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

     call self%aij(T, ai, daidt, daidt2, aij, daijdt, daijdt2)

      D = 0
      dDdT = 0
      dDdT2 = 0

      if (solving_volume) then
         D = dot_product(n, matmul(n, aij))
      else
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
      end if
   end subroutine Dmix

   subroutine Bmix(self, n, bi, B, dBi, dBij)
      !! Mixture repulsive parameter.
      !!
      !! Calculate the mixture's repulsive parameter and it's derivatives
      !! with respect to composition:
      !!
      !! \[
      !!    nB = \sum_i \sum_j n_i n_j \frac{b_i + b_j}{2} (1 - l_{ij})
      !! \]
      !!
      class(QMR), intent(in) :: self !! Mixing rule object.
      real(pr), intent(in) :: n(:) !! Moles vector.
      real(pr), intent(in) :: bi(:) !! Pure components repulsive parameters.
      real(pr), intent(out) :: B !! Mixture repulsive parameter.
      real(pr), intent(out) :: dBi(:) !! \(\frac{dB}{dn_i}\)
      real(pr), intent(out) :: dBij(:, :) !!\(\frac{d^2B}{dn_{ij}}\)
      call bmix_qmr(n, bi, self%l, b, dbi, dbij)
   end subroutine Bmix

   subroutine RKPR_D1mix(self, n, d1i, D1, dD1i, dD1ij)
      use yaeos__models_ar_cubic_mixing_base, only: d1mix_rkpr
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
      class(QMR), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: d1i(:)
      real(pr), intent(out) :: D1
      real(pr), intent(out) :: dD1i(:)
      real(pr), intent(out) :: dD1ij(:, :)
      call d1mix_rkpr(n, d1i, d1, dd1i, dd1ij)
   end subroutine RKPR_D1mix

   subroutine kij_constant(&
      self, T, a, dadt, dadt2, &
      aij, daijdt, daijdt2 &
      )
      !! Combining rule that uses constant \(k_{ij}\) values.
      !!
      !! \[
      !!  a_{ij} = \sqrt{a_i a_j} (1 - k_{ij})
      !! ]
      class(QMR), intent(in) :: self
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: a(:) !! Pure components attractive parameters (\a_i\)
      real(pr), intent(in) :: dadt(:) !! \(\frac{da_i}{dT}\)
      real(pr), intent(in) :: dadt2(:) !! \(\frac{d^2a_i}{dT^2}\)
      real(pr), intent(out) :: aij(:, :) !! \(a_{ij}\) Matrix
      real(pr), intent(out) :: daijdt(:, :) !! \(\frac{da_{ij}{dT}\)
      real(pr), intent(out) :: daijdt2(:, :)!! \(\frac{d^2a_{ij}{dT^2}\)

      integer :: i, j
      real(pr) :: sqrt_aii_ajj
      real(pr) :: inner_sum

      do i=1, size(a)
         aij(i, i) = a(i)
         daijdt(i, i) = dadt(i)
         daijdt2(i, i) = dadt2(i)

         do j=i+1, size(a)
            sqrt_aii_ajj = sqrt(a(i) * a(j))

            aij(i, j)    = sqrt_aii_ajj * (1 - self%k(i, j))

            inner_sum = a(i) * dadt(j) + a(j) * dadt(i)
            daijdt(i, j) = 0.5_pr * aij(i, j) * (inner_sum) / (a(i)*a(j))

            daijdt2(i, j) = &
               (1 - self%k(i, j))*(dadT(j)*dadT(i)/sqrt(a(i)*a(j)) &
               + sqrt(a(i)/a(j))*(dadT2(j) - dadT(j)**2/(2*a(j))) &
               + sqrt(a(j)/a(i))*(dadT2(i) - dadT(i)**2/(2*a(i))))/2

            aij(j, i) = aij(i, j)
            daijdt(j, i) = daijdt(i, j)
            daijdt2(j, i) = daijdt2(i, j)
         end do
      end do
   end subroutine kij_constant

   subroutine kij_exp_tdep(&
      self, T, a, dadt, dadt2, &
      aij, daijdt, daijdt2 &
      )
      !! # kij_exp_tdep
      !!
      !! Combining rule that uses temperature dependant \(k_{ij}\) values.
      !! With the following expression:
      !! \[
      !! k_{ij}(T) = k_{ij}^0 + k_{ij}^\infty \exp\left(\frac{-T}{T^*}\right)
      !!  \]
      !!
      !! \[
      !!  a_{ij} = \sqrt{a_i a_j} (1 - k_{ij})
      !! ]
      use hyperdual_mod
      class(QMRTD), intent(in) :: self
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: a(:) !! Pure components attractive parameters (\a_i\)
      real(pr), intent(in) :: dadt(:) !! \(\frac{da_i}{dT}\)
      real(pr), intent(in) :: dadt2(:) !! \(\frac{d^2a_i}{dT^2}\)
      real(pr), intent(out) :: aij(:, :) !! \(a_{ij}\) Matrix
      real(pr), intent(out) :: daijdt(:, :) !! \(\frac{da_{ij}{dT}\)
      real(pr), intent(out) :: daijdt2(:, :)!! \(\frac{d^2a_{ij}{dT^2}\)

      real(pr) :: k0(size(a), size(a))
      real(pr) :: kinf(size(a), size(a))
      real(pr) :: Tstar(size(a), size(a))

      type(hyperdual) :: aij_hd(size(a), size(a)), kij_hd(size(a), size(a)), T_hd
      type(hyperdual) :: a_hd(size(a))

      integer :: i, j, nc

      T_hd = T
      T_hd%f1 = 1
      T_hd%f2 = 1

      k0 = self%k0
      kinf = self%k

      Tstar = self%Tref

      kij_hd = kinf + k0 * exp(-T_hd / Tstar)

      a_hd = a
      
      ! Inject the already calculated derivatives
      a_hd%f1 = dadt
      a_hd%f2 = dadt
      a_hd%f12 = dadt2

      nc = size(a)

      do i=1,size(a)
         aij_hd(i, i) = sqrt(a_hd(i) * a_hd(i))
         do j=i+1,size(a)
            aij_hd(i, j) = sqrt(a_hd(i) * a_hd(j)) * (1._pr - kij_hd(i, j))
            aij_hd(j, i) = aij_hd(i, j)
         end do
      end do

      aij = aij_hd%f0
      daijdt = aij_hd%f1
      daijdt2 = aij_hd%f12
   end subroutine kij_exp_tdep
end module yaeos__models_ar_cubic_quadratic_mixing
