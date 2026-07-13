module yaeos__models_ar_cubic_quadratic_mixing
   !! Cubic Mixing Rules for Cubic EoS.
   use yaeos__constants, only: pr, solving_volume
   use yaeos__substance, only: substances
   use yaeos__models_ar_genericcubic, only: CubicMixRule
   use yaeos__models_ar_cubic_mixing_base, only: bmix_qmr
   implicit none

   type, extends(CubicMixRule) :: CMR
      !! Cubic Mixing Rule (CMR) derived type. Classic Van der Waals mixing
      !! rules.
      !!
      !! QMR depends on binary interaction parameters, on a Cubic EoS
      !! the mixture is obtained by the combination of an attractive and
      !! repulsive parameter matrices.
      !!
      !! By default the attractive parameter matrix is calculated with:
      !! \[a_{ijk} = \sqrt[3]{a_i a_j a_k}(1 - k_{ijk})\]
      !! generating the \(a_{ijk}\) matrix, but this procedure can be overriden
      !! replacing the `aijk` pointer procedure.
      real(pr), allocatable :: k(:, :, :) !! Attractive Binary Interatction parameter matrix
      real(pr), allocatable :: l(:, :, :) !! Repulsive Binary Interatction parameter matrix
   contains
      procedure :: aijk => kijk_constant
      !! Default attractive parameter combining rule
      procedure :: Dmix
      !! Attractive parameter mixing rule
      procedure :: Bmix
      !! Repulsive parameter mixing rule
      procedure :: D1mix => RKPR_D1mix
   end type CMR

   ! type, extends(CMR) :: CMRTD
   !    real(pr), allocatable :: k0(:, :)
   !    real(pr), allocatable :: Tref(:, :)
   ! contains
   !    procedure :: aijk => kij_exp_tdep
   ! end type CMRTD

   abstract interface
      subroutine get_aijk(&
         self, T, &
         ai, daidt, daidt2, &
         a, dadt, dadt2 &
         )
         !! Combining rule for the attractive parameter.
         !!
         !! From previously calculated attractive parameters calculate the
         !! \(a_{ij}\) matrix and it's corresponding derivatives.
         import pr, CMR
         class(CMR), intent(in) :: self
         real(pr), intent(in) :: T
         real(pr), intent(in) :: ai(:), daidt(:), daidt2(:)
         real(pr), intent(out):: a(:, :), dadt(:, :), dadt2(:, :)
      end subroutine get_aijk
   end interface

contains

   subroutine Dmix(self, n, V, T, &
      ai, daidt, daidt2, &
      D, &
      dDdV, dDdT, dDdV2, dDdT2, dDi, dDdTV, dDidV, dDidT, dDij &
      )
      !! Attractive parameter mixing rule with quadratic mix.
      !!
      !! Takes the all the pure components attractive parameters and their
      !! derivatives with respect to temperature and mix them with the
      !! Van der Waals quadratic mixing rule:
      !!
      !! \[
      !!   D = \sum_i \sum_j \sum_k n_i n_j n_k a_{ijk} = n^3 a_{mix} / n
      !! \]
      !!
      !! Inside the routine the \(a_{ijk}\) matrix is calculated using the
      !! procedure contained in the `CMR` object, this procedures defaults
      !! to the common combining rule:
      !! \(a_{ijk} = \sqrt[3]{a_i a_j a_k} (1 - k_{ijk}) \)
      !!
      !! The procedure can be overloaded by a common one that respects the
      !! interface [[get_aijk(interface)]]
      !!
      !! ```fortran
      !! type(CMR) :: my_mixing_rule
      !! my_mixing_rule%aij => new_aij_procedure
      !! ```
      class(CMR), intent(in) :: self !! Mixing rule object.
      real(pr), intent(in) :: V !! Volume [L] (unused)
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: n(:) !! Moles vector [mol]
      real(pr), intent(in) :: ai(:) !! Pure components attractive parameters \(a_i\)
      real(pr), intent(in) :: daidt(:) !! \(\frac{da_i}{dT}\)
      real(pr), intent(in) :: daidt2(:) !! \(\frac{d^2a_i}{dT^2}\)

      real(pr), intent(out) :: D !! Mixture attractive parameter \(n^2a_{mix}\)
      real(pr), intent(out) :: dDdV !! \(\frac{dD}{dT}\)
      real(pr), intent(out) :: dDdT !! \(\frac{dD}{dV}\)
      real(pr), intent(out) :: dDdT2 !! \(\frac{d^2D}{dT^2}\)
      real(pr), intent(out) :: dDdV2 !! \(\frac{d^2D}{dV^2}\)
      real(pr), intent(out) :: dDdTV !! \(\frac{d^2D}{dTV\)
      real(pr), intent(out) :: dDi(:) !! \(\frac{dD}{dn_i}\)
      real(pr), intent(out) :: dDidV(:) !! \(\frac{d^2D}{dVn_i}\)
      real(pr), intent(out) :: dDidT(:) !! \(\frac{d^2D}{dTn_i}\)
      real(pr), intent(out) :: dDij(:, :)!! \(\frac{d^2D}{dn_{ij}}\)

      integer :: i, j, nc
      real(pr) :: a(size(ai), size(ai))
      real(pr) :: dadt(size(ai), size(ai))
      real(pr) :: dadt2(size(ai), size(ai))

      nc = size(ai)

      call self%aijk(T, ai, daidt, daidt2, a, dadt, dadt2)

   end subroutine Dmix

   subroutine Bmix(self, n, bi, B, dBi, dBij)
      !! Mixture repulsive parameter.
      !!
      !! Calculate the mixture's repulsive parameter and it's derivatives
      !! with respect to composition:
      !!
      !! \[
      !!    n^2B = \sum_i \sum_j \sum_k n_i n_j n_k
      !!           \frac{b_i + b_j + b_k}{3} (1 - l_{ijk})
      !! \]
      !!
      use yaeos__models_ar_cubic_mixing_base, only: CMR_Bmix
      class(CMR), intent(in) :: self !! Mixing rule object.
      real(pr), intent(in) :: n(:) !! Moles vector.
      real(pr), intent(in) :: bi(:) !! Pure components repulsive parameters.
      real(pr), intent(out) :: B !! Mixture repulsive parameter.
      real(pr), intent(out) :: dBi(:) !! \(\frac{dB}{dn_i}\)
      real(pr), intent(out) :: dBij(:, :) !!\(\frac{d^2B}{dn_{ij}}\)

      real(pr) :: bijk(size(n), size(n), size(n))

      integer :: i, j, l, nc

      nc = size(n)

      do i=1,nc
         do j=i,nc
            do l=i,nc
               bijk(i, j, l) = (bi(i) + bi(j) + bi(l))/3._pr * (1 - self%l(i, j, l))
            end do
         end do
      end do

      call CMR_Bmix(n, bijk, b, dbi, dbij)
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
      class(CMR), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: d1i(:)
      real(pr), intent(out) :: D1
      real(pr), intent(out) :: dD1i(:)
      real(pr), intent(out) :: dD1ij(:, :)
      call d1mix_rkpr(n, d1i, d1, dd1i, dd1ij)
   end subroutine RKPR_D1mix

   subroutine kijk_constant(&
      self, T, ai, daidt, daidt2, &
      a, dadt, dadt2 &
      )
      !! Combining rule that uses constant \(k_{ijk}\) values.
      !!
      !! \[
      !!  a_{ijk} = \sqrt[3]{a_i a_j a_k} (1 - k_{ijk})
      !! ]
      use yaeos__models_ar_cubic_mixing_base, only: CMR_aijk
      class(CMR), intent(in) :: self
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: ai(:) !! Pure components attractive parameters (\a_i\)
      real(pr), intent(in) :: daidt(:) !! \(\frac{da_i}{dT}\)
      real(pr), intent(in) :: daidt2(:) !! \(\frac{d^2a_i}{dT^2}\)
      real(pr), intent(out) :: a(:, :) !! \(a_{ij}\) Matrix
      real(pr), intent(out) :: dadt(:, :) !! \(\frac{da_{ij}{dT}\)
      real(pr), intent(out) :: dadt2(:, :)!! \(\frac{d^2a_{ij}{dT^2}\)

      integer :: i, j, l

      real(pr) :: k(size(ai), size(ai), size(ai))
      real(pr) :: zeros(size(ai), size(ai), size(ai))

      k = self%k
      zeros = 0
      call CMR_aijk(ai, daidt, daidt2, k, zeros, zeros, a, dadt, dadt2)
   end subroutine kijk_constant

   ! subroutine kijk_exp_tdep(&
   !    self, T, a, dadt, dadt2, &
   !    aij, daijdt, daijdt2 &
   !    )
   !    !! # kij_exp_tdep
   !    !!
   !    !! Combining rule that uses temperature dependant \(k_{ij}\) values.
   !    !! With the following expression:
   !    !! \[
   !    !! k_{ij}(T) = k_{ij}^0 + k_{ij}^\infty \exp\left(\frac{-T}{T^*}\right)
   !    !!  \]
   !    !!
   !    !! \[
   !    !!  a_{ij} = \sqrt{a_i a_j} (1 - k_{ij})
   !    !! ]
   !    use hyperdual_mod
   !    class(QMRTD), intent(in) :: self
   !    real(pr), intent(in) :: T !! Temperature [K]
   !    real(pr), intent(in) :: a(:) !! Pure components attractive parameters (\a_i\)
   !    real(pr), intent(in) :: dadt(:) !! \(\frac{da_i}{dT}\)
   !    real(pr), intent(in) :: dadt2(:) !! \(\frac{d^2a_i}{dT^2}\)
   !    real(pr), intent(out) :: aij(:, :) !! \(a_{ij}\) Matrix
   !    real(pr), intent(out) :: daijdt(:, :) !! \(\frac{da_{ij}{dT}\)
   !    real(pr), intent(out) :: daijdt2(:, :)!! \(\frac{d^2a_{ij}{dT^2}\)

   !    real(pr) :: k0(size(a), size(a))
   !    real(pr) :: kinf(size(a), size(a))
   !    real(pr) :: Tstar(size(a), size(a))

   !    type(hyperdual) :: aij_hd(size(a), size(a)), kij_hd(size(a), size(a)), T_hd
   !    type(hyperdual) :: a_hd(size(a))

   !    integer :: i, j, nc

   !    T_hd = T
   !    T_hd%f1 = 1
   !    T_hd%f2 = 1

   !    k0 = self%k0
   !    kinf = self%k

   !    Tstar = self%Tref

   !    kij_hd = kinf + k0 * exp(-T_hd / Tstar)

   !    a_hd = a

   !    ! Inject the already calculated derivatives
   !    a_hd%f1 = dadt
   !    a_hd%f2 = dadt
   !    a_hd%f12 = dadt2

   !    nc = size(a)

   !    do i=1,size(a)
   !       aij_hd(i, i) = sqrt(a_hd(i) * a_hd(i))
   !       do j=i+1,size(a)
   !          aij_hd(i, j) = sqrt(a_hd(i) * a_hd(j)) * (1._pr - kij_hd(i, j))
   !          aij_hd(j, i) = aij_hd(i, j)
   !       end do
   !    end do

   !    aij = aij_hd%f0
   !    daijdt = aij_hd%f1
   !    daijdt2 = aij_hd%f12
   ! end subroutine kij_exp_tdep
end module yaeos__models_ar_cubic_quadratic_mixing
