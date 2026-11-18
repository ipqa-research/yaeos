module yaeos__models_ar_cubic_cubic_mixing
   !! # Cubic mixing rules
   !! Three-component combining rules for cubic equations of state.
   !!
   !! # Description
   !! `CMR` combines pure-component parameters and ternary interaction
   !! coefficients. `CMRTD` adds temperature dependence to its attractive
   !! interaction parameters. These types implement the `CubicMixRule` API.
   !!
   !! # Examples
   !! ```fortran
   !! rule = CMR(k=kijk, l=l3)
   !! call model%set_mixrule(rule)
   !! ```
   use yaeos__constants, only: pr, solving_volume
   use yaeos__substance, only: substances
   use yaeos__models_ar_genericcubic, only: CubicMixRule
   use yaeos__models_ar_cubic_mixing_base, only: bmix_qmr
   implicit none

   private

   public :: CMR
   public :: CMRTD
   public :: kijk_exp_tdep

   type, extends(CubicMixRule) :: CMR
      !! # CMR
      !! Cubic mixing rule with ternary component interactions.
      !!
      !! # Description
      !! CMR depends on ternary interaction parameters; on a Cubic EoS
      !! the mixture is obtained by the combination of attractive and
      !! repulsive parameter tensors.
      !!
      !! By default the attractive parameter matrix is calculated with:
      !! \[a_{ijk} = \sqrt[3]{a_i a_j a_k}(1 - k_{ijk})\]
      !! generating the \(a_{ijk}\) tensor, but this procedure can be overridden
      !! in an extended type by rebinding the `aijk` type-bound procedure.
      !!
      !! # Examples
      !! ```fortran
      !! rule = CMR(k=kijk, l=l3)
      !! call model%set_mixrule(rule)
      !! ```
      real(pr), allocatable :: k(:, :, :) !! Attractive ternary interaction parameters \(k_{ijk}\)
      real(pr), allocatable :: l(:, :, :) !! Repulsive ternary interaction parameters \(l_{ijk}\)
   contains
      procedure :: aijk => kijk_constant
      !! Default attractive parameter combining rule
      procedure :: Dmix
      !! Attractive parameter mixing rule
      procedure :: Bmix
      !! Repulsive parameter mixing rule
      procedure :: D1mix => RKPR_D1mix
   end type CMR

   type, extends(CMR) :: CMRTD
      !! # CMRTD
      !! Cubic mixing rule with temperature-dependent ternary interactions.
      !!
      !! # Description
      !! Extends `CMR` with temperature-dependent three-component attractive
      !! interaction parameters,
      !! \(k_{ijk}(T) = k_{ijk} + k^0_{ijk}\exp(-T/T^{ref}_{ijk})\).
      !!
      !! # Examples
      !! ```fortran
      !! rule = CMRTD(k=k, k0=k0, Tref=Tref, l=l)
      !! ```
      real(pr), allocatable :: k0(:, :, :) !! Temperature-dependent ternary interaction amplitudes.
      real(pr), allocatable :: Tref(:, :, :) !! Reference temperatures for component triplets [K].
   contains
      procedure :: aijk => kijk_exp_tdep
   end type CMRTD

   abstract interface
      subroutine get_aijk(&
         self, T, &
         ai, daidt, daidt2, &
         a, dadt, dadt2 &
         )
         !! Combining rule for the attractive parameter.
         !!
         !! From previously calculated attractive parameters calculate the
         !! \(a_{ijk}\) tensor and its corresponding derivatives.
         import pr, CMR
         class(CMR), intent(in) :: self !! Mixing rule object.
         real(pr), intent(in) :: T !! Temperature [K]
         real(pr), intent(in) :: ai(:) !! Pure-component attractive parameters \(a_i\)
         real(pr), intent(in) :: daidt(:) !! \(da_i/dT\)
         real(pr), intent(in) :: daidt2(:) !! \(d^2a_i/dT^2\)
         real(pr), intent(out) :: a(:, :) !! \(a_{ijk}\) values
         real(pr), intent(out) :: dadt(:, :) !! \(da_{ijk}/dT\)
         real(pr), intent(out) :: dadt2(:, :) !! \(d^2a_{ijk}/dT^2\)
      end subroutine get_aijk
   end interface

contains

   subroutine Dmix(self, n, V, T, &
      ai, daidt, daidt2, &
      D, &
      dDdV, dDdT, dDdV2, dDdT2, dDi, dDdTV, dDidV, dDidT, dDij &
      )
      !! # Dmix
      !!
      !! # Description
      !! Attractive parameter mixing rule with cubic mix.
      !!
      !! Takes the all the pure components attractive parameters and their
      !! derivatives with respect to temperature and mix them with the
      !! cubic (three-body) mixing rule:
      !!
      !! \[
      !!   D = \frac{1}{n_T}\sum_i \sum_j \sum_k n_i n_j n_k a_{ijk} = n_T^2 a_{mix}
      !! \]
      !!
      !! Inside the routine the \(a_{ijk}\) matrix is calculated using the
      !! procedure contained in the `CMR` object, this procedures defaults
      !! to the common combining rule:
      !! \(a_{ijk} = \sqrt[3]{a_i a_j a_k} (1 - k_{ijk}) \)
      !!
      !! The combining rule can be replaced by extending `CMR` and rebinding
      !! `aijk` to a procedure with the same arguments as `kijk_constant`:
      !!
      !! ```fortran
      !! type, extends(CMR) :: MyCMR
      !! contains
      !!    procedure :: aijk => my_aijk
      !! end type MyCMR
      !! ```
      use yaeos__models_ar_cubic_mixing_base, only: CMR_Dmix
      class(CMR), intent(in) :: self !! Mixing rule object.
      real(pr), intent(in) :: V !! Volume [L] (unused)
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: n(:) !! Moles vector [mol]
      real(pr), intent(in) :: ai(:) !! Pure components attractive parameters \(a_i\)
      real(pr), intent(in) :: daidt(:) !! \(\frac{da_i}{dT}\)
      real(pr), intent(in) :: daidt2(:) !! \(\frac{d^2a_i}{dT^2}\)

      real(pr), intent(out) :: D !! Mixture attractive parameter \(n^2a_{mix}\)
      real(pr), intent(out) :: dDdV !! \(\frac{dD}{dV}\) (zero)
      real(pr), intent(out) :: dDdT !! \(\frac{dD}{dT}\)
      real(pr), intent(out) :: dDdT2 !! \(\frac{d^2D}{dT^2}\)
      real(pr), intent(out) :: dDdV2 !! \(\frac{d^2D}{dV^2}\)
      real(pr), intent(out) :: dDdTV !! \(\frac{d^2D}{dTV}\)
      real(pr), intent(out) :: dDi(:) !! \(\frac{dD}{dn_i}\)
      real(pr), intent(out) :: dDidV(:) !! \(\frac{d^2D}{dVn_i}\)
      real(pr), intent(out) :: dDidT(:) !! \(\frac{d^2D}{dTn_i}\)
      real(pr), intent(out) :: dDij(:, :) !! \(\frac{d^2D}{dn_{ij}}\)

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.
      integer :: nc !! Number of components.
      real(pr) :: a(size(ai), size(ai), size(ai))
      real(pr) :: dadt(size(ai), size(ai), size(ai))
      real(pr) :: dadt2(size(ai), size(ai), size(ai))

      nc = size(ai)

      call self%aijk(T, ai, daidt, daidt2, a, dadt, dadt2)
      call CMR_Dmix(&
         n, V, T, &
         a, dadt, dadt2, &
         D=D, &
         dDdV=dDdV, &
         dDdT=dDdT, &
         dDdT2=dDdT2, &
         dDdV2=dDdV2, &
         dDdTV=dDdTV, &
         dDi=dDi, &
         dDidV=dDidV, &
         dDidT=dDidT, &
         dDij=dDij &
         )
   end subroutine Dmix

   subroutine Bmix(self, n, bi, B, dBi, dBij)
      !! # Bmix
      !!
      !! # Description
      !! Mixture repulsive parameter.
      !!
      !! Calculate the mixture's repulsive parameter and it's derivatives
      !! with respect to composition:
      !!
      !! \[
      !!    n_T^2 B = \sum_i \sum_j \sum_k n_i n_j n_k
      !!           \frac{b_i + b_j + b_k}{3} (1 - l_{ijk})
      !! \]
      !!
      use yaeos__models_ar_cubic_mixing_base, only: CMR_Bmix
      class(CMR), intent(in) :: self !! Mixing rule object.
      real(pr), intent(in) :: n(:) !! Moles vector.
      real(pr), intent(in) :: bi(:) !! Pure-component covolumes [L/mol]
      real(pr), intent(out) :: B !! Mixture covolume \(B\) [L]
      real(pr), intent(out) :: dBi(:) !! \(\frac{dB}{dn_i}\)
      real(pr), intent(out) :: dBij(:, :) !! \(\frac{d^2B}{dn_{ij}}\)

      real(pr) :: bijk(size(n), size(n), size(n))

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.
      integer :: l !! Loop index `l`.
      integer :: nc !! Number of components.

      nc = size(n)

      do i=1,nc
         do j=i,nc
            do l=i,nc
               bijk(i, j, l) = (bi(i) + bi(j) + bi(l))/3._pr * (1 - self%l(i, j, l))
               bijk(i, l, j) = bijk(i, j, l)
               
               bijk(j, i, l) = bijk(i, j, l)
               bijk(j, l, i) = bijk(i, j, l)
               
               bijk(l, i, j) = bijk(i, j, l)
               bijk(l, j, i) = bijk(i, j, l)
            end do
         end do
      end do

      call CMR_Bmix(n, bijk, b, dbi, dbij)
   end subroutine Bmix

   subroutine RKPR_D1mix(self, n, d1i, D1, dD1i, dD1ij)
      !! # RKPR_D1mix
      !!
      !! # Description
      !! RKPR \(\delta_1\) parameter mixing rule.
      !!
      !! The RKPR EoS doesn't have a constant \(\delta_1\) value for each
      !! component, so a proper mixing rule should be provided. A linear
      !! combination is used.
      !!
      !! \[
      !!     \Delta_1 = \frac{\sum_i^N n_i \delta_{1i}}{\sum_i^N n_i}
      !! \]
      !!
      use yaeos__models_ar_cubic_mixing_base, only: d1mix_rkpr
      class(CMR), intent(in) :: self !! Mixing rule object.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol]
      real(pr), intent(in) :: d1i(:) !! Pure-component \(\delta_{1i}\) values
      real(pr), intent(out) :: D1 !! Mixture \(\Delta_1\)
      real(pr), intent(out) :: dD1i(:) !! \(\partial \Delta_1/\partial n_i\)
      real(pr), intent(out) :: dD1ij(:, :) !! \(\partial^2 \Delta_1/\partial n_i \partial n_j\)
      call d1mix_rkpr(n, d1i, d1, dd1i, dd1ij)
   end subroutine RKPR_D1mix

   subroutine kijk_constant(&
      self, T, ai, daidt, daidt2, &
      a, dadt, dadt2 &
      )
      !! # kijk_constant
      !!
      !! # Description
      !! Combining rule that uses constant \(k_{ijk}\) values.
      !!
      !! \[
      !!  a_{ijk} = \sqrt[3]{a_i a_j a_k} (1 - k_{ijk})
      !! \]
      use yaeos__models_ar_cubic_mixing_base, only: CMR_aijk
      class(CMR), intent(in) :: self !! Mixing rule object.
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: ai(:) !! Pure components attractive parameters \(a_i\)
      real(pr), intent(in) :: daidt(:) !! \(\frac{da_i}{dT}\)
      real(pr), intent(in) :: daidt2(:) !! \(\frac{d^2a_i}{dT^2}\)
      real(pr), intent(out) :: a(:, :, :) !! \(a_{ijk}\) Matrix
      real(pr), intent(out) :: dadt(:, :, :) !! \(\frac{da_{ijk}}{dT}\)
      real(pr), intent(out) :: dadt2(:, :, :) !! \(\frac{d^2a_{ijk}}{dT^2}\)

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.
      integer :: l !! Loop index `l`.

      real(pr) :: k(size(ai), size(ai), size(ai)) !! Ternary interaction parameters \(k_{ijk}\)
      real(pr) :: zeros(size(ai), size(ai), size(ai))

      k = self%k
      zeros = 0
      call CMR_aijk(ai, daidt, daidt2, k, zeros, zeros, a, dadt, dadt2)
   end subroutine kijk_constant

   subroutine kijk_exp_tdep(&
      self, T, ai, daidt, daidt2, &
      a, dadt, dadt2 &
      )
      !! # kijk_exp_tdep
      !!
      !! # Description
      !! Combining rule that uses \(k_{ijk}\) as a function of temperature.
      !!
      !! \[
      !!  k_{ijk}(T) = k_{ijk} + k_{ijk}^0 \exp \left(-T/T^{*}_{ijk}\right)
      !! \]
      !!
      !! \[
      !!  a_{ijk} = \sqrt[3]{a_i a_j a_k} (1 - k_{ijk})
      !! \]
      use yaeos__models_ar_cubic_mixing_base, only: CMR_aijk
      class(CMRTD), intent(in) :: self !! Mixing rule object.
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: ai(:) !! Pure components attractive parameters \(a_i\)
      real(pr), intent(in) :: daidt(:) !! \(\frac{da_i}{dT}\)
      real(pr), intent(in) :: daidt2(:) !! \(\frac{d^2a_i}{dT^2}\)
      real(pr), intent(out) :: a(:, :, :) !! \(a_{ijk}\) Matrix
      real(pr), intent(out) :: dadt(:, :, :) !! \(\frac{da_{ijk}}{dT}\)
      real(pr), intent(out) :: dadt2(:, :, :) !! \(\frac{d^2a_{ijk}}{dT^2}\)

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.
      integer :: l !! Loop index `l`.
      integer :: nc !! Number of components.

      real(pr) :: c(size(ai), size(ai), size(ai))
      real(pr) :: k(size(ai), size(ai), size(ai)) !! Ternary interaction parameters \(k_{ijk}\)
      real(pr) :: dkdt(size(ai), size(ai), size(ai))
      real(pr) :: dkdt2(size(ai), size(ai), size(ai))

      nc = size(ai)
      k = 0
      dkdt = 0
      dkdt2 = 0

      where(self%k0 /= 0 .or. self%k /= 0)
         c = self%k0 * exp(-T/self%Tref)
         k = self%k + c
         dkdt = -c / self%Tref
         dkdt2 = c / self%Tref**2
      end where

      call CMR_aijk(ai, daidt, daidt2, k, dkdt, dkdt2, a, dadt, dadt2)
   end subroutine kijk_exp_tdep
end module yaeos__models_ar_cubic_cubic_mixing
