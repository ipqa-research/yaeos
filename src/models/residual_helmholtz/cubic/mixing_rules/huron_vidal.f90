module yaeos__models_cubic_mixing_rules_huron_vidal
   !! # Huron-Vidal cubic mixing rules
   !! Cubic-EOS mixing rules that couple an excess-Gibbs model to the attractive
   !! mixture parameter.
   !!
   !! # Description
   !! This module implements Huron-Vidal, Michelsen-Huron-Vidal (MHV), and an
   !! NRTL-HV variant. These rules use an excess-Gibbs model to represent
   !! non-ideal liquid interactions while retaining the cubic equation of state
   !! for volumetric behavior.
   !!
   !! # Examples
   !! ```fortran
   !! ge = NRTL(A, B, alpha)
   !! mixrule = MHV(ge=ge, b=b, q=-0.593_pr)
   !! call model%set_mixrule(mixrule)
   !! ```
   !!
   !! # References
   !! - Huron, M.-J., and J. Vidal, *Fluid Phase Equilibria* 3 (1979) 255-271.
   !! - Michelsen, M. L., *Fluid Phase Equilibria* 60 (1990) 213-219.
   use yaeos__constants, only: pr, R, solving_volume
   use yaeos__models_ar_genericcubic, only: CubicMixRule
   use yaeos__models_ar_cubic_quadratic_mixing, only: QMR
   use yaeos__models_ar_cubic_mixing_base, only: bmix_qmr
   use yaeos__models_ge, only: GeModel
   use yaeos__models_ge_nrtlhv, only: NRTLHV
   implicit none

   private

   public :: HV
   public :: MHV
   public :: DmixMHV
   public :: HV_NRTL, init_hvnrtl

   type, extends(CubicMixRule) :: HV
      !! # HV
      !! Huron-Vidal cubic mixing rule.
      !!
      !! # Description
      !! Couples an excess-Gibbs model with a cubic EOS through the
      !! infinite-pressure-limit relation. Component covolumes (`bi`) and
      !! \(\delta_1\) values (`del1`) are stored in the type; there is no
      !! custom constructor, so use the structure constructor.
      !!
      !! # Examples
      !! The initialized variants exported here are `MHV` and `HV_NRTL`.
      class(GeModel), allocatable :: ge !! Excess-Gibbs model used for mixture non-ideality.
      real(pr), allocatable :: del1(:) !! Pure-component cubic \(\delta_1\) parameters.
      real(pr), allocatable :: bi(:) !! Pure-component covolumes [L/mol].
   contains
      procedure :: Bmix => BmixHV
      procedure :: D1Mix => D1MixHV
      procedure :: Dmix => DmixHV
   end type HV

   type, extends(CubicMixRule) :: MHV
      !! # MHV
      !! Michelsen's modified Huron-Vidal cubic mixing rule.
      !!
      !! # Description
      !! A zero-pressure-limit mixing rule coupling a cubic EOS with an
      !! excess-Gibbs model. The parameter `q` selects the MHV variant, and
      !! `l` can hold repulsive binary interactions.
      !! This mixing rule is based on the approximate zero-pressure limit
      !! of a cubic equation of state. At the approximate zero-pressure limit the
      !! attractive parameter can be expressed as:
      !!
      !! \[
      !! \frac{D}{RTB}(n, T) = \sum_i n_i \frac{a_i(T)}{RT b_i} + \frac{1}{q}
      !!  \left(\frac{G^E(n, T)}{RT} + \sum_i n_i \ln \frac{B}{n_T b_i} \right)
      !! \]
      !! Where \(q\) is a weak function of temperature. In this implementation
      !! it is a constant chosen for the cubic EOS (e.g. \(-0.593\) for SRK, as
      !! in MHV1, and \(-0.64663\) in PSRK).
      !!
      !! # Examples
      !! To use the modified Huron-Vidal mixing rule it is necessary to define
      !! a `CubicEoS` and replace its original mixing rule with the one generated
      !! by the user.
      !! ```fortran
      !! type(MHV) :: mixrule
      !! type(NRTL) :: ge_model
      !! type(CubicEoS) :: model
      !!
      !! ! Define the Ge model to be used and the CubicEoS
      !! ge_model = NRTL(a, b, c)
      !! model = SoaveRedlichKwong(tc, pc, w)
      !!
      !! ! Use the initialization function to setup
      !! mixrule = MHV(ge=ge_model, q=-0.593_pr, b=model%b)
      !!
      !! ! Replace the original mixrule on the previously defined model
      !! model%mixrule = mixrule
      !!
      !! ! Ready to do calculations
      !! call model%pressure(n, V, T, P)
      !! ```
      !!
      !! # References
      !! Michelsen, M. L., *Fluid Phase Equilibria* 60 (1990) 213-219.
      real(pr), allocatable :: l(:, :) !! Repulsive binary interaction matrix.
      real(pr), private, allocatable :: bi(:) !! Cached pure-component covolumes.
      real(pr), private, allocatable :: B !! Mixture covolume.
      real(pr), private, allocatable :: dBi(:) !! Mole-number gradient of the mixture covolume.
      real(pr), private, allocatable :: dBij(:, :) !! Mole-number Hessian of the mixture covolume.
      class(GeModel), allocatable :: ge !! Excess-Gibbs model in the MHV relation.
      real(pr) :: q !! MHV \(q\) parameter defining the zero-pressure relation.
   contains
      procedure :: Bmix => BmixMHV
      procedure :: D1Mix => D1MixMHV
      procedure :: Dmix => DmixMHV
   end type MHV

   type, extends(CubicMixRule) :: HV_NRTL
      !! # HV_NRTL
      !! Huron-Vidal cubic mixing rule with an NRTL-HV excess-Gibbs model.
      !!
      !! # Description
      !! This variant evaluates the Huron-Vidal relation using the modified
      !! NRTL expression. Under the limiting parameter choices below, its
      !! attractive mixing expression reduces to the quadratic form.
      !!
      !! \[
      !!  \alpha_{ji} = 0
      !! \]
      !!
      !! \[
      !!  g_{ii} = -\frac{a_i}{b_i} \lambda
      !! \]
      !!
      !! \[
      !!  g_{ji} = -2\frac{\sqrt{b_i b_j}}{b_i + b_j} \sqrt{g_{ii}g_{jj}}
      !!           \left(1 - k_{ij}\right)
      !! \]
      !!
      !! # Examples
      !! ```fortran
      !! mixrule = HV_NRTL(b, del1, alpha, gji0, gjiT, use_kij, kij)
      !! call model%set_mixrule(mixrule)
      !! ```
      type(NRTLHV) :: ge !! NRTL-HV excess-Gibbs model.
      real(pr), allocatable :: del1(:) !! Pure-component cubic \(\delta_1\) parameters.
      real(pr), allocatable :: bi(:) !! Pure-component covolumes [L/mol].
      logical, allocatable :: use_kij(:, :) !! Pair mask: where true, the NRTL parameters are replaced so that the rule reduces to the quadratic one with `kij`.
      real(pr), allocatable :: kij(:, :) !! Binary attractive interaction parameters.
   contains
      procedure :: Bmix => BmixHVNRTL
      procedure :: D1Mix => D1MixHVNRTL
      procedure :: Dmix => DmixHVNRTL
   end type HV_NRTL

   interface MHV
      module procedure :: init_mhv
   end interface MHV

   interface HV_NRTL
      module procedure :: init_hvnrtl
   end interface HV_NRTL

contains

   ! ===========================================================================
   ! Huron-Vidal Mixing rule
   ! ---------------------------------------------------------------------------
   subroutine BmixHV(self, n, bi, B, dBi, dBij)
      !! # Repulsive parameter \(B\) mixing rule
      !! Linear mixing rule for the repulsive parameter.
      !!
      !! # Description
      !! \[B = \sum_i n_i b_i\]
      use yaeos__models_ar_cubic_mixing_base, only: bmix_linear
      class(HV), intent(in) :: self !! Mixing rule object.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol]
      real(pr), intent(in) :: bi(:) !! Pure-component covolumes [L/mol]
      real(pr), intent(out) :: B !! Mixture covolume \(B\) [L]
      real(pr), intent(out) :: dBi(:) !! \(\partial B/\partial n_i\)
      real(pr), intent(out) :: dBij(:, :) !! \(\partial^2 B/\partial n_i \partial n_j\)
      call bmix_linear(n, bi, b, dbi, dbij)
   end subroutine BmixHV

   subroutine D1MixHV(self, n, d1i, D1, dD1i, dD1ij)
      !! # D1MixHV
      !! Mole-fraction-weighted mixture \(\delta_1\) (see `d1mix_rkpr`).
      use yaeos__models_ar_cubic_mixing_base, only: d1mix_rkpr
      class(HV), intent(in) :: self !! Mixing rule object.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol]
      real(pr), intent(in) :: d1i(:) !! Pure-component \(\delta_{1i}\) values
      real(pr), intent(out) :: D1 !! Mixture \(\delta_1\)
      real(pr), intent(out) :: dD1i(:) !! \(\partial \delta_1/\partial n_i\)
      real(pr), intent(out) :: dD1ij(:, :) !! \(\partial^2 \delta_1/\partial n_i \partial n_j\)
      call d1mix_rkpr(n, d1i, D1, dD1i, dD1ij)
   end subroutine D1MixHV

   subroutine DmixHV(self, n, V, T, &
      ai, daidt, daidt2, &
      D, &
      dDdV, dDdT, dDdV2, dDdT2, dDi, dDdTV, dDidV, dDidT, dDij &
      )
      !! # DmixHV
      !! Huron-Vidal attractive parameter using the stored `GeModel`.
      !!
      !! # Description
      !! Evaluates \(B\), \(\delta_1\), and \(G^E\) and its derivatives, then
      !! applies \(D = B(\sum_i n_i a_i/b_i - G^E/\Lambda)\) through the
      !! module-level `DmixHV` kernel. Volume derivatives are not assigned.
      use yaeos__models_ar_cubic_mixing_base, only: lamdba_hv, mix => DMixHV
      class(HV), intent(in) :: self !! Mixing rule object.
      real(pr), intent(in) :: V !! Volume [L] (unused)
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: n(:) !! Moles vector [mol]
      real(pr), intent(in) :: ai(:) !! Pure components attractive parameters \(a_i\)
      real(pr), intent(in) :: daidt(:) !! \(\frac{da_i}{dT}\)
      real(pr), intent(in) :: daidt2(:) !! \(\frac{d^2a_i}{dT^2}\)

      real(pr), intent(out) :: D !! Mixture attractive parameter \(n^2a_{mix}\)
      real(pr), intent(out) :: dDdV !! \(\frac{dD}{dV}\) (not assigned)
      real(pr), intent(out) :: dDdT !! \(\frac{dD}{dT}\)
      real(pr), intent(out) :: dDdT2 !! \(\frac{d^2D}{dT^2}\)
      real(pr), intent(out) :: dDdV2 !! \(\frac{d^2D}{dV^2}\)
      real(pr), intent(out) :: dDdTV !! \(\frac{d^2D}{dTV}\)
      real(pr), intent(out) :: dDi(:) !! \(\frac{dD}{dn_i}\)
      real(pr), intent(out) :: dDidV(:) !! \(\frac{d^2D}{dVn_i}\)
      real(pr), intent(out) :: dDidT(:) !! \(\frac{d^2D}{dTn_i}\)
      real(pr), intent(out) :: dDij(:, :) !! \(\frac{d^2D}{dn_{ij}}\)

      real(pr) :: b
      real(pr) :: bi(size(n))
      real(pr) :: dbi(size(n))
      real(pr) :: dbij(size(n), size(n))
      real(pr) :: del1(size(n))
      real(pr) :: del2(size(n))
      real(pr) :: d1
      real(pr) :: d1i(size(n))
      real(pr) :: dd1i(size(n))
      real(pr) :: dd1ij(size(n), size(n))
      real(pr) :: Ge !! Excess Gibbs energy.
      real(pr) :: GeT !! Temperature derivative of excess Gibbs energy.
      real(pr) :: GeT2 !! Second temperature derivative of excess Gibbs energy.
      real(pr) :: Gen(size(n)) !! Mole-number gradient of excess Gibbs energy.
      real(pr) :: GeTn(size(n)) !! Mixed temperature/mole-number derivative of excess Gibbs energy.
      real(pr) :: Gen2(size(n), size(n)) !! Mole-number Hessian of excess Gibbs energy.

      real(pr) :: totn !! Total number of moles

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.
      integer :: nc !! Number of components.
      real(pr) :: L !! Infinite-pressure parameter \(\Lambda\)
      real(pr) :: dL(size(n))
      real(pr) :: dL2(size(n), size(n))

      nc = size(n)
      totn = sum(n)

      del1 = self%del1
      bi = self%bi

      call self%Bmix(n, bi, B, dBi, dBij)
      call self%D1Mix(n, del1, D1, dD1i, dD1ij)
      call lamdba_hv(nc, D1, dD1i, dD1ij, L, dL, dL2)

      call self%ge%excess_gibbs( &
         n, T, Ge=Ge, GeT=GeT, GeT2=GeT2, Gen=Gen, GeTn=GeTn, Gen2=Gen2 &
         )
      call mix(n, T,&
         bi, B, dBi, dBij, &
         D1, dD1i, dD1ij, &
         ai, daidt, daidt2, &
         Ge, GeT, GeT2, Gen, GeTn, Gen2,&
         D, dDdT, dDdT2, dDi, dDidT, dDij)
   end subroutine DmixHV

   ! ===========================================================================
   ! Huron-Vidal Mixing rule with Huron-Vidal NRTL
   ! ---------------------------------------------------------------------------
   type(HV_NRTL) function init_hvnrtl(b, del1, alpha, gji0, gjiT, use_kij, kij) result(mixrule)
      !! # Huron-Vidal NRTL mixing rule
      !! This is the Huron-Vidal mixing rule that includes the NRTL model
      !! modified by Huron and Vidal.
      !!
      !! # Description
      !! This is the Huron-Vidal mixing rule that includes the NRTL model
      !! modified by Huron and Vidal. It is a mixing rule that allows to
      !! use the NRTL model as an excess Gibbs energy model and can be
      !! simplified to the classic quadratic mixing rules when the parameters
      !! are set to:
      !! \[
      !!  \alpha_{ji} = 0
      !! \]
      !!
      !! \[
      !!  g_{ii} = -\frac{a_i}{b_i} \lambda
      !! \]
      !!
      !! \[
      !!  g_{ji} = -2\frac{\sqrt{b_i b_j}}{b_i + b_j} \sqrt{g_{ii}g_{jj}}
      !!           \left(1 - k_{ij}\right)
      !! \]
      !!
      !! # Examples
      !! ```fortran
      !! mixrule = HV_NRTL(b, del1, alpha, gji0, gjiT, use_kij, kij)
      !! ```
      use yaeos__models_ge_nrtlhv, only: NRTLHV
      real(pr), intent(in) :: b(:) !! Pure-component covolumes [L/mol]
      real(pr), intent(in) :: del1(:) !! Pure-component \(\delta_1\) values
      real(pr), intent(in) :: alpha(:, :) !! NRTL non-randomness parameters \(\alpha_{ji}\)
      real(pr), intent(in) :: gji0(:, :) !! Temperature-independent NRTL-HV energy parameters
      real(pr), intent(in) :: gjiT(:, :) !! Temperature coefficients of the NRTL-HV energy parameters
      logical, intent(in) :: use_kij(:, :) !! Pairs where the quadratic-limit parameters with `kij` are used
      real(pr), intent(in) :: kij(:, :) !! Binary attractive interaction parameters

      integer :: i !! Loop/index variable.
      integer :: nc !! Number of components.

      nc = size(b)

      mixrule%ge = NRTLHV(b=b, alpha=alpha, gji0=gji0, gjiT=gjiT)
      mixrule%bi = b
      mixrule%del1 = del1
      mixrule%use_kij = use_kij
      mixrule%kij = kij
   end function init_hvnrtl

   subroutine BmixHVNRTL(self, n, bi, B, dBi, dBij)
      !! # Repulsive parameter \(B\) mixing rule
      !! Linear mixing rule for the repulsive parameter.
      !!
      !! # Description
      !! \[B = \sum_i n_i b_i\]
      use yaeos__models_ar_cubic_mixing_base, only: bmix_linear
      class(HV_NRTL), intent(in) :: self !! Mixing rule object.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol]
      real(pr), intent(in) :: bi(:) !! Pure-component covolumes [L/mol]
      real(pr), intent(out) :: B !! Mixture covolume \(B\) [L]
      real(pr), intent(out) :: dBi(:) !! \(\partial B/\partial n_i\)
      real(pr), intent(out) :: dBij(:, :) !! \(\partial^2 B/\partial n_i \partial n_j\)
      call bmix_linear(n, bi, b, dbi, dbij)
   end subroutine BmixHVNRTL

   subroutine D1MixHVNRTL(self, n, d1i, D1, dD1i, dD1ij)
      !! # D1MixHVNRTL
      !! Mole-fraction-weighted mixture \(\delta_1\) (see `d1mix_rkpr`).
      use yaeos__models_ar_cubic_mixing_base, only: d1mix_rkpr
      class(HV_NRTL), intent(in) :: self !! Mixing rule object.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol]
      real(pr), intent(in) :: d1i(:) !! Pure-component \(\delta_{1i}\) values
      real(pr), intent(out) :: D1 !! Mixture \(\delta_1\)
      real(pr), intent(out) :: dD1i(:) !! \(\partial \delta_1/\partial n_i\)
      real(pr), intent(out) :: dD1ij(:, :) !! \(\partial^2 \delta_1/\partial n_i \partial n_j\)
      call d1mix_rkpr(n, d1i, D1, dD1i, dD1ij)
   end subroutine D1MixHVNRTL

   subroutine DmixHVNRTL(self, n, V, T, &
      ai, daidt, daidt2, &
      D, &
      dDdV, dDdT, dDdV2, dDdT2, dDi, dDdTV, dDidV, dDidT, dDij &
      )
      !! # DmixHVNRTL
      !! Huron-Vidal attractive parameter with the NRTL-HV model.
      !!
      !! # Description
      !! For pairs flagged in `use_kij`, the NRTL-HV parameters are replaced
      !! (\(\alpha_{ij}=0\), `gjiT`=0, and `gji0` from \(g_{ii} = -\Lambda a_i/b_i\)
      !! and `kij`) before evaluating \(G^E\) and the Huron-Vidal relation.
      use yaeos__models_ar_cubic_mixing_base, only: lamdba_hv, DmixHV
      use yaeos__models_ge_nrtlhv, only: NRTLHV
      class(HV_NRTL), intent(in) :: self !! Mixing rule object.
      real(pr), intent(in) :: V !! Volume [L] (unused)
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: n(:) !! Moles vector [mol]
      real(pr), intent(in) :: ai(:) !! Pure components attractive parameters \(a_i\)
      real(pr), intent(in) :: daidt(:) !! \(\frac{da_i}{dT}\)
      real(pr), intent(in) :: daidt2(:) !! \(\frac{d^2a_i}{dT^2}\)

      real(pr), intent(out) :: D !! Mixture attractive parameter \(n^2a_{mix}\)
      real(pr), intent(out) :: dDdV !! \(\frac{dD}{dV}\) (not assigned)
      real(pr), intent(out) :: dDdT !! \(\frac{dD}{dT}\)
      real(pr), intent(out) :: dDdT2 !! \(\frac{d^2D}{dT^2}\)
      real(pr), intent(out) :: dDdV2 !! \(\frac{d^2D}{dV^2}\)
      real(pr), intent(out) :: dDdTV !! \(\frac{d^2D}{dTV}\)
      real(pr), intent(out) :: dDi(:) !! \(\frac{dD}{dn_i}\)
      real(pr), intent(out) :: dDidV(:) !! \(\frac{d^2D}{dVn_i}\)
      real(pr), intent(out) :: dDidT(:) !! \(\frac{d^2D}{dTn_i}\)
      real(pr), intent(out) :: dDij(:, :) !! \(\frac{d^2D}{dn_{ij}}\)

      real(pr) :: Ge
      real(pr) :: GeT
      real(pr) :: GeT2 !! Second-order intermediate quantity used by this procedure.
      real(pr) :: Gen(size(n))
      real(pr) :: GeTn(size(n))
      real(pr) :: Gen2(size(n), size(n)) !! Second-order intermediate quantity used by this procedure.

      real(pr) :: B
      real(pr) :: dBi(size(n)) !! Mole-number gradient of the mixture covolume.
      real(pr) :: dBij(size(n), size(n)) !! Mole-number Hessian of the mixture covolume.
      real(pr) :: D1
      real(pr) :: dD1i(size(n)) !! Mole-number gradient of the cubic denominator parameter.
      real(pr) :: dD1ij(size(n), size(n)) !! Mole-number Hessian of the cubic denominator parameter.
      real(pr) :: L !! Infinite-pressure parameter \(\Lambda\)

      type(NRTLHV) :: ge_model

      real(pr) :: gii(size(n))
      real(pr) :: gji(size(n), size(n))
      real(pr) :: bi(size(n))

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.
      integer :: nc !! Number of components.

      ge_model = self%ge

      nc = size(n)
      bi = self%bi
      call self%Bmix(n, bi, B, dBi, dBij)
      call self%D1Mix(n, self%del1, D1, dD1i, dD1ij)
      call lamdba_hv(nc, D1, L=L)

      gii = - ai/bi * L

      do i=1,nc
         do j=1,nc
            if (self%use_kij(i, j)) then
               ge_model%alpha(i, j) = 0
               ge_model%gji0(i, j) = -2 * sqrt(bi(i) * bi(j)) / (bi(i) + bi(j)) * &
                  sqrt(gii(i) * gii(j)) * (1 - self%kij(i, j)) - gii(j)
               ge_model%gjiT(i, j) = 0
            end if
         end do
      end do

      call ge_model%excess_gibbs( &
         n, T, Ge=Ge, GeT=GeT, GeT2=GeT2, Gen=Gen, GeTn=GeTn, Gen2=Gen2 &
         )

      call DMixHV(n, T,&
         bi, B, dBi, dBij, &
         D1, dD1i, dD1ij, &
         ai, daidt, daidt2, &
         Ge, GeT, GeT2, Gen, GeTn, Gen2,&
         D, dDdT, dDdT2, dDi, dDidT, dDij)
   end subroutine DmixHVNRTL

   ! ===========================================================================
   ! Michelsen's Modified Huron-Vidal 1
   ! ---------------------------------------------------------------------------
   type(MHV) function init_mhv(ge, b, q, lij) result(mixrule)
      !! # MHV
      !! Construct a Michelsen modified Huron-Vidal mixing rule.
      !!
      !! # Examples
      !! ```fortran
      !! mixrule = MHV(ge=ge_model, b=model%b, q=-0.593_pr)
      !! ```
      class(GeModel), intent(in) :: Ge !! Excess Gibbs energy model
      real(pr), intent(in) :: b(:) !! Pure-component covolumes [L/mol]
      real(pr), intent(in) :: q !! MHV \(q\) parameter
      real(pr), optional, intent(in) :: lij(:, :) !! Binary repulsive interaction parameters (default zero)

      integer :: i !! Loop/index variable.
      integer :: nc !! Number of components.

      nc = size(b)

      mixrule%q = q
      mixrule%bi = b
      mixrule%Ge = ge
      if (present(lij)) then
         mixrule%l = lij
      else
         mixrule%l = reshape([(0, i=1, nc**2)], [nc, nc])
      end if
   end function init_mhv

   subroutine BmixMHV(self, n, bi, B, dBi, dBij)
      !! # Repulsive parameter \(B\) mixing rule
      !! Quadratic mixing rule for the repulsive parameter, using
      !! \( b_{ij} = \frac{b_i + b_j}{2} (1 - l_{ij}) \) as a combining rule.
      !!
      !! # Description
      !! Michelsen's modified Huron-Vidal mixing rule assumes a linear mix of
      !! the repulsive parameter.
      !!
      !! \[B = \sum_i n_i b_i\]
      !!
      !! In this implementation the most known crossed combining rule is used:
      !! \[n_T B = \sum_i \sum_j n_i n_j \frac{b_i + b_j}{2} (1 - l_{ij})\]
      !! to provide versatility to the used model.
      !!
      !! @warning
      !! This mixing rule is intended to use only with a linear combining
      !! rule, using \(l_{ij}\) could negatively affect the thermodynamic
      !! consistency of the model.
      !! @endwarning
      !!
      !! # Examples
      !! Called through `CubicEoS%mixrule%Bmix` once the rule is installed.
      use yaeos__models_ar_cubic_mixing_base, only: bmix_linear
      class(MHV), intent(in) :: self !! Mixing rule object.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol]
      real(pr), intent(in) :: bi(:) !! Pure-component covolumes [L/mol]
      real(pr), intent(out) :: B !! Mixture covolume \(B\) [L]
      real(pr), intent(out) :: dBi(:) !! \(\partial B/\partial n_i\)
      real(pr), intent(out) :: dBij(:, :) !! \(\partial^2 B/\partial n_i \partial n_j\)
      call bmix_qmr(n, bi, self%l, b, dbi, dbij)
      ! call bmix_linear(n, bi, b, dbi, dbij)
   end subroutine BmixMHV

   subroutine DmixMHV(self, n, V, T, &
      ai, daidt, daidt2, &
      D, &
      dDdV, dDdT, dDdV2, dDdT2, dDi, dDdTV, dDidV, dDidT, dDij &
      )
      !! # Michelsen Modified Huron-Vidal mixing rule.
      !! Mixing rule at the approximate zero-pressure limit as defined in the
      !! book of Michelsen and Møllerup.
      !!
      !! # Description
      !! At the approximate zero-pressure limit of a cubic equation of state it is
      !! possible to relate the mixing rule for the attractive term with an excess Gibbs energy
      !! model like NRTL with the expression:
      !!
      !! \[
      !! \frac{D}{RTB}(n, T) = \sum_i n_i \frac{a_i(T)}{RT b_i} + \frac{1}{q}
      !!  \left(\frac{G^E(n, T)}{RT} + \sum_i n_i \ln \frac{B}{n_T b_i} \right)
      !! \]
      !!
      !! # References
      !! Michelsen, M. L., *Fluid Phase Equilibria* 60 (1990) 213-219.
      class(MHV), intent(in) :: self !! Mixing rule object.
      real(pr), intent(in) :: V !! Volume [L] (unused)
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: n(:) !! Moles vector [mol]
      real(pr), intent(in) :: ai(:) !! Pure components attractive parameters \(a_i\)
      real(pr), intent(in) :: daidt(:) !! \(\frac{da_i}{dT}\)
      real(pr), intent(in) :: daidt2(:) !! \(\frac{d^2a_i}{dT^2}\)

      real(pr), intent(out) :: D !! Mixture attractive parameter \(n^2a_{mix}\)
      real(pr), intent(out) :: dDdV !! \(\frac{dD}{dV}\) (not assigned)
      real(pr), intent(out) :: dDdT !! \(\frac{dD}{dT}\)
      real(pr), intent(out) :: dDdT2 !! \(\frac{d^2D}{dT^2}\)
      real(pr), intent(out) :: dDdV2 !! \(\frac{d^2D}{dV^2}\)
      real(pr), intent(out) :: dDdTV !! \(\frac{d^2D}{dTV}\)
      real(pr), intent(out) :: dDi(:) !! \(\frac{dD}{dn_i}\)
      real(pr), intent(out) :: dDidV(:) !! \(\frac{d^2D}{dVn_i}\)
      real(pr), intent(out) :: dDidT(:) !! \(\frac{d^2D}{dTn_i}\)
      real(pr), intent(out) :: dDij(:, :) !! \(\frac{d^2D}{dn_{ij}}\)
      real(pr) :: f !! Residual or scalar equation value.
      real(pr) :: fdt
      real(pr) :: fdt2 !! Second-order intermediate quantity used by this procedure.
      real(pr) :: fdi(size(n))
      real(pr) :: fdit(size(n))
      real(pr) :: fdij(size(n), size(n))

      real(pr) :: b
      real(pr) :: bi(size(n))
      real(pr) :: dbi(size(n))
      real(pr) :: dbij(size(n), size(n))
      real(pr) :: Ge !! Excess Gibbs energy.
      real(pr) :: GeT !! Temperature derivative of excess Gibbs energy.
      real(pr) :: GeT2 !! Second temperature derivative of excess Gibbs energy.
      real(pr) :: Gen(size(n)) !! Mole-number gradient of excess Gibbs energy.
      real(pr) :: GeTn(size(n)) !! Mixed temperature/mole-number derivative of excess Gibbs energy.
      real(pr) :: Gen2(size(n), size(n)) !! Mole-number Hessian of excess Gibbs energy.

      real(pr) :: totn !! Total number of moles
      real(pr) :: dot_n_logB_nbi
      real(pr) :: logB_nbi(size(n)) !! \(\ln \frac{B}{n b_i}\)
      real(pr) :: dlogBi_nbi(size(n))
      real(pr) :: d2logBi_nbi(size(n), size(n))

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.
      integer :: l !! Loop index `l`.
      integer :: nc !! Number of components.
      real(pr) :: q

      nc = size(n)
      totn = sum(n)

      q = self%q
      bi = self%bi

      if (.not. solving_volume) then
         call self%ge%excess_gibbs( &
            n, T, Ge=Ge, GeT=GeT, GeT2=GeT2, Gen=Gen, GeTn=GeTn, Gen2=Gen2 &
            )
      else
         call self%ge%excess_gibbs( &
            n, T, Ge=Ge &
            )
      end if
      call self%Bmix(n, bi, B, dBi, dBij)
      logb_nbi = log(B/(totn*bi))
      dot_n_logB_nbi = dot_product(n, logB_nbi)

      do i = 1, nc
         dlogBi_nbi(i) = logB_nbi(i) + sum(n*dBi(i))/B - 1
      end do

      if (.not. solving_volume) then
         do i = 1, nc
            do j = 1, nc
               !TODO: Need to figure out this derivative
               d2logBi_nbi(i, j) = dlogBi_nbi(j) &
                  + (sum(n*dBij(i, j)) + dBi(i))/B &
                  - totn*dBi(i)*dBi(j)/B**2
            end do
         end do

         autodiff: block
            !! Autodiff injection until we can decipher this derivative
            use hyperdual_mod
            type(hyperdual) :: hB
            type(hyperdual) :: hdot_ln_B_nbi
            type(hyperdual) :: hn(nc)

            integer :: ii
            integer :: jj
            hn = n

            do i = 1, nc
               do j = i, nc
                  hn = n
                  hn(i)%f1 = 1
                  hn(j)%f2 = 1

                  hB = 0._pr
                  do ii=1,nc
                     do jj=1,nc
                        hB = hB &
                           + (hn(ii)*hn(jj)) &
                           * 0.5_pr * (bi(ii) + bi(jj)) * (1._pr - self%l(ii, jj))
                     end do
                  end do
                  hB = hB/sum(hn)

                  hdot_ln_B_nbi = sum(hn*log(hB/(sum(hn)*bi)))

                  d2logBi_nbi(i, j) = hdot_ln_B_nbi%f12
                  d2logBi_nbi(j, i) = hdot_ln_B_nbi%f12
               end do
            end do
         end block autodiff
      end if

      f = sum(n*ai/bi) + (Ge + R*T*dot_n_logB_nbi)/q
      fdt = sum(n*daidt/bi) + (GeT + R*dot_n_logB_nbi)/q
      fdt2 = sum(n*daidt2/bi) + (GeT2)/q

      fdi = ai/bi + (1._pr/q)*(GeN + R*T*(dlogBi_nbi))
      fdit = daidt/bi + (1._pr/q)*(GeTn + R*(dlogBi_nbi))

      do i = 1, nc
         do j = 1, nc
            fdij(i, j) = R*T*(d2logBi_nbi(i, j))
            fdij(i, j) = 1/q*(fdij(i, j) + GeN2(i, j))
            fdij(i, j) = &
               dBi(j)*fdi(i) + B*fdij(i, j) + fdi(j)*dBi(i) + f*dBij(i, j)
         end do
      end do

      dDi = B*fdi + f*dBi
      dDidT = B*fdiT + fdT*dBi

      D = f*B
      dDdT = fdT*B
      dDdT2 = fdT2*B
      dDij = fdij

   end subroutine DmixMHV

   subroutine D1MixMHV(self, n, d1i, D1, dD1i, dD1ij)
      !! # D1MixMHV
      !! Calculate the MHV mixture \(\delta_1\) parameter and mole derivatives.
      !!
      !! # Description
      !! Delegates to the RKPR mole-weighted mixing relation for the supplied
      !! component values.
      use yaeos__models_ar_cubic_mixing_base, only: d1mix_rkpr
      class(MHV), intent(in) :: self !! Mixing rule object.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol]
      real(pr), intent(in) :: d1i(:) !! Pure-component \(\delta_{1i}\) values
      real(pr), intent(out) :: D1 !! Mixture \(\delta_1\)
      real(pr), intent(out) :: dD1i(:) !! \(\partial \delta_1/\partial n_i\)
      real(pr), intent(out) :: dD1ij(:, :) !! \(\partial^2 \delta_1/\partial n_i \partial n_j\)
      call d1mix_rkpr(n, d1i, D1, dD1i, dD1ij)
   end subroutine D1MixMHV
end module yaeos__models_cubic_mixing_rules_huron_vidal
