module yaeos__models_cubic_mixing_rules_huron_vidal
   !! # Huron-Vidal (like) mixing rules module
   !! This module contains the mixing rules that are based/similar to the
   !! mixing rules defined by Huron-Vidal
   !!
   !! # Description
   !! Huron-Vidal presented a way to link a \(G^E\) model with a Cubic EoS
   !! mixing rule. This makes it possible to make good predictions on
   !! polar compounds containing mixtures.
   !!
   !! # Examples
   !!
   !! ```fortran
   !!  A basic code example
   !! ```
   !!
   !! # References
   !!
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
      class(GeModel), allocatable :: ge
      real(pr), allocatable :: del1(:)
      real(pr), allocatable :: bi(:)
   contains
      procedure :: Bmix => BmixHV
      procedure :: D1Mix => D1MixHV
      procedure :: Dmix => DmixHV
   end type HV

   type, extends(CubicMixRule) :: MHV
      !! # Michelsen's modified Huron-Vidal mixing rule
      !! Mixing rule at zero-pressure which allows for the inclusion of an
      !! excess-gibbs model.
      !!
      !! # Description
      !! This mixing rule is based on the aproximate zero-pressure limit
      !!  of a cubic equation of state. At the aproximate zero-pressure limit the
      !! attractive parameter can be expressed as:
      !!
      !! \[
      !! \frac{D}{RTB}(n, T) = \sum_i n_i \frac{a_i(T)}{b_i} + \frac{1}{q}
      !!  \left(\frac{G^E(n, T)}{RT} + \sum_i n_i \ln \frac{B}{nb_i} \right)
      !! \]
      !! Where \(q\) is a weak function of temperature. In the case of `MHV`
      !! and simplicity it is considered that depends on the model used.
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
      !! mixrule = MHV(ge=ge_model, q=-0.593_pr, bi=model%b)
      !!
      !! ! Replace the original mixrule on the previously defined model
      !! model%mixrule = mixrule
      !!
      !! ! Ready to do calculations
      !! call pressure(model, n, v, T)
      !! ```
      !!
      !! # References
      !!
      real(pr), allocatable :: l(:, :)
      real(pr), private, allocatable :: bi(:)
      real(pr), private, allocatable :: B, dBi(:), dBij(:, :)
      class(GeModel), allocatable :: ge
      real(pr) :: q
   contains
      procedure :: Bmix => BmixMHV
      procedure :: D1Mix => D1MixMHV
      procedure :: Dmix => DmixMHV
   end type MHV

   type, extends(CubicMixRule) :: HV_NRTL
      !! # HV_NRTL
      !! Huron-Vidal mixing rule including the NRTL model modified by Huron
      !! and Vidal.
      !!
      !! # Description
      !! This is the Huron-Vidal mixing rule that includes the NRTL model
      !! modified by Huron and Vidal. It is a mixing rule that allows to
      !! use the NRTL model as an excess Gibbs energy model and can. be
      !! simplified to the classic Quatratic mixing rules when the parameters
      !! are set to:
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
      !!           \left(1 - k_{ij})\right)
      !! \]
      !!
      !! # Examples
      !!
      !! # References
      type(NRTLHV) :: ge
      real(pr), allocatable :: del1(:)
      real(pr), allocatable :: bi(:)
      logical, allocatable :: use_kij(:, :)
      real(pr), allocatable :: kij(:, :)
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
      !! Quadratinc mixing rule for the repulsive parameter.
      !!
      !! # Description
      !! \[B = \sum_i n_i b_i\]
      use yaeos__models_ar_cubic_mixing_base, only: bmix_linear
      class(HV), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: bi(:)
      real(pr), intent(out) :: B, dBi(:), dBij(:, :)
      call bmix_linear(n, bi, b, dbi, dbij)
   end subroutine BmixHV

   subroutine D1MixHV(self, n, d1i, D1, dD1i, dD1ij)
      use yaeos__models_ar_cubic_mixing_base, only: d1mix_rkpr
      class(HV), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: d1i(:)
      real(pr), intent(out) :: D1
      real(pr), intent(out) :: dD1i(:)
      real(pr), intent(out) :: dD1ij(:, :)
      call d1mix_rkpr(n, d1i, D1, dD1i, dD1ij)
   end subroutine D1MixHV

   subroutine DmixHV(self, n, T, &
      ai, daidt, daidt2, &
      D, dDdT, dDdT2, dDi, dDidT, dDij &
      )
      use yaeos__models_ar_cubic_mixing_base, only: lamdba_hv, mix => DMixHV
      class(HV), intent(in) :: self
      real(pr), intent(in) :: T, n(:)
      real(pr), intent(in) :: ai(:), daidt(:), daidt2(:)
      real(pr), intent(out) :: D, dDdT, dDdT2, dDi(:), dDidT(:), dDij(:, :)

      real(pr) :: b, bi(size(n)), dbi(size(n)), dbij(size(n), size(n))
      real(pr) :: del1(size(n)), del2(size(n))
      real(pr) :: d1, d1i(size(n)), dd1i(size(n)), dd1ij(size(n), size(n))
      real(pr) :: Ge, GeT, GeT2, Gen(size(n)), GeTn(size(n)), Gen2(size(n), size(n))

      real(pr) :: totn !! Total number of moles

      integer :: i, j, nc
      real(pr) :: L, dL(size(n)), dL2(size(n), size(n))

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
      !! use the NRTL model as an excess Gibbs energy model and can. be
      !! simplified to the classic Quatratic mixing rules when the parameters
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
      !!           \left(1 - k_{ij})\right)
      !! \]
      !!
      !! # Examples
      !!
      use yaeos__models_ge_nrtlhv, only: NRTLHV
      real(pr), intent(in) :: b(:)
      real(pr), intent(in) :: del1(:)
      real(pr), intent(in) :: alpha(:, :)
      real(pr), intent(in) :: gji0(:, :)
      real(pr), intent(in) :: gjiT(:, :)
      logical, intent(in) :: use_kij(:, :)
      real(pr), intent(in) :: kij(:, :)

      integer :: i, nc

      nc = size(b)

      mixrule%ge = NRTLHV(b=b, alpha=alpha, gji0=gji0, gjiT=gjiT)
      mixrule%bi = b
      mixrule%del1 = del1
      mixrule%use_kij = use_kij
      mixrule%kij = kij
   end function init_hvnrtl

   subroutine BmixHVNRTL(self, n, bi, B, dBi, dBij)
      !! # Repulsive parameter \(B\) mixing rule
      !! Quadratinc mixing rule for the repulsive parameter.
      !!
      !! # Description
      !! \[B = \sum_i n_i b_i\]
      use yaeos__models_ar_cubic_mixing_base, only: bmix_linear
      class(HV_NRTL), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: bi(:)
      real(pr), intent(out) :: B, dBi(:), dBij(:, :)
      call bmix_linear(n, bi, b, dbi, dbij)
   end subroutine BmixHVNRTL

   subroutine D1MixHVNRTL(self, n, d1i, D1, dD1i, dD1ij)
      use yaeos__models_ar_cubic_mixing_base, only: d1mix_rkpr
      class(HV_NRTL), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: d1i(:)
      real(pr), intent(out) :: D1
      real(pr), intent(out) :: dD1i(:)
      real(pr), intent(out) :: dD1ij(:, :)
      call d1mix_rkpr(n, d1i, D1, dD1i, dD1ij)
   end subroutine D1MixHVNRTL

   subroutine DmixHVNRTL(self, n, T, &
      ai, daidt, daidt2, &
      D, dDdT, dDdT2, dDi, dDidT, dDij &
      )
      use yaeos__models_ar_cubic_mixing_base, only: lamdba_hv, DmixHV
      use yaeos__models_ge_nrtlhv, only: NRTLHV
      class(HV_NRTL), intent(in) :: self
      real(pr), intent(in) :: T, n(:)
      real(pr), intent(in) :: ai(:), daidt(:), daidt2(:)
      real(pr), intent(out) :: D, dDdT, dDdT2, dDi(:), dDidT(:), dDij(:, :)

      real(pr) :: Ge, GeT, GeT2
      real(pr) :: Gen(size(n)), GeTn(size(n)), Gen2(size(n), size(n))

      real(pr) :: B, dBi(size(n)), dBij(size(n), size(n))
      real(pr) :: D1, dD1i(size(n)), dD1ij(size(n), size(n))
      real(pr) :: L

      type(NRTLHV) :: ge_model

      real(pr) :: gii(size(n)), gji(size(n), size(n))
      real(pr) :: bi(size(n))

      integer :: i, j, nc

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
      class(GeModel), intent(in) :: Ge
      real(pr), intent(in) :: b(:)
      real(pr), intent(in) :: q
      real(pr), optional, intent(in) :: lij(:, :)

      integer :: i, nc

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
      !! Quadratinc mixing rule for the repulsive parameter, using
      !! \( b_{ij} = \frac{b_i + b_j}{2} (1 - l_{ij}) \) as a combining rule.
      !!
      !! # Description
      !! Michelsen's modified Huron-Vidal mixing rule assumes a linear mix of
      !! the repulsive parameter.
      !!
      !! \[B = \sum_i n_i b_i\]
      !!
      !! In this implementation the most known crossed combining rule is used:
      !! \[nB = \sum_i \sum_j \frac{b_i + b_j}{2} (1 - l_{ij})\]
      !! to provide versatility to the used model.
      !!
      !! @warning
      !! This mixing rule is intended to use only with a linear combining
      !! rule, using \(l_{ij}\) could negatively affect the thermodynamic
      !! consistency of the model.
      !! @endwarning
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  A basic code example
      !! ```
      !!
      !! # References
      !!
      use yaeos__models_ar_cubic_mixing_base, only: bmix_linear
      class(MHV), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: bi(:)
      real(pr), intent(out) :: B, dBi(:), dBij(:, :)
      call bmix_qmr(n, bi, self%l, b, dbi, dbij)
      ! call bmix_linear(n, bi, b, dbi, dbij)
   end subroutine BmixMHV

   subroutine DmixMHV(self, n, T, &
      ai, daidt, daidt2, &
      D, dDdT, dDdT2, dDi, dDidT, dDij &
      )
      !! # Michelsen Modified Huron-Vidal mixing rule.
      !! Mixing rule at infinite pressure as defined in the book of Michelsen and
      !! Møllerup.
      !!
      !! # Description
      !! At the infinite pressure limit of a cubic equation of state it is possible to
      !! relate teh mixing rule for the attractive term with a excess Gibbs energy
      !! model like NRTL with the expression:
      !!
      !! \[
      !! \frac{D}{RTB}(n, T) = \sum_i n_i \frac{a_i(T)}{b_i} + \frac{1}{q}
      !!  \left(\frac{G^E(n, T)}{RT} + \sum_i n_i \ln \frac{B}{nb_i} \right)
      !! \]
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  type(CubicEoS)
      !! ```
      !!
      !! # References
      !!
      class(MHV), intent(in) :: self
      real(pr), intent(in) :: T, n(:)
      real(pr), intent(in) :: ai(:), daidt(:), daidt2(:)
      real(pr), intent(out) :: D, dDdT, dDdT2, dDi(:), dDidT(:), dDij(:, :)
      real(pr) :: f, fdt, fdt2, fdi(size(n)), fdit(size(n)), fdij(size(n), size(n))

      real(pr) :: b, bi(size(n)), dbi(size(n)), dbij(size(n), size(n))
      real(pr) :: Ge, GeT, GeT2, Gen(size(n)), GeTn(size(n)), Gen2(size(n), size(n))

      real(pr) :: totn !! Total number of moles
      real(pr) :: dot_n_logB_nbi
      real(pr) :: logB_nbi(size(n)) !! \(\ln \frac{B}{n b_i}\)
      real(pr) :: dlogBi_nbi(size(n))
      real(pr) :: d2logBi_nbi(size(n), size(n))

      integer :: i, j, l, nc
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

            integer :: ii, jj
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
      use yaeos__models_ar_cubic_mixing_base, only: d1mix_rkpr
      class(MHV), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: d1i(:)
      real(pr), intent(out) :: D1
      real(pr), intent(out) :: dD1i(:)
      real(pr), intent(out) :: dD1ij(:, :)
      call d1mix_rkpr(n, d1i, D1, dD1i, dD1ij)
   end subroutine D1MixMHV
end module yaeos__models_cubic_mixing_rules_huron_vidal
