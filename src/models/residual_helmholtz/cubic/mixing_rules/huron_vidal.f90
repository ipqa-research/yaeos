module yaeos__models_cubic_mixing_rules_huron_vidal
   use yaeos__constants, only: pr, R
   use yaeos__models_ar_genericcubic, only: CubicMixRule
   use yaeos__models_ar_cubic_mixing_base, only: bmix_qmr
   use yaeos__models_ge, only: GeModel
   implicit none

   type, extends(CubicMixRule) :: HV
      real(pr), allocatable :: l(:, :)
      real(pr), private, allocatable :: bi(:)
      real(pr), private, allocatable :: B, dBi(:), dBij(:, :)
      class(GeModel), allocatable :: ge
      real(pr) :: q
   contains
      procedure :: Bmix => BmixHV
      procedure :: D1Mix => D1Mix_constantHV
      procedure :: Dmix => DmixHV
   end type

   interface HV
      module procedure :: init
   end interface

contains

   type(HV) function init(ge, b, lij) result(mixrule)
      class(GeModel), intent(in) :: Ge
      real(pr), intent(in) :: b(:)
      real(pr), optional, intent(in) :: lij(:, :)

      integer :: i, nc

      nc = size(b)

      mixrule%bi = b
      mixrule%Ge = ge
      if (present(lij)) then
         mixrule%l = lij
      else
         mixrule%l = reshape([(0, i=1, nc**2)], [nc, nc])
      end if
   end function

   subroutine BmixHV(self, n, bi, B, dBi, dBij)
      class(HV), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: bi(:)
      real(pr), intent(out) :: B, dBi(:), dBij(:, :)
      call bmix_qmr(n, bi, self%l, b, dbi, dbij)
   end subroutine

   subroutine DmixHV(self, n, T, &
                     ai, daidt, daidt2, &
                     D, dDdT, dDdT2, dDi, dDidT, dDij &
                     )
      class(HV), intent(in) :: self
      real(pr), intent(in) :: T, n(:)
      real(pr), intent(in) :: ai(:), daidt(:), daidt2(:)
      real(pr), intent(out) :: D, dDdT, dDdT2, dDi(:), dDidT(:), dDij(:, :)

      real(pr) :: b, bi(size(n)), dbi(size(n)), dbij(size(n), size(n))
      real(pr) :: Ge, GeT, GeT2, Gen(size(n)), GeTn(size(n)), Gen2(size(n), size(n))

      real(pr) :: totn !! Total number of moles
      real(pr) :: logB_nbi(size(n)) !! \(\ln \frac{B}{n b_i}\)
      real(pr) :: dlogBi_nbi(size(n), size(n))

      integer :: i, j, l, nc
      real(pr) :: q

      nc = size(n)
      totn = sum(n)

      q = self%q
      bi = self%bi

      call self%Bmix(n, bi, B, dBi, dBij)
      call self%ge%excess_gibbs( &
         n, T, Ge=Ge, GeT=GeT, GeT2=GeT2, Gen=Gen, GeTn=GeTn, Gen2=Gen2 &
         )

      logb_nbi = log(B/(totn*bi))
      D = (sum(n*ai/bi) + 1/q*(Ge + R*T*sum(n*logB_nbi)))
      dDdT = (sum(n*daidt/bi) + 1/q*(GeT + R*sum(n*logB_nbi)))
      dDdT2 = (sum(n*daidt2/bi) + 1/q*(GeT2 + R*sum(n*logB_nbi)))
      dDi = ai/bi + 1._pr/q*(GeN + R*T*(logB_nbi + sum(n*(dBi/B))))

      do i = 1, nc
         do j = 1, nc
            dDij(i, j) = R*T* (sum(n*(dBij(:, j)/B * (-dBi(j)*dBi(i))/B**2)) + dBi(j)/B)
            dDij(i, j) = 1._pr/q*(dDij(i, j) + GeN2(i, j))
            dDij(i, j) = dDi(i)*dBi(j) + dDij(i, j)*dBi(i) + dBi(i)*dDi(j) + dBij(i, j)*D
         end do
      end do

      dDi = B*dDi + D*dBi
      D = D*B
      dDdT = dDdT*B
      dDdT2 = dDdT2*B
   end subroutine

   subroutine D1Mix_constantHV(self, n, d1i, D1, dD1i, dD1ij)
      class(HV), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: d1i(:)
      real(pr), intent(out) :: D1
      real(pr), intent(out) :: dD1i(:)
      real(pr), intent(out) :: dD1ij(:, :)

      D1 = d1i(1)
      dD1i = 0
      dD1ij = 0
   end subroutine D1Mix_constantHV
end module yaeos__models_cubic_mixing_rules_huron_vidal
