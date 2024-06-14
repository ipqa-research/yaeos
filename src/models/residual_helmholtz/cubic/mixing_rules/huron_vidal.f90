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
        mixrule%l = reshape([(0, i=1,nc**2)], [nc,nc])
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

      totn = sum(n)
      call self%Bmix(n, self%bi, B, dBi, dBij)
      call self%ge%excess_gibbs(n, T, Ge=Ge, GeT=GeT, GeT2=GeT2, Gen=Gen, GeTn=GeTn, Gen2=Gen2)

      D = sum(n*ai/self%bi) + 1/self%q*(Ge + R*T*sum(n*log(B/(totn*self%bi))))
      D = D*B
   end subroutine

   subroutine D1Mix_constantHV(self, n, d1i, D1, dD1i, dD1ij)
      class(HV), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: d1i(:)
      real(pr), intent(out) :: D1
      real(pr), intent(out) :: dD1i(:)
      real(pr), intent(out) :: dD1ij(:, :)
   end subroutine D1Mix_constantHV
end module yaeos__models_cubic_mixing_rules_huron_vidal
