module yaeos__models_cubic_mixing_rules_huron_vidal
   use yaeos__constants, only: pr
   use yaeos__models_ar_genericcubic, only: CubicMixRule
   use yaeos__models_ar_cubic_quadratic_mixing, only: &
      D1mix_constant, bmix_lij
   use yaeos__models_ge, only: GeModel
   implicit none

   type, extends(CubicMixRule) :: HV
      real(pr), allocatable :: l(:, :)
      real(pr), private, allocatable :: bi(:)
      real(pr), private, allocatable :: B, dBi(:), dBij(:, :)
      class(GeModel), allocatable :: ge
   contains
      procedure :: Bmix => BmixHV
      procedure :: D1Mix => D1Mix_constantHV
      procedure :: Dmix => DmixHV
   end type

contains

   type(HV) function init(ge, b, lij) result(mixrule)
      class(GeModel), intent(in) :: Ge
      real(pr), intent(in) :: b(:)
      real(pr), intent(in) :: lij(:, :)

      mixrule%bi = b
      mixrule%l = lij
      mixrule%Ge = ge
   end function

   subroutine BmixHV(self, n, bi, B, dBi, dBij)
      class(HV), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: bi(:)
      real(pr), intent(out) :: B, dBi(:), dBij(:, :)
      call bmix_lij(n, bi, self%l, b, dbi, dbij)
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
      call self%Bmix(n, self%bi, B, dBi, dBij)
   end subroutine

end module yaeos__models_cubic_mixing_rules_huron_vidal
