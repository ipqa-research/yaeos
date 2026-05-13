module yaeos__models_ar_cubic_sddlc
   !! # Cubic EoS with segmented-DDLC mixing rule module
   !!
   !! This module contains the implementation of a generic CubicEoS using the
   !! segmented Density Dependent Local Compositions mixing rule.
   use yaeos__constants, only: pr
   use yaeos__autodiff, only: hyperdual, ArModelAdiff

   type, extends(ArModelAdiff) :: CubicEOSsDDLC
      real(pr), allocatable :: ac(:)
      real(pr), allocatable :: b(:)
      real(pr), allocatable :: del1(:)
   contains
      procedure :: Ar => Ar
      procedure :: get_v0 => get_v0
   end type CubicEOSsDDLC

contains

   type(hyperdual) function Ar(self, n, V, T)
      class(CubicEOSsDDLC) :: self !! Model
      type(hyperdual), intent(in) :: n(:) !! Vector of moles
      type(hyperdual), intent(in) :: V !! Volume [L]
      type(hyperdual), intent(in) :: T !! Temperature [K]
   end function Ar

   real(pr) function get_v0(self, n, P, T)
      class(CubicEOSsDDLC) :: self !! Model
      real(pr), intent(in) :: n(:) !! Vector of moles
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
   end function get_v0

end module yaeos__models_ar_cubic_sddlc
