module yaeos_models_ge
   !! Excess Gibbs Models.
   use yaeos_constants, only: pr, R
   use yaeos_models_base, only: BaseModel
   implicit none

   type, extends(BaseModel), abstract :: GeModel
      !! Excess Gibbs energy model.
   contains
      procedure(excess_gibbs), deferred :: excess_gibbs
      procedure :: ln_activity_coefficient => ln_activity_coefficient
   end type

   abstract interface
      subroutine excess_gibbs(self, n, t, Ge, GeT, GeT2, Gen, GeTn, Gen2)
         !! Excess Gibbs and derivs procedure
         import pr, GeModel
         class(GeModel), intent(in) :: self !! Model
         real(pr), intent(in) ::n(:) !! Moles vector
         real(pr), intent(in) :: t !! Temperature [K]
         real(pr), optional, intent(out) :: Ge !! Excess Gibbs
         real(pr), optional, intent(out) :: GeT !! \(\frac{dG^E}{dT}\)
         real(pr), optional, intent(out) :: GeT2 !! \(\frac{d^2G^E}{dT^2}\)
         real(pr), optional, intent(out) :: Gen(size(n))
         real(pr), optional, intent(out) :: GeTn(size(n))
         real(pr), optional, intent(out) :: Gen2(size(n), size(n))
      end subroutine
   end interface

contains

   subroutine ln_activity_coefficient(self, n, T, lngamma)
      class(GeModel), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: T
      real(pr), intent(out) :: lngamma(:)

      real(pr) :: ge, dgedn(size(n))

      call self%excess_gibbs(n, t, ge=ge, gen=dgedn)
      lngamma = dgedn/(R*T)
   end subroutine
end module
