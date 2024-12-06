module yaeos__models_ge
   !! Excess Gibbs Models.
   use yaeos__constants, only: pr, R
   use yaeos__models_base, only: BaseModel
   implicit none

   type, extends(BaseModel), abstract :: GeModel
      !! Excess Gibbs energy model.
   contains
      procedure(excess_gibbs), deferred :: excess_gibbs
      procedure :: ln_activity_coefficient => ln_activity_coefficient
      procedure :: excess_enthalpy => excess_enthalpy
      procedure :: excess_entropy => excess_entropy
   end type

   abstract interface
      subroutine excess_gibbs(self, n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)
         !! Calculate Excess Gibbs and its derivatives.
         !!
         import pr, GeModel
         class(GeModel), intent(in) :: self !! Model
         real(pr), intent(in) ::n(:) !! Moles vector
         real(pr), intent(in) :: T !! Temperature [K]
         real(pr), optional, intent(out) :: Ge !! Excess Gibbs free energy
         real(pr), optional, intent(out) :: GeT !! \(\frac{dG^E}{dT}\)
         real(pr), optional, intent(out) :: GeT2 !! \(\frac{d^2G^E}{dT^2}\)
         real(pr), optional, intent(out) :: Gen(size(n)) !! \(\frac{dG^E}{dn}\)
         real(pr), optional, intent(out) :: GeTn(size(n)) !! \(\frac{d^2G^E}{dTdn}\)
         real(pr), optional, intent(out) :: Gen2(size(n), size(n)) !! \(\frac{d^2G^E}{dn^2}\)
      end subroutine
   end interface

contains

   subroutine ln_activity_coefficient(self, n, T, lngamma, dlngammadT, dlngammadn)
      !! Calculate natural logarithm of activity coefficients.
      !!
      !! \[
      !!  \ln \gamma_i = \frac{1}{RT} \frac{\partial G^E}{\partial n_i}
      !! \]
      !!
      class(GeModel), intent(in) :: self !! Model
      real(pr), intent(in) :: n(:) !! Moles vector
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), optional, intent(out) :: lngamma(:)
      !! Natural logarithm of activity coefficients
      real(pr), optional, intent(out) :: dlngammadT(size(n))
      !! \(\frac{d\ln \gamma_i}{dT}\)
      real(pr), optional, intent(out) :: dlngammadn(size(n),size(n))
      !! \(\frac{d\ln \gamma_i}{dn_j}\)

      real(pr) :: Ge, Gen(size(n)), GeTn(size(n)), Gen2(size(n), size(n))

      ! little if present optmization
      if (.not. present(dlngammadn)) then
         if (.not. present(dlngammadT)) then
            call self%excess_gibbs(n, T, Ge=Ge, Gen=Gen)
         else
            call self%excess_gibbs(n, T, Ge=Ge, Gen=Gen, GeTn=GeTn)
         end if
      else
         call self%excess_gibbs(n, T, Ge=Ge, Gen=Gen, GeTn=GeTn, Gen2=Gen2)
      end if

      if (present(lngamma)) lngamma = GeN / (R * T)
      if (present(dlngammadT)) dlngammadT = (GeTn - GeN / T) / (R * T)
      if (present(dlngammadn)) dlngammadn = Gen2 / (R * T)
   end subroutine

   subroutine excess_enthalpy(self, n, T, He, HeT, Hen)
      !! Calculate Excess enthalpy and its derivatives.
      !!
      !! \[
      !! H^E = G^E - T \frac{\partial G^E}{\partial T}
      !! \]
      !!
      !! ## References
      !! [1] https://en.wikipedia.org/wiki/Gibbs%E2%80%93Helmholtz_equation
      !!
      class(GeModel), intent(in) :: self !! Model
      real(pr), intent(in) :: n(:) !! Moles vector
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), optional, intent(out) :: He !! Excess enthalpy
      real(pr), optional, intent(out) :: HeT !! \(\frac{dH^E}{dT}\)
      real(pr), optional, intent(out) :: Hen(:) !! \(\frac{dH^E}{dn}\)

      real(pr) :: Ge, GeT, GeT2, Gen(size(n)), GeTn(size(n))

      call self%excess_gibbs(&
         n, T, Ge=Ge, GeT=GeT, GeT2=GeT2, Gen=Gen, GeTn=GeTn &
         )

      if (present(He)) He = Ge - T*GeT
      if (present(HeT)) HeT = -T * GeT2
      if (present(Hen)) Hen = Gen - T*GeTn
   end subroutine excess_enthalpy

   subroutine excess_entropy(self, n, T, Se, SeT, Sen)
      !! Calculate Excess entropy and its derivatives.
      !!
      !! \[
      !! S^E = \frac{H^E - G^E}{T}
      !! \]
      !!
      class(GeModel), intent(in) :: self !! Model
      real(pr), intent(in) :: n(:) !! Moles vector
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), optional, intent(out) :: Se !! Excess entropy
      real(pr), optional, intent(out) :: SeT !! \(\frac{dS^E}{dT}\)
      real(pr), optional, intent(out) :: Sen(:) !! \(\frac{dS^E}{dn}\)

      real(pr) :: Ge, GeT, Gen(size(n))
      real(pr) :: He, HeT, Hen(size(n))

      call self%excess_gibbs(n, T, Ge=Ge, GeT=GeT, Gen=Gen)
      call self%excess_enthalpy(n, T, He=He, HeT=HeT, Hen=Hen)

      if (present(Se)) Se = (He - Ge) / T
      if (present(SeT)) SeT = ((HeT - GeT) * T - He + Ge) / T**2
      if (present(Sen)) Sen = (Hen - Gen) / T
   end subroutine excess_entropy
end module
