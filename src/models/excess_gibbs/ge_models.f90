module yaeos__models_ge
   !! # Excess-Gibbs models
   !! Abstract API and thermodynamic property relations for activity models.
   !!
   !! # Description
   !! A `GeModel` implementation supplies the excess Gibbs energy and its
   !! temperature/composition derivatives. The base type derives activity
   !! coefficients and excess enthalpy, entropy, and heat capacity from those
   !! derivatives. Energies are extensive and expressed in the units of
   !! \(nRT\), i.e. bar L with the library's `R` [bar L / (mol K)].
   !!
   !! # Examples
   !! ```fortran
   !! model = NRTL(A, B, alpha)
   !! call model%ln_activity_coefficient(n, T, ln_gamma)
   !! ```
   use yaeos__constants, only: pr, R
   use yaeos__models_base, only: BaseModel
   implicit none

   type, extends(BaseModel), abstract :: GeModel
      !! # GeModel
      !! Abstract base for models defined by an excess Gibbs energy.
      !!
      !! # Description
      !! Derived types implement the `excess_gibbs` deferred procedure. The
      !! inherited methods calculate activity coefficients and excess
      !! thermodynamic properties using the Gibbs-energy identities.
      !!
      !! # Examples
      !! See `NRTL`, `UNIQUAC`, and `UNIFAC` for concrete models.
   contains
      procedure(excess_gibbs), deferred :: excess_gibbs
      procedure :: ln_activity_coefficient => ln_activity_coefficient
      procedure :: excess_enthalpy => excess_enthalpy
      procedure :: excess_entropy => excess_entropy
      procedure :: excess_Cp => excess_Cp
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
         real(pr), optional, intent(out) :: GeTn(size(n))
         !! \(\frac{d^2G^E}{dTdn}\)
         real(pr), optional, intent(out) :: Gen2(size(n), size(n))
         !! \(\frac{d^2G^E}{dn^2}\)
      end subroutine
   end interface

contains

   subroutine ln_activity_coefficient(self, n, T, lngamma, dlngammadT, dlngammadn)
      !! # ln_activity_coefficient
      !!
      !! # Description
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

      real(pr) :: Ge
      real(pr) :: Gen(size(n))
      real(pr) :: GeTn(size(n))
      real(pr) :: Gen2(size(n), size(n))

      logical :: tt
      logical :: dt
      logical :: dn

      tt = present(lngamma)
      dt = present(dlngammadT)
      dn = present(dlngammadn)

      if (tt .and. .not. dt .and. .not. dn) then
         call self%excess_gibbs(n, T, Ge=Ge, Gen=Gen)
      else if (.not. dn) then
         call self%excess_gibbs(n, T, Ge=Ge, Gen=Gen, GeTn=GeTn)
      else
         call self%excess_gibbs(n, T, Ge=Ge, Gen=Gen, GeTn=GeTn, Gen2=Gen2)
      end if

      if (tt) lngamma = Gen / (R * T)
      if (dt) dlngammadT = (GeTn - Gen / T) / (R * T)
      if (dn) dlngammadn = Gen2 / (R * T)
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

      real(pr) :: Ge
      real(pr) :: GeT
      real(pr) :: GeT2
      real(pr) :: Gen(size(n))
      real(pr) :: GeTn(size(n))

      call self%excess_gibbs(&
         n, T, Ge=Ge, GeT=GeT, GeT2=GeT2, Gen=Gen, GeTn=GeTn &
         )

      if (present(He)) He = Ge - T*GeT
      if (present(HeT)) HeT = -T * GeT2
      if (present(Hen)) Hen = Gen - T*GeTn
   end subroutine excess_enthalpy

   subroutine excess_entropy(self, n, T, Se, SeT, Sen)
      !! # excess_entropy
      !!
      !! # Description
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

      real(pr) :: Ge
      real(pr) :: GeT
      real(pr) :: GeT2
      real(pr) :: GeTn(size(n))

      call self%excess_gibbs(n, T, Ge=Ge, GeT=GeT, GeT2=GeT2, GeTn=GeTn)

      if (present(Se)) Se = -GeT
      if (present(SeT)) SeT = -GeT2
      if (present(Sen)) Sen = -GeTn
   end subroutine excess_entropy

   subroutine excess_Cp(self, n, T, Cpe)
      !! # excess_Cp
      !!
      !! # Description
      !! Calculate Excess heat capacity.
      !!
      !! \[
      !! C_p^E = -T \frac{\partial^2 G^E}{\partial T^2}
      !! \]
      !!
      class(GeModel), intent(in) :: self !! Model
      real(pr), intent(in) :: n(:) !! Moles vector
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(out) :: Cpe !! Excess heat capacity [bar L / K]

      real(pr) :: Ge
      real(pr) :: GeT
      real(pr) :: GeT2
      real(pr) :: Gen(size(n))
      real(pr) :: GeTn(size(n))

      call self%excess_gibbs(&
         n, T, Ge=Ge, GeT=GeT, GeT2=GeT2, Gen=Gen, GeTn=GeTn &
         )

      Cpe = -T * GeT2
   end subroutine excess_Cp
end module
