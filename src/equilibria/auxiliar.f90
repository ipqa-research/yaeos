module yaeos__equilibria_auxiliar
   !! # Phase-equilibrium initialization utilities
   !! Wilson correlations for initial K-factors and pressure estimates.
   !!
   !! # Description
   !! `k_wilson` provides component volatility ratios for initializing flash
   !! calculations; `p_wilson` estimates the pressure from a composition and
   !! temperature. These are empirical initial estimates, not equilibrium
   !! solutions.
   !!
   !! # Examples
   !! ```fortran
   !! K0 = k_wilson(model, T, P)
   !! P0 = p_wilson(model, z, T)
   !! ```
   use yaeos__constants, only: pr
   use yaeos__models_base, only: BaseModel
   use yaeos__models_ar, only: ArModel
   implicit none

contains

   function k_wilson(model, T, P) result(K)
      !! # K_wilson
      !!
      !! ## Description
      !! K-factors regression done by Wilson, used for initialization.
      !!
      !! # Examples
      !! ```fortran
      !! K0 = k_wilson(model, T, P)
      !! ```
      class(BaseModel), intent(in) :: model !! Model supplying component critical constants.
      real(pr), intent(in) :: T !! Temperature [K].
      real(pr), intent(in) :: P !! Pressure [bar].
      real(pr)  :: K(size(model%components%pc))

      K = (model%components%Pc/P) &
         * exp(5.373_pr*(1 + model%components%w)&
         * (1 - model%components%Tc/T))
   end function k_wilson

   real(pr) function P_wilson(model, z, T) result(P)
      !! # P_wilson
      !!
      !! ## Description
      !! Calculate the pressure at a given T of a mixture using the Wilson
      !! equation.
      class(BaseModel), intent(in) :: model !! Model of the mixture.
      real(pr), intent(in) :: z(:) !! Overall composition (mole fractions).
      real(pr), intent(in) :: T !! Temperature [K].

      P = 1.0_pr/sum(&
         z*model%components%Pc &
         * exp(5.373_pr &
         * (1 + model%components%w)*(1 - model%components%Tc/T)))
   end function P_wilson

end module yaeos__equilibria_auxiliar
