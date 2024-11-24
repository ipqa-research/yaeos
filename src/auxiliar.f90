module yaeos__equilibria_auxiliar
   !! Auxiliar functions used for phase-equilibria calculation.
   use yaeos__constants, only: pr
   use yaeos__models_base, only: BaseModel
   implicit none

contains

   function k_wilson(model, T, P) result(K)
      !! # K_wilson
      !!
      !! ## Description
      !! K-factors regression done by Wilson, used for initialization.
      class(BaseModel), intent(in) :: model
      real(pr), intent(in) :: T
      real(pr), intent(in) :: P
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
      real(pr), intent(in) :: z(:) !! Mole fractions of the components.
      real(pr), intent(in) :: T !! Temperature [K].

      P = 1.0_pr/sum(&
         z*model%components%Pc &
         * exp(5.373_pr &
         * (1 + model%components%w)*(1 - model%components%Tc/T)))
   end function P_wilson
end module yaeos__equilibria_auxiliar
