module yaeos__phase_equilibria_auxiliar
   !! Auxiliar functions used for phase-equilibria calculation.
   use yaeos_constants, only: pr
   use yaeos_models_base, only: BaseModel
   implicit none

contains
   function k_wilson(model, T, P) result(K)
      !! K-factors regressión done by Wilson, used for initialization.
      class(BaseModel), intent(in) :: model
      real(pr), intent(in) :: T
      real(pr), intent(in) :: P
      real(pr)  :: K(size(model%components%pc))

      K = (model%components%Pc/P) &
         * exp(5.373_pr*(1 + model%components%w)&
         * (1 - model%components%Tc/T))
   end function k_wilson
end module yaeos__phase_equilibria_auxiliar
