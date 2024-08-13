module yaeos__equilibria_auxiliar
   !! Auxiliar functions used for phase-equilibria calculation.
   use yaeos__constants, only: pr
   use yaeos__models_base, only: BaseModel
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

   real(pr) function P_wilson(model, z, T) result(P)
      class(BaseModel), intent(in) :: model
      real(pr), intent(in) :: z(:)
      real(pr), intent(in) :: T

      P = 1.0_pr/sum(&
         z*model%components%Pc &
         * exp(5.373_pr &
         * (1 + model%components%w)*(1 - model%components%Tc/T)))
   end function P_wilson
end module yaeos__equilibria_auxiliar
