module yaeos__models_ar_volume_translations_constant
   !! # Constant volume translation
   !! Temperature-independent component shifts with a linear combining rule.
   !!
   !! # Description
   !! `ConstantVolumeTranslation` wraps an `ArModel` with
   !! \(C(n)=\sum_i n_i c_i\). The inherited `VolumeTranslatedModel` supplies
   !! the translated Helmholtz energy, properties, and volume solver.
   !!
   !! # Examples
   !! Constructor fragment, after defining the base EOS and shifts in L/mol:
   !!
   !! ```fortran
   !! model = ConstantVolumeTranslation(base, c)
   !! ```
   !!
   !! # References
   !! Péneloux, A., Rauzy, E., & Fréze, R. (1982). A consistent correction for
   !! Redlich-Kwong-Soave volumes. Fluid phase equilibria, 8(1), 7-23.
   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel
   use yaeos__models_ar_volume_translated, only: VolumeTranslatedModel
   implicit none
   private

   public :: ConstantVolumeTranslation

   type, extends(VolumeTranslatedModel) :: ConstantVolumeTranslation
      !! # ConstantVolumeTranslation
      !! An EOS translated by constant molar shifts for each component.
      !!
      !! # Description
      !! The coefficients have the base model's component order and units
      !! of L/mol. With physical volume \(V\), the base EOS volume is evaluated
      !! at \(W=V+\sum_i n_i c_i\); positive shifts reduce the physical volume
      !! at fixed pressure and temperature. The constructor copies the base
      !! model, its component metadata, and the shift coefficients.
      !!
      !! # Examples
      !! Constructor fragment with illustrative, not fitted, shifts:
      !!
      !! ```fortran
      !! model = ConstantVolumeTranslation(base, c=[0.01_pr, 0.02_pr])
      !! ```
      real(pr), allocatable :: c(:)
      !! Temperature-independent molar shifts [L/mol], one per component.
   contains
      procedure :: translation => constant_translation
   end type ConstantVolumeTranslation

   interface ConstantVolumeTranslation
      module procedure :: constant_volume_translation
   end interface ConstantVolumeTranslation

contains

   function constant_volume_translation(base, c) result(model)
      !! # constant_volume_translation
      !! Construct a volume-translated EOS from a base model and molar shifts.
      !!
      !! # Description
      !! `c` must contain one coefficient per component in `base`, in the
      !! same order. Coefficients are supplied directly, not estimated from
      !! critical properties. Zero and negative coefficients are allowed.
      !!
      !! # Examples
      !! Constructor fragment:
      !! ```fortran
      !! model = ConstantVolumeTranslation(base, c)
      !! ```
      class(ArModel), intent(in) :: base
      !! Untranslated EOS to copy.
      real(pr), intent(in) :: c(:)
      !! Constant component shifts in the base model's order [L/mol].
      type(ConstantVolumeTranslation) :: model
      !! Initialized EOS with the inherited residual-property interface.

      call model%init(base)
      model%c = c
   end function constant_volume_translation

   subroutine constant_translation(self, n, T, C, CT, CT2, Cn, CTn, Cn2)
      !! # constant_translation
      !! Evaluate the linear extensive shift and its requested derivatives.
      !!
      !! # Description
      !! \(C=\sum_i n_i c_i\) and \(C_{n_i}=c_i\). Every temperature
      !! derivative and the amount Hessian are zero. Amounts need not sum
      !! to one; `C` is a total volume, not a molar volume.
      !!
      !! # Examples
      !! Evaluation fragment:
      !! ```fortran
      !! call model%translation(n, T, C, Cn=Cn)
      !! ```
      class(ConstantVolumeTranslation), intent(in) :: self
      !! Model containing the component shifts.
      real(pr), intent(in) :: n(:)
      !! Component mole amounts [mol].
      real(pr), intent(in) :: T
      !! Temperature [K]; unused because the shifts are constant.
      real(pr), intent(out) :: C
      !! Extensive volume translation [L].
      real(pr), optional, intent(out) :: CT
      !! Zero temperature derivative [L/K].
      real(pr), optional, intent(out) :: CT2
      !! Zero second temperature derivative [L/K^2].
      real(pr), optional, intent(out) :: Cn(size(n))
      !! Component shifts, the amount gradient [L/mol].
      real(pr), optional, intent(out) :: CTn(size(n))
      !! Zero temperature/amount derivatives [L/(K mol)].
      real(pr), optional, intent(out) :: Cn2(size(n), size(n))
      !! Zero amount Hessian [L/mol^2].

      C = dot_product(n, self%c)
      if (present(CT)) CT = 0._pr
      if (present(CT2)) CT2 = 0._pr
      if (present(Cn)) Cn = self%c
      if (present(CTn)) CTn = 0._pr
      if (present(Cn2)) Cn2 = 0._pr
   end subroutine constant_translation

end module yaeos__models_ar_volume_translations_constant
