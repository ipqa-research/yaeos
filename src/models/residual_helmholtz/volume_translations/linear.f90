module yaeos__models_ar_volume_translations_linear
   !! # Linear temperature-dependent volume translation
   !! Component molar shifts \(c_i(T)=a_i+b_iT\).
   !!
   !! # Description
   !! Uses \(C(n,T)=\sum_i n_i(a_i+b_iT)\) with the inherited Helmholtz
   !! transformation. No restriction is imposed on the signs of the parameters.
   !!
   !! # References
   !! J. Shi and H. A. Li, Fluid Phase Equilibria (2016),
   !! doi:10.1016/j.fluid.2016.09.017.
   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel
   use yaeos__models_ar_volume_translated, only: VolumeTranslatedModel
   implicit none
   private

   public :: LinearVolumeTranslation

   type, extends(VolumeTranslatedModel) :: LinearVolumeTranslation
      !! Volume translation with a linear temperature-dependent shift per component.
      real(pr), allocatable :: a(:)
      !! Intercepts in component order [L/mol].
      real(pr), allocatable :: b(:)
      !! Temperature slopes in component order [L/(mol K)].
   contains
      procedure :: translation => linear_translation
   end type LinearVolumeTranslation

   interface LinearVolumeTranslation
      module procedure :: linear_volume_translation
   end interface LinearVolumeTranslation

contains

   function linear_volume_translation(base, a, b) result(model)
      !! Construct a translated copy of the base model.
      !! Arrays `a` and `b` must match the base model's component order and size.
      class(ArModel), intent(in) :: base
      !! Base model to copy.
      real(pr), intent(in) :: a(:)
      !! Component intercepts [L/mol].
      real(pr), intent(in) :: b(:)
      !! Component temperature slopes [L/(mol K)].
      type(LinearVolumeTranslation) :: model
      !! Initialized translated model.

      call model%init(base)
      model%a = a
      model%b = b
   end function linear_volume_translation

   subroutine linear_translation(self, n, T, C, CT, CT2, Cn, CTn, Cn2)
      !! Evaluate the extensive shift and its requested derivatives.
      class(LinearVolumeTranslation), intent(in) :: self
      !! Translation model.
      real(pr), intent(in) :: n(:)
      !! Component amounts [mol].
      real(pr), intent(in) :: T
      !! Temperature [K].
      real(pr), intent(out) :: C
      !! Extensive shift [L].
      real(pr), optional, intent(out) :: CT
      !! Temperature derivative [L/K].
      real(pr), optional, intent(out) :: CT2
      !! Second temperature derivative [L/K²].
      real(pr), optional, intent(out) :: Cn(size(n))
      !! Amount gradient [L/mol].
      real(pr), optional, intent(out) :: CTn(size(n))
      !! Temperature/amount derivatives [L/(mol K)].
      real(pr), optional, intent(out) :: Cn2(size(n), size(n))
      !! Amount Hessian [L/mol²].

      C = dot_product(n, self%a + self%b*T)
      if (present(CT)) CT = dot_product(n, self%b)
      if (present(CT2)) CT2 = 0._pr
      if (present(Cn)) Cn = self%a + self%b*T
      if (present(CTn)) CTn = self%b
      if (present(Cn2)) Cn2 = 0._pr
   end subroutine linear_translation
end module yaeos__models_ar_volume_translations_linear
