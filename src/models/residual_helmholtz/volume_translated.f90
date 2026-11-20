module yaeos__models_ar_volume_translated
   !! # Volume-translated residual Helmholtz models
   !! Abstract adapter for volume translations independent of volume and pressure.
   !!
   !! # Description
   !! With physical volume \(V\), base-model volume \(W=V+C(n,T)\), and
   !! \(N=\sum_i n_i\), the translated residual energy is
   !! \[
   !! A^r_{\mathrm{tr}}(n,V,T)
   !! = A^r_{\mathrm{base}}(n,W,T)+NRT\ln(V/W).
   !! \]
   !! Positive translations therefore reduce the physical volume at fixed
   !! pressure and temperature. The logarithmic term accounts for the change
   !! in the ideal-gas reference volume; no ideal caloric model is required.
   !! Extensions supply only the extensive translation and its derivatives.
   !!
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use yaeos__constants, only: pr, R
   use yaeos__models_ar, only: ArModel
   implicit none
   private

   public :: VolumeTranslatedModel

   type, abstract, extends(ArModel) :: VolumeTranslatedModel
      !! # VolumeTranslatedModel
      !! Residual Helmholtz model with a configurable volume-translation law.
      !!
      !! # Description
      !! Implement `translation` in a concrete extension and call `init` with
      !! an unmodified `ArModel`. The translation \(C(n,T)\) must be homogeneous
      !! of degree one in mole amounts, twice differentiable, and independent
      !! of volume and pressure. Component ordering must match the base model.
      !! Translations depending on a cubic alpha function can inspect `base`
      !! in their implementation; such extensions must enforce their required
      !! base-model type in their own constructors.
      !!
      !! Energy evaluation requires positive finite physical and base volumes.
      !! Volume solves delegate to `base` and subtract the translation;
      !! failed or nonpositive/nonfinite translated roots return `-1`.
      !!
      !! # Examples
      !! ```fortran
      !! ! model is an instance of a concrete extension.
      !! call model%init(PengRobinson76(Tc, Pc, omega))
      !! call model%pressure(n, V, T, P)
      !! ```
      class(ArModel), allocatable :: base
      !! Owned copy of the untranslated model, available to translation laws.
   contains
      procedure :: init
      procedure(translation_interface), deferred :: translation
      procedure :: residual_helmholtz => translated_helmholtz
      procedure :: get_v0 => translated_v0
      procedure :: volume => translated_volume
   end type VolumeTranslatedModel

   abstract interface
      subroutine translation_interface(self, n, T, C, CT, CT2, Cn, CTn, Cn2)
         !! # translation_interface
         !! Extensive volume shift and independent amount/temperature derivatives.
         !!
         !! # Description
         !! `C` is always requested. Every other output is optional and must be
         !! computed when present, including when requested individually.
         !! Derivatives hold independent mole amounts fixed, not normalized
         !! mole fractions. A constant-per-component law has
         !! \(C=\sum_i n_i c_i\), `Cn=c`, and zero remaining derivatives.
         !!
         !! # Examples
         !! ```fortran
         !! call model%translation(n, T, C, CT=CT)
         !! ```
         import VolumeTranslatedModel, pr
         class(VolumeTranslatedModel), intent(in) :: self
         !! Translated model.
         real(pr), intent(in) :: n(:)
         !! Component mole amounts [mol].
         real(pr), intent(in) :: T
         !! Temperature [K].
         real(pr), intent(out) :: C
         !! Total volume translation [L].
         real(pr), optional, intent(out) :: CT
         !! Temperature derivative [L/K].
         real(pr), optional, intent(out) :: CT2
         !! Second temperature derivative [L/K^2].
         real(pr), optional, intent(out) :: Cn(size(n))
         !! Amount gradient [L/mol].
         real(pr), optional, intent(out) :: CTn(size(n))
         !! Temperature/amount derivatives [L/(K mol)].
         real(pr), optional, intent(out) :: Cn2(size(n), size(n))
         !! Amount Hessian [L/mol^2].
      end subroutine translation_interface
   end interface

contains

   subroutine init(self, base)
      !! # init
      !! Store an owned copy of the base EOS and its component metadata.
      !!
      !! # Description
      !! Translation parameters belong to the concrete extension and are not
      !! modified here. Reinitialization replaces the base model and metadata.
      !!
      !! # Examples
      !! ```fortran
      !! call model%init(base)
      !! ```
      class(VolumeTranslatedModel), intent(inout) :: self
      !! Model being initialized.
      class(ArModel), intent(in) :: base
      !! Untranslated model in the same component order.

      self%base = base
      self%components = base%components
      if (allocated(base%name)) then
         self%name = base%name // " (volume translated)"
      else
         self%name = "VolumeTranslatedModel"
      end if
   end subroutine init

   subroutine translated_helmholtz(&
      self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2 &
      )
      !! # translated_helmholtz
      !! Transform the base residual Helmholtz energy and requested derivatives.
      !!
      !! # Description
      !! Applies the chain rule at \(W=V+C(n,T)\), including the ideal-reference
      !! correction. Scalar volume derivatives are evaluated together; the
      !! amount-volume derivatives and amount Hessians are requested only
      !! when needed. Call `init` before evaluation. Invalid physical/base
      !! volumes cause `error stop`; a call without outputs is a no-op.
      !!
      !! # Examples
      !! ```fortran
      !! call model%residual_helmholtz(n, V, T, ArT2=ArT2)
      !! ```
      class(VolumeTranslatedModel), intent(in) :: self !! Translated model.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
      real(pr), intent(in) :: v !! Physical total volume [L].
      real(pr), intent(in) :: t !! Temperature [K].
      real(pr), optional, intent(out) :: Ar
      !! Residual Helmholtz energy [bar L].
      real(pr), optional, intent(out) :: ArV
      !! Volume derivative [bar].
      real(pr), optional, intent(out) :: ArT
      !! Temperature derivative [bar L/K].
      real(pr), optional, intent(out) :: ArTV
      !! Temperature/volume derivative [bar/K].
      real(pr), optional, intent(out) :: ArV2
      !! Second volume derivative [bar/L].
      real(pr), optional, intent(out) :: ArT2
      !! Second temperature derivative [bar L/K^2].
      real(pr), optional, intent(out) :: Arn(size(n))
      !! Amount gradient [bar L/mol].
      real(pr), optional, intent(out) :: ArVn(size(n))
      !! Volume/amount derivatives [bar/mol].
      real(pr), optional, intent(out) :: ArTn(size(n))
      !! Temperature/amount derivatives [bar L/(K mol)].
      real(pr), optional, intent(out) :: Arn2(size(n), size(n))
      !! Amount Hessian [bar L/mol^2].

      ! Base residual Helholtz Evaluations
      real(pr) :: ArV_base
      real(pr) :: ArV2_base
      real(pr) :: ArTV_base
      real(pr) :: ArVn_base(size(n))
      real(pr) :: Av ! Base total Helmholtz volume derivative [bar].
      real(pr) :: Avv ! Base total Helmholtz second volume derivative [bar/L].
      real(pr) :: Avt ! Base total Helmholtz volume/temperature derivative [bar/K].
      real(pr) :: Avn(size(n)) ! Base total Helmholtz volume/amount derivatives [bar/mol].

      ! Translations
      real(pr) :: C
      real(pr) :: W ! Base-model volume V+C
      real(pr) :: CT
      real(pr) :: CT2
      real(pr) :: Cn(size(n))
      real(pr) :: CTn(size(n))
      real(pr) :: Cn2(size(n), size(n))

      ! Auxiliars
      real(pr) :: NRT
      real(pr) :: NR
      real(pr) :: log_ratio ! Logarithm of physical/base volume ratio.
      real(pr) :: inverse_volume_difference ! Stable evaluation of 1/V-1/W.

      integer :: i, j

      ! Gather the shift and its derivatives. Only composition-dependent
      if (present(Arn2)) then
         call self%translation(n, t, C, CT=CT, CT2=CT2, Cn=Cn, CTn=CTn, Cn2=Cn2)
      else if (present(ArTn)) then
         call self%translation(n, t, C, CT=CT, CT2=CT2, Cn=Cn, CTn=CTn)
      else if (present(Arn) .or. present(ArVn)) then
         call self%translation(n, t, C, CT=CT, CT2=CT2, Cn=Cn)
      else
         call self%translation(n, t, C, CT=CT, CT2=CT2)
      end if

      ! Volume of the cubic with the inverse translation
      W = v + C

      ! Evaluate the base EOS once at W. TODO: This could be more efficient
      if (present(ArVn) .or. present(ArTn) .or. present(Arn2)) then
         call self%base%residual_helmholtz(&
            n, W, t, Ar=Ar, ArV=ArV_base, ArT=ArT, ArTV=ArTV_base, &
            ArV2=ArV2_base, ArT2=ArT2, Arn=Arn, ArVn=ArVn_base, &
            ArTn=ArTn, Arn2=Arn2 &
            )
      else
         call self%base%residual_helmholtz(&
            n, W, t, Ar=Ar, ArV=ArV_base, ArT=ArT, ArTV=ArTV_base, &
            ArV2=ArV2_base, ArT2=ArT2, Arn=Arn &
            )
      end if

      NR = sum(n)*R
      NRT = NR*t
      log_ratio = log(v/W)
      inverse_volume_difference = (C/v)/W
      ! Retain the base residual derivatives directly to avoid subtracting and
      ! re-adding the ideal term, especially for zero shifts at low density.
      if (present(ArV)) ArV = ArV_base + NRT*inverse_volume_difference
      if (present(ArV2)) ArV2 = ArV2_base - NRT*inverse_volume_difference*(1._pr/v + 1._pr/W)

      ! These total Helmholtz derivatives collect the chain-rule terms due
      ! to C(n,T). The ideal-gas terms here are evaluated at the base volume W.
      Av = ArV_base - NRT/W
      Avv = ArV2_base + NRT/W**2
      Avt = ArTV_base - NR/W
      if (present(ArVn) .or. present(ArTn) .or. present(Arn2)) Avn = ArVn_base - R*t/W

      if (present(Ar)) Ar = Ar + NRT*log_ratio
      if (present(ArT)) ArT = ArT + Av*CT + NR*log_ratio
      if (present(ArTV)) ArTV = ArTV_base + NR*inverse_volume_difference + Avv*CT
      if (present(ArT2)) ArT2 = ArT2 + 2._pr*Avt*CT + Avv*CT**2 + Av*CT2
      if (present(Arn)) Arn = Arn + Av*Cn + R*t*log_ratio
      if (present(ArVn)) ArVn = ArVn_base + R*t*inverse_volume_difference + Avv*Cn
      if (present(ArTn)) ArTn = ArTn + Avt*Cn + Avn*CT + Avv*CT*Cn + Av*CTn + R*log_ratio
      if (present(Arn2)) then
         do j = 1, size(n)
            do i = 1, size(n)
               Arn2(i, j) = Arn2(i, j) + Avn(i)*Cn(j) + Avn(j)*Cn(i) &
                  + Avv*Cn(i)*Cn(j) + Av*Cn2(i, j)
            end do
         end do
      end if
   end subroutine translated_helmholtz

   function translated_v0(self, n, p, t) result(v0)
      !! # translated_v0
      !! Translate the base model's liquid-volume initializer.
      !!
      !! # Description
      !! Returns `base%get_v0(n,p,t)-C(n,t)`. Invalid base estimates or
      !! nonpositive/nonfinite translated estimates return `-1`. This is an
      !! initializer, not an exact volume root; a nonpositive translated
      !! estimate does not imply that all translated roots are inadmissible.
      !!
      !! # Examples
      !! ```fortran
      !! V0 = model%get_v0(n, P, T)
      !! ```
      class(VolumeTranslatedModel), intent(in) :: self !! Translated model.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
      real(pr), intent(in) :: p !! Pressure [bar].
      real(pr), intent(in) :: t !! Temperature [K].
      real(pr) :: v0 !! Translated initial volume [L], or -1 on failure.
      real(pr) :: C !! Total volume translation [L].

      v0 = self%base%get_v0(n, p, t)

      ! Do not translate a failed base estimate into an apparently valid one.
      if (.not. ieee_is_finite(v0)) then
         v0 = -1._pr
         return
      end if

      if (v0 <= 0._pr) then
         v0 = -1._pr
         return
      end if

      call self%translation(n, t, C)
      v0 = v0 - C

      ! The shifted estimate must still be a positive physical volume.
      if (.not. ieee_is_finite(v0)) then
         v0 = -1._pr
      else if (v0 <= 0._pr) then
         v0 = -1._pr
      end if
   end function translated_v0

   subroutine translated_volume(eos, n, P, T, V, root_type)
      !! # translated_volume
      !! Preserve the base model's volume solver and translate its result.
      !!
      !! # Description
      !! At fixed `n,T`, pressure is the base pressure at `V+C(n,T)`.
      !! Root labels are forwarded unchanged. Gibbs energies of corresponding
      !! roots differ by the root-independent term `-P*C`, so stable-root
      !! selection is preserved when the corresponding roots are admissible.
      !! Invalid base roots or nonpositive/nonfinite physical roots return `-1`;
      !! no alternative-root search is attempted for an inadmissible result.
      !!
      !! # Examples
      !! ```fortran
      !! call model%volume(n, P, T, V, root_type="liquid")
      !! ```
      class(VolumeTranslatedModel), intent(in) :: eos
      !! Translated model.
      real(pr), intent(in) :: n(:)
      !! Component mole amounts [mol].
      real(pr), intent(in) :: P
      !! Pressure [bar].
      real(pr), intent(in) :: T
      !! Temperature [K].
      real(pr), intent(out) :: V
      !! Physical total volume [L], or -1 on failure.
      character(len=*), intent(in) :: root_type
      !! Base solver root selection.

      real(pr) :: C

      call eos%base%volume(n, P, T, V, root_type)

      ! Preserve a failed solve before applying the volume correction.
      if (.not. ieee_is_finite(V)) then
         V = -1._pr
         return
      end if
      if (V <= 0._pr) then
         V = -1._pr
         return
      end if

      call eos%translation(n, T, C)
      V = V - C

      ! A successful base root can still give an inadmissible physical volume.
      if (.not. ieee_is_finite(V)) then
         V = -1._pr
      else if (V <= 0._pr) then
         V = -1._pr
      end if
   end subroutine translated_volume

end module yaeos__models_ar_volume_translated
