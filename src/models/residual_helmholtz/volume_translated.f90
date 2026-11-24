module yaeos__models_ar_volume_translated
   !! # Volume-translated residual Helmholtz models
   !! Abstract adapter for volume translations independent of volume and
   !! pressure.
   !!
   !! # Description
   !! With physical volume \(V\), base-model volume \(W=V+C(n,T)\), and
   !!
   !! \(N=\sum_i n_i\), the translated residual energy is
   !!
   !! \[
   !! A^r_{\mathrm{tr}}(n,V,T)
   !! = A^r_{\mathrm{base}}(n,W,T)+NRT\ln(V/W).
   !! \]
   !!
   !! Positive translations therefore reduce the physical volume at fixed
   !! pressure and temperature.
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
      !! Volume solves delegate to `base` and subtract the translation.
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
         !! Volume shift and independent amount/temperature derivatives.
         !!
         !! # Description
         !! `C` is always requested. Every other output is optional and must be
         !! computed when present, including when requested individually.
         !! Derivatives hold independent mole amounts fixed, not normalized
         !! mole fractions.
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
      self, n, V, T, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2 &
      )
      !! # translated_helmholtz
      !! Transform the base residual Helmholtz energy and requested derivatives.
      !!
      !! # Description
      !! Applies the chain rule at \(W=V+C(n,T)\), including the
      !! ideal-reference correction. Base derivatives are requested according
      !! to the outputs and their chain-rule dependencies.
      !!
      !! # Examples
      !! ```fortran
      !! call model%residual_helmholtz(n, V, T, ArT2=ArT2)
      !! ```
      class(VolumeTranslatedModel), intent(in) :: self !! Translated model.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
      real(pr), intent(in) :: V !! Physical total volume [L].
      real(pr), intent(in) :: T !! Temperature [K].
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

      ! Derivatives of the base residual and total Helmholtz energies at W.
      real(pr) :: ArV_base !! Base residual volume derivative [bar].
      real(pr) :: ArV2_base !! Base residual second volume derivative [bar/L].
      real(pr) :: ArTV_base
      !! Base residual temperature/volume derivative [bar/K].
      real(pr) :: ArVn_base(size(n))
      !! Base residual volume/amount gradient [bar/mol].

      real(pr) :: Av !! Base total Helmholtz volume derivative [bar].
      real(pr) :: Avv !! Base total Helmholtz second volume derivative [bar/L].
      real(pr) :: Avt
      !! Base total Helmholtz temperature/volume derivative [bar/K].
      real(pr) :: Avn(size(n))
      !! Base total Helmholtz volume/amount gradient [bar/mol].

      ! Translations
      real(pr) :: C !! Extensive volume translation [L].
      real(pr) :: W !! Base-model total volume V+C [L].
      real(pr) :: CT !! Translation temperature derivative [L/K].
      real(pr) :: CT2 !! Translation second temperature derivative [L/K^2].
      real(pr) :: Cn(size(n)) !! Translation amount gradient [L/mol].
      real(pr) :: CTn(size(n))
      !! Translation temperature/amount gradient [L/(K mol)].
      real(pr) :: Cn2(size(n), size(n)) !! Translation amount Hessian [L/mol^2].

      ! Auxiliars
      real(pr) :: NRT !! Total amount times RT [bar L].
      real(pr) :: NR !! Total amount times R [bar L/K].
      real(pr) :: log_ratio
      !! Logarithm of physical/base total volume ratio, log(V/W).
      real(pr) :: inverse_volume_difference !! Stable difference 1/V-1/W [1/L].
      integer :: i !! Row index of the amount Hessian.
      integer :: j !! Column index of the amount Hessian.

      ! Requested derivative outputs.
      logical :: energy_requested !! Residual energy Ar was requested.
      logical :: dV !! Volume derivative ArV was requested.
      logical :: dV2 !! Second volume derivative ArV2 was requested.
      logical :: dT !! Temperature derivative ArT was requested.
      logical :: dT2 !! Second temperature derivative ArT2 was requested.
      logical :: dn !! Amount gradient Arn was requested.
      logical :: dn2 !! Amount Hessian Arn2 was requested.
      logical :: dVn !! Volume/amount gradient ArVn was requested.
      logical :: dTn !! Temperature/amount gradient ArTn was requested.
      logical :: dTV !! Temperature/volume derivative ArTV was requested.

      ! ========================================================================
      ! Logicals for clarity
      ! ------------------------------------------------------------------------
      energy_requested = present(Ar)
      dV = present(ArV)
      dV2 = present(ArV2)
      dT = present(ArT)
      dT2 = present(ArT2)
      dn = present(Arn)
      dn2 = present(Arn2)
      dVn = present(ArVn)
      dTn = present(ArTn)
      dTV = present(ArTV)

      ! ========================================================================
      ! Calculate translation an V_eos (W)
      ! ask translation derivatives as efficient as possible
      ! ------------------------------------------------------------------------
      if (.not. (dT .or. dTV .or. dT2 .or. dTn)) then
         if (dn2) then
            call self%translation(n, T, C, Cn=Cn, Cn2=Cn2)
         else if (dn .or. dVn) then
            call self%translation(n, T, C, Cn=Cn)
         else
            call self%translation(n, T, C)
         end if
      else if (dn2) then
         call self%translation(n, T, C, CT=CT, CT2=CT2, Cn=Cn, CTn=CTn, Cn2=Cn2)
      else if (dTn) then
         call self%translation(n, T, C, CT=CT, CT2=CT2, Cn=Cn, CTn=CTn)
      else if (dn .or. dVn) then
         call self%translation(n, T, C, CT=CT, CT2=CT2, Cn=Cn)
      else
         call self%translation(n, T, C, CT=CT, CT2=CT2)
      end if

      ! Volume of the cubic with the inverse translation
      W = V + C

      ! ========================================================================
      ! Call base Ar as efficient as possible (normal call cases)
      ! ------------------------------------------------------------------------
      if (dTn .or. dn2 .or. (dVn .and. (dT2 .or. dTV))) then
         call self%base%residual_helmholtz(&
            n, W, t, Ar=Ar, ArV=ArV_base, ArT=ArT, ArTV=ArTV_base, &
            ArV2=ArV2_base, ArT2=ArT2, Arn=Arn, ArVn=ArVn_base, &
            ArTn=ArTn, Arn2=Arn2 &
            )
      else if (dVn) then
         call self%base%residual_helmholtz(&
            n, W, t, Ar=Ar, ArV=ArV_base, ArT=ArT, ArTV=ArTV, &
            ArV2=ArV2_base, ArT2=ArT2, Arn=Arn, ArVn=ArVn_base &
            )
      else if (dT2 .or. dTV) then
         call self%base%residual_helmholtz(&
            n, W, t, Ar=Ar, ArV=ArV_base, ArT=ArT, ArTV=ArTV_base, &
            ArV2=ArV2_base, ArT2=ArT2, Arn=Arn &
            )
      else if (dn .or. dT) then
         if (dV2) then
            call self%base%residual_helmholtz(&
               n, W, t, Ar=Ar, ArV=ArV_base, ArT=ArT, ArV2=ArV2_base, Arn=Arn &
               )
         else
            call self%base%residual_helmholtz(&
               n, W, t, Ar=Ar, ArV=ArV_base, ArT=ArT, Arn=Arn &
               )
         end if
      else if (dV2) then
         if (dV) then
            call self%base%residual_helmholtz(&
               n, W, t, Ar=Ar, ArV=ArV_base, ArV2=ArV2_base &
               )
         else
            call self%base%residual_helmholtz(n, W, t, Ar=Ar, ArV2=ArV2_base)
         end if
      else if (dV) then
         call self%base%residual_helmholtz(n, W, t, Ar=Ar, ArV=ArV_base)
      else
         call self%base%residual_helmholtz(n, W, t, Ar=Ar)
      end if

      ! ========================================================================
      ! Volume translated Helmholtz evaluations
      ! ------------------------------------------------------------------------
      NR = sum(n)*R
      NRT = NR*t
      log_ratio = log(V/W)
      inverse_volume_difference = (C/V)/W

      ! Volume derivatives
      if (dV) ArV = ArV_base + NRT*inverse_volume_difference

      if (dV2) then
         ArV2 = ArV2_base - NRT*inverse_volume_difference*(1._pr/V + 1._pr/W)
      end if

      ! These total Helmholtz derivatives collect the chain-rule terms due
      ! to C(n,T). The ideal-gas terms here are evaluated at the base volume W.
      if (dT .or. dT2 .or. dn .or. dTn .or. dn2) Av = ArV_base - NRT/W
      if (dTV .or. dT2 .or. dVn .or. dTn .or. dn2) Avv = ArV2_base + NRT/W**2
      if (dT2 .or. dTn) Avt = ArTV_base - NR/W
      if (dTn .or. dn2) Avn = ArVn_base - R*T/W

      ! Rest of derivatives
      if (energy_requested) Ar = Ar + NRT*log_ratio

      if (dT) ArT = ArT + Av*CT + NR*log_ratio

      if (dTV) ArTV = ArTV_base + NR*inverse_volume_difference + Avv*CT

      if (dT2) ArT2 = ArT2 + 2._pr*Avt*CT + Avv*CT**2 + Av*CT2

      if (dn) Arn = Arn + Av*Cn + R*T*log_ratio

      if (dVn) ArVn = ArVn_base + R*T*inverse_volume_difference + Avv*Cn

      if (dTn) ArTn = ArTn + Avt*Cn + Avn*CT + Avv*CT*Cn + Av*CTn + R*log_ratio

      if (dn2) then
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
      !! Returns `base%get_v0(n,p,t)-C(n,t)`.
      !!
      !! # Examples
      !! ```fortran
      !! V0 = model%get_v0(n, P, T)
      !! ```
      class(VolumeTranslatedModel), intent(in) :: self !! Translated model.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
      real(pr), intent(in) :: p !! Pressure [bar].
      real(pr), intent(in) :: t !! Temperature [K].
      real(pr) :: v0 !! Translated initial volume [L].

      real(pr) :: C
      !! Total volume translation at the specified amounts and temperature [L].

      call self%translation(n, t, C)
      v0 = self%base%get_v0(n, p, t) - C
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
      !! Physical total volume [L].
      character(len=*), intent(in) :: root_type
      !! Desired root-type to solve. Options are:
      !! `["liquid", "vapor", "stable"]

      real(pr) :: C
      !! Total volume translation at the specified amounts and temperature [L].

      call eos%base%volume(n, P, T, V, root_type)
      call eos%translation(n, T, C)

      V = V - C
   end subroutine translated_volume

end module yaeos__models_ar_volume_translated
