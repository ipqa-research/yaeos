program test_constant_volume_translation
   use yaeos, only: pr, R, ArModel, CubicEoS, PengRobinson76, SoaveRedlichKwong, &
      ConstantVolumeTranslation, numeric_ar_derivatives, ar_consistency
   use testing_aux, only: assert, test_title
   implicit none

   type(CubicEoS) :: base !! Binary methane/ethane EOS.
   type(ConstantVolumeTranslation) :: model !! Constant volume translation of the base EOS.
   type(ConstantVolumeTranslation) :: zero_shift !! Adapter with zero component shifts.
   class(ArModel), allocatable :: polymorphic !! Public constructor result stored as an ArModel.
   real(pr), parameter :: n(2) = [0.8_pr, 1.5_pr] !! Nonnormalized mole amounts [mol].
   real(pr), parameter :: shifts(2) = [0.012_pr, 0.021_pr] !! Illustrative component shifts [L/mol].
   real(pr) :: coefficients(2) !! Mutable constructor input [L/mol].
   real(pr) :: reference(16) !! Packed base energy and derivatives.
   real(pr) :: values(16) !! Packed translated energy and derivatives.
   real(pr) :: pressure !! Pressure through the concrete adapter [bar].
   real(pr) :: polymorphic_pressure !! Pressure through the ArModel interface [bar].

   write(*, *) test_title("ConstantVolumeTranslation")
   base = PengRobinson76([190.564_pr, 305.32_pr], [45.99_pr, 48.72_pr], [0.01142_pr, 0.0995_pr])
   coefficients = shifts
   model = ConstantVolumeTranslation(base, c=coefficients)
   coefficients = 0._pr
   call assert(all(model%c == shifts), "Constructor owns its component coefficients")
   call assert(all(model%components%Tc == base%components%Tc), "Constructor preserves component order and metadata")
   call check_translation(model)
   call check_derivatives(model)
   call check_properties(base, model)

   zero_shift = ConstantVolumeTranslation(base, c=[0._pr, 0._pr])
   call evaluate(base, reference)
   call evaluate(zero_shift, values)
   call assert(all(abs(values-reference) <= 1.e-12_pr*max(1.e-5_pr, abs(reference))), &
      "Zero shifts recover every base Helmholtz derivative")
   call check_properties(base, zero_shift)

   model = ConstantVolumeTranslation(base, c=-shifts)
   call check_derivatives(model)
   call check_properties(base, model)

   ! The same implementation works with other EOS and via the abstract API.
   base = SoaveRedlichKwong([190.564_pr, 305.32_pr], [45.99_pr, 48.72_pr], [0.01142_pr, 0.0995_pr])
   model = ConstantVolumeTranslation(base, shifts)
   polymorphic = ConstantVolumeTranslation(base, shifts)
   call model%pressure(n, 0.8_pr, 330._pr, pressure)
   call polymorphic%pressure(n, 0.8_pr, 330._pr, polymorphic_pressure)
   call assert(abs(pressure-polymorphic_pressure) < 1.e-12_pr, "Constructor result works as an ordinary ArModel")
   call check_derivatives(model)
   call check_properties(base, model)
   call check_pure_component()
   write(*, *) ""

contains

   subroutine check_translation(eos)
      type(ConstantVolumeTranslation), intent(in) :: eos !! Model with the positive illustrative coefficients.
      real(pr) :: C !! Total translation [L].
      real(pr) :: CT !! Temperature derivative [L/K].
      real(pr) :: CT2 !! Second temperature derivative [L/K^2].
      real(pr) :: Cn(2) !! Amount gradient [L/mol].
      real(pr) :: CTn(2) !! Mixed temperature/amount derivatives [L/(K mol)].
      real(pr) :: Cn2(2, 2) !! Amount Hessian [L/mol^2].
      real(pr) :: expected !! Independent value for the specified amounts [L].
      integer :: i !! Temperature case.

      expected = 0.0411_pr
      do i = 1, 3
         call eos%translation(n, real(i, pr)*200._pr, C, CT, CT2, Cn, CTn, Cn2)
         call assert(abs(C-expected) < 1.e-15_pr, "C is the linear total-volume shift and is independent of T")
         call assert(all(Cn == shifts), "Amount derivatives equal the constant molar shifts")
         call assert(CT == 0._pr .and. CT2 == 0._pr, "Temperature derivatives are zero")
         call assert(all(CTn == 0._pr) .and. all(Cn2 == 0._pr), "Mixed derivatives and amount Hessian are zero")
      end do

      ! Every optional derivative is also usable without requesting the others.
      call eos%translation(n, 300._pr, C)
      call assert(abs(C-expected) < 1.e-15_pr, "Value-only translation call")
      call eos%translation(n, 300._pr, C, CT=CT)
      call assert(CT == 0._pr, "Individual CT call")
      call eos%translation(n, 300._pr, C, CT2=CT2)
      call assert(CT2 == 0._pr, "Individual CT2 call")
      call eos%translation(n, 300._pr, C, Cn=Cn)
      call assert(all(Cn == shifts), "Individual Cn call")
      call eos%translation(n, 300._pr, C, CTn=CTn)
      call assert(all(CTn == 0._pr), "Individual CTn call")
      call eos%translation(n, 300._pr, C, Cn2=Cn2)
      call assert(all(Cn2 == 0._pr), "Individual Cn2 call")

      call eos%translation(4._pr*n, 300._pr, C, Cn=Cn)
      call assert(abs(C-4._pr*expected) < 1.e-15_pr .and. all(Cn == shifts), &
         "C scales with amounts while its amount gradient is intensive")
      call eos%translation([0._pr, 2._pr], 300._pr, C)
      call assert(abs(C-0.042_pr) < 1.e-15_pr, "A zero component amount contributes no shift")
   end subroutine check_translation

   subroutine evaluate(eos, values)
      class(ArModel), intent(in) :: eos !! Model used through the shared residual interface.
      real(pr), intent(out) :: values(16) !! Six scalars, three gradients, and a flattened Hessian.
      real(pr) :: hessian(2, 2) !! Amount Hessian [bar L/mol^2].

      call eos%residual_helmholtz(n, 0.8_pr, 330._pr, &
         Ar=values(1), ArV=values(2), ArT=values(3), ArTV=values(4), &
         ArV2=values(5), ArT2=values(6), Arn=values(7:8), &
         ArVn=values(9:10), ArTn=values(11:12), Arn2=hessian)
      values(13:16) = reshape(hessian, [4])
   end subroutine evaluate

   subroutine check_derivatives(eos)
      type(ConstantVolumeTranslation), intent(in) :: eos !! Constant translated EOS.
      real(pr) :: analytic(16) !! Implemented energy and derivatives.
      real(pr) :: numeric(16) !! Independent energy finite differences.
      real(pr) :: hessian(2, 2) !! Numerical amount Hessian [bar L/mol^2].
      real(pr) :: eq31 !! Gibbs/fugacity identity residual.
      real(pr) :: eq33(2, 2) !! Fugacity amount-derivative symmetry residual.
      real(pr) :: eq34(2) !! Gibbs-Duhem residual.
      real(pr) :: eq36 !! Pressure/fugacity identity residual.
      real(pr) :: eq37 !! Enthalpy/fugacity identity residual.

      call evaluate(eos, analytic)
      call numeric_ar_derivatives(eos, n, 0.8_pr, 330._pr, d_n=2.e-4_pr, d_v=1.e-4_pr, d_t=0.05_pr, &
         Ar=numeric(1), ArV=numeric(2), ArT=numeric(3), ArTV=numeric(4), &
         ArV2=numeric(5), ArT2=numeric(6), Arn=numeric(7:8), &
         ArVn=numeric(9:10), ArTn=numeric(11:12), Arn2=hessian)
      numeric(13:16) = reshape(hessian, [4])
      call assert(all(abs(analytic-numeric) <= 3.e-5_pr*max(1.e-5_pr, abs(numeric))), &
         "Every inherited Helmholtz derivative agrees with energy finite differences")
      call ar_consistency(eos, n, 0.8_pr, 330._pr, eq31, eq33, eq34, eq36, eq37)
      call assert(maxval(abs([eq31, reshape(eq33, [4]), eq34, eq36, eq37])) < 1.e-12_pr, &
         "Constant translation satisfies the residual thermodynamic identities")
   end subroutine check_derivatives

   subroutine check_properties(base, eos)
      class(ArModel), intent(in) :: base !! Untranslated EOS.
      type(ConstantVolumeTranslation), intent(in) :: eos !! Its constant-shift counterpart.
      character(len=6), parameter :: roots(3) = ["liquid", "vapor ", "stable"]
      !! Requested homogeneous root labels.
      real(pr) :: C !! Total translation [L].
      real(pr) :: W !! Corresponding base volume [L].
      real(pr) :: V !! Translated physical volume [L].
      real(pr) :: pressure !! Pressure reconstructed from the translated energy [bar].
      real(pr) :: base_pressure !! Pressure at W from the base EOS [bar].
      real(pr) :: energy !! Translated residual Helmholtz energy [bar L].
      real(pr) :: base_energy !! Base residual Helmholtz energy [bar L].
      real(pr) :: lnphi(2) !! Translated log fugacity coefficients.
      real(pr) :: base_lnphi(2) !! Base log fugacity coefficients at W.
      real(pr) :: enthalpy !! Translated residual enthalpy [bar L].
      real(pr) :: base_enthalpy !! Base residual enthalpy at W [bar L].
      real(pr) :: Cp !! Translated residual heat capacity [bar L/K].
      real(pr) :: base_Cp !! Base residual heat capacity at W [bar L/K].
      integer :: i !! Requested root index.

      C = dot_product(n, eos%c)
      W = 0.8_pr + C
      call base%pressure(n, W, 330._pr, base_pressure)
      call eos%pressure(n, 0.8_pr, 330._pr, pressure)
      call assert(abs(pressure-base_pressure) < 1.e-11_pr, "Translated pressure is the base pressure at V+C")
      call base%residual_helmholtz(n, W, 330._pr, Ar=base_energy)
      call eos%residual_helmholtz(n, 0.8_pr, 330._pr, Ar=energy)
      call assert(abs(energy-base_energy-sum(n)*R*330._pr*log(0.8_pr/W)) < 1.e-12_pr, &
         "Residual Helmholtz contains the logarithmic reference correction")
      call base%lnphi_vt(n, W, 330._pr, lnPhi=base_lnphi)
      call eos%lnphi_vt(n, 0.8_pr, 330._pr, lnPhi=lnphi)
      call assert(maxval(abs(lnphi-base_lnphi+pressure*eos%c/(R*330._pr))) < 1.e-12_pr, &
         "Log fugacity correction is -P*c_i/(RT)")
      call base%enthalpy_residual_vt(n, W, 330._pr, Hr=base_enthalpy)
      call eos%enthalpy_residual_vt(n, 0.8_pr, 330._pr, Hr=enthalpy)
      call assert(abs(enthalpy-base_enthalpy+pressure*C) < 1.e-11_pr, &
         "Constant translation changes residual enthalpy by -P*C")
      call base%Cp_residual_vt(n, W, 330._pr, base_Cp)
      call eos%Cp_residual_vt(n, 0.8_pr, 330._pr, Cp)
      call assert(abs(Cp-base_Cp) < 1.e-12_pr, "Constant translation preserves Cp at corresponding P,T states")

      do i = 1, 3
         call base%volume(n, 1._pr, 200._pr, W, trim(roots(i)))
         call eos%volume(n, 1._pr, 200._pr, V, trim(roots(i)))
         call assert(V > 0._pr .and. abs(V-(W-C)) < 1.e-12_pr*W, &
            "Inherited volume solver translates liquid, vapor, and stable roots")
      end do
      call assert(abs(eos%get_v0(n, 1._pr, 200._pr) - (base%get_v0(n, 1._pr, 200._pr)-C)) < 1.e-14_pr, &
         "Inherited liquid initializer uses the extensive constant shift")
   end subroutine check_properties

   subroutine check_pure_component()
      type(CubicEoS) :: base !! Pure methane EOS.
      type(ConstantVolumeTranslation) :: model !! One-component constant translation.
      real(pr) :: C !! Shift for a nonunit pure-component amount [L].

      base = PengRobinson76([190.564_pr], [45.99_pr], [0.01142_pr])
      model = ConstantVolumeTranslation(base, c=[0.015_pr])
      call model%translation([2.5_pr], 300._pr, C)
      call assert(abs(C-0.0375_pr) < 1.e-15_pr, "Pure-component construction uses molar coefficients and total amounts")
   end subroutine check_pure_component

end program test_constant_volume_translation
