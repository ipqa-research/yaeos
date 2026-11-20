module volume_translation_test_models
   !! Test-only translation laws and a non-cubic model with prescribed roots.
   use yaeos, only: pr, ArModel, VolumeTranslatedModel
   implicit none
   private
    public :: TestTranslation, ZeroResidualModel, TrackingResidualModel, base_requested
    logical :: base_requested(10) = .false.

   type, extends(VolumeTranslatedModel) :: TestTranslation
      real(pr) :: c0(2) = [0.012_pr, 0.021_pr] !! Constant molar shifts [L/mol].
      real(pr) :: c1(2) = 0._pr !! Linear temperature coefficients [L/(mol K)].
      real(pr) :: c2(2) = 0._pr !! Quadratic temperature coefficients [L/(mol K^2)].
      real(pr) :: k(3) = 0._pr !! Coefficients of a temperature-dependent composition term.
   contains
      procedure :: translation => test_translation
   end type TestTranslation

    type, extends(ArModel) :: ZeroResidualModel
      real(pr) :: initial_volume = 0.2_pr !! Prescribed volume initializer [L].
      real(pr) :: root_volume = 1._pr !! Prescribed volume root or failure sentinel [L].
   contains
      procedure :: residual_helmholtz => zero_helmholtz
      procedure :: get_v0 => zero_v0
      procedure :: volume => prescribed_volume
    end type ZeroResidualModel

    type, extends(ZeroResidualModel) :: TrackingResidualModel
    contains
       procedure :: residual_helmholtz => tracked_helmholtz
    end type TrackingResidualModel

contains

   subroutine test_translation(self, n, T, C, CT, CT2, Cn, CTn, Cn2)
      !! Homogeneous nonlinear law: C = n.c(T) + k(T)*n1*n2/sum(n).
      class(TestTranslation), intent(in) :: self !! Translation parameters and base EOS.
      real(pr), intent(in) :: n(:) !! Component amounts [mol].
      real(pr), intent(in) :: T !! Temperature [K].
      real(pr), intent(out) :: C !! Extensive volume shift [L].
      real(pr), optional, intent(out) :: CT !! Temperature derivative [L/K].
      real(pr), optional, intent(out) :: CT2 !! Second temperature derivative [L/K^2].
      real(pr), optional, intent(out) :: Cn(size(n)) !! Amount gradient [L/mol].
      real(pr), optional, intent(out) :: CTn(size(n)) !! Mixed derivatives [L/(K mol)].
      real(pr), optional, intent(out) :: Cn2(size(n), size(n)) !! Amount Hessian [L/mol^2].
      real(pr) :: dt !! Temperature offset from 300 K [K].
      real(pr) :: total !! Total mole amount [mol].
      real(pr) :: f !! Nonlinear extensive composition factor [mol].
      real(pr) :: fn(2) !! Gradient of the composition factor.
      real(pr) :: fnn(2, 2) !! Hessian of the composition factor [1/mol].
      real(pr) :: ci(2) !! Temperature-dependent pure-component shifts [L/mol].
      real(pr) :: ciT(2) !! Pure-component temperature derivatives [L/(mol K)].
      real(pr) :: kt !! Temperature-dependent nonlinear coefficient [L/mol].
      real(pr) :: dkt !! Its first temperature derivative [L/(mol K)].

      dt = T - 300._pr
      total = sum(n)
      f = n(1)*n(2)/total
      fn = [n(2)**2, n(1)**2]/total**2
      fnn(1, 1) = -2._pr*n(2)**2/total**3
      fnn(2, 2) = -2._pr*n(1)**2/total**3
      fnn(1, 2) = 2._pr*n(1)*n(2)/total**3
      fnn(2, 1) = fnn(1, 2)
      ci = self%c0 + self%c1*dt + self%c2*dt**2
      ciT = self%c1 + 2._pr*self%c2*dt
      kt = self%k(1) + self%k(2)*dt + self%k(3)*dt**2
      dkt = self%k(2) + 2._pr*self%k(3)*dt

      C = dot_product(n, ci) + kt*f
      if (present(CT)) CT = dot_product(n, ciT) + dkt*f
      if (present(CT2)) CT2 = 2._pr*(dot_product(n, self%c2) + self%k(3)*f)
      if (present(Cn)) Cn = ci + kt*fn
      if (present(CTn)) CTn = ciT + dkt*fn
      if (present(Cn2)) Cn2 = kt*fnn
   end subroutine test_translation

    subroutine zero_helmholtz(self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2)
      class(ZeroResidualModel), intent(in) :: self !! Unused; the residual energy is identically zero.
      real(pr), intent(in) :: n(:) !! Amounts defining the output shapes [mol].
      real(pr), intent(in) :: v !! Unused physical volume [L].
      real(pr), intent(in) :: t !! Unused temperature [K].
      real(pr), optional, intent(out) :: Ar !! Zero residual energy [bar L].
      real(pr), optional, intent(out) :: ArV !! Zero volume derivative [bar].
      real(pr), optional, intent(out) :: ArT !! Zero temperature derivative [bar L/K].
      real(pr), optional, intent(out) :: ArTV !! Zero mixed derivative [bar/K].
      real(pr), optional, intent(out) :: ArV2 !! Zero second volume derivative [bar/L].
      real(pr), optional, intent(out) :: ArT2 !! Zero second temperature derivative [bar L/K^2].
      real(pr), optional, intent(out) :: Arn(size(n)) !! Zero amount gradient [bar L/mol].
      real(pr), optional, intent(out) :: ArVn(size(n)) !! Zero volume/amount derivatives [bar/mol].
      real(pr), optional, intent(out) :: ArTn(size(n)) !! Zero temperature/amount derivatives [bar L/(K mol)].
      real(pr), optional, intent(out) :: Arn2(size(n), size(n)) !! Zero amount Hessian [bar L/mol^2].

      if (present(Ar)) Ar = 0._pr
      if (present(ArV)) ArV = 0._pr
      if (present(ArT)) ArT = 0._pr
      if (present(ArTV)) ArTV = 0._pr
      if (present(ArV2)) ArV2 = 0._pr
      if (present(ArT2)) ArT2 = 0._pr
      if (present(Arn)) Arn = 0._pr
      if (present(ArVn)) ArVn = 0._pr
      if (present(ArTn)) ArTn = 0._pr
      if (present(Arn2)) Arn2 = 0._pr
    end subroutine zero_helmholtz

    subroutine tracked_helmholtz(self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2)
       class(TrackingResidualModel), intent(in) :: self
       real(pr), intent(in) :: n(:)
       real(pr), intent(in) :: v
       real(pr), intent(in) :: t
       real(pr), optional, intent(out) :: Ar
       real(pr), optional, intent(out) :: ArV
       real(pr), optional, intent(out) :: ArT
       real(pr), optional, intent(out) :: ArTV
       real(pr), optional, intent(out) :: ArV2
       real(pr), optional, intent(out) :: ArT2
       real(pr), optional, intent(out) :: Arn(size(n))
       real(pr), optional, intent(out) :: ArVn(size(n))
       real(pr), optional, intent(out) :: ArTn(size(n))
       real(pr), optional, intent(out) :: Arn2(size(n), size(n))

       base_requested = [present(Ar), present(ArV), present(ArT), present(ArTV), &
          present(ArV2), present(ArT2), present(Arn), present(ArVn), present(ArTn), present(Arn2)]
       call zero_helmholtz(self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2)
    end subroutine tracked_helmholtz

   function zero_v0(self, n, p, t) result(v0)
      class(ZeroResidualModel), intent(in) :: self !! Model carrying the prescribed estimate.
      real(pr), intent(in) :: n(:) !! Unused amounts [mol].
      real(pr), intent(in) :: p !! Unused pressure [bar].
      real(pr), intent(in) :: t !! Unused temperature [K].
      real(pr) :: v0 !! Prescribed initial volume [L].

      v0 = self%initial_volume
   end function zero_v0

   subroutine prescribed_volume(eos, n, P, T, V, root_type)
      class(ZeroResidualModel), intent(in) :: eos !! Model carrying the prescribed root.
      real(pr), intent(in) :: n(:) !! Unused amounts [mol].
      real(pr), intent(in) :: P !! Unused pressure [bar].
      real(pr), intent(in) :: T !! Unused temperature [K].
      real(pr), intent(out) :: V !! Prescribed volume or failure sentinel [L].
      character(len=*), intent(in) :: root_type !! Unused requested root label.

      V = eos%root_volume
   end subroutine prescribed_volume
end module volume_translation_test_models

program test_volume_translated
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use yaeos, only: pr, R, ArModel, CubicEoS, PengRobinson76, SoaveRedlichKwong, &
      numeric_ar_derivatives, ar_consistency, size, EquilibriumState, flash
    use volume_translation_test_models, only: TestTranslation, ZeroResidualModel, TrackingResidualModel, base_requested
   use testing_aux, only: assert, test_title
   implicit none

   type(CubicEoS) :: base !! Binary methane/ethane cubic EOS.
   type(TestTranslation) :: model !! Concrete test extension of the abstract adapter.
   real(pr), parameter :: n(2) = [0.8_pr, 1.5_pr] !! Nonnormalized component amounts [mol].
   real(pr), parameter :: V = 0.8_pr !! Physical test volume [L].
   real(pr), parameter :: T = 330._pr !! Test temperature [K].
   real(pr) :: C !! Total shift at the test state [L].
   real(pr) :: translated_values(16) !! Packed energy and all first/second derivatives.
   real(pr) :: base_values(16) !! Packed untranslated reference derivatives.
   character(len=40) :: mode !! Optional properties-only mode or expected-failure scenario.

   call get_command_argument(1, mode)
   if (len_trim(mode) > 0 .and. mode /= "properties-only") call reject_invalid_state(trim(mode))

   write(*, *) test_title("VolumeTranslatedModel")
   base = PengRobinson76([190.564_pr, 305.32_pr], [45.99_pr, 48.72_pr], [0.01142_pr, 0.0995_pr])
   call model%init(base)
   call assert(size(model) == 2, "Adapter preserves component count")
   call assert(all(model%components%Tc == base%components%Tc), "Adapter preserves critical temperatures")
   call assert(all(model%components%Pc == base%components%Pc), "Adapter preserves critical pressures")
   call assert(all(model%components%w == base%components%w), "Adapter preserves acentric factors")
   call assert(index(model%name, base%name) > 0, "Adapter name identifies the base model")
   base%ac = 0._pr
   select type (owned => model%base)
    type is (CubicEoS)
      call assert(all(owned%ac > 0._pr), "Initialization owns an independent copy of the base model")
   end select

   ! Zero translation must recover every residual derivative.
   model%c0 = 0._pr
   call evaluate(model, n, V, T, translated_values)
   call evaluate(model%base, n, V, T, base_values)
   call assert(close(translated_values, base_values, 1.e-12_pr), "Zero shift recovers the base energy and derivatives")
   call evaluate(model, n, 1.e6_pr, T, translated_values)
   call evaluate(model%base, n, 1.e6_pr, T, base_values)
   call assert(all(translated_values == base_values), "Zero shift does not lose dilute-gas residual derivatives")
   call check_optional_outputs(model)
   call check_derivatives(model)

   ! Constant, negative, temperature-dependent, and nonlinear composition laws.
   model%c0 = [0.012_pr, 0.021_pr]
   call check_derivatives(model)
   call check_optional_outputs(model)
   call check_roots(model)
   ! Allow property diagnostics independently of the existing equilibrium solvers.
   if (mode /= "properties-only") call check_flash(model)
   model%c0 = -model%c0
   call check_derivatives(model)
   call check_roots(model)
   model%c0 = -model%c0
   model%c1 = [2.e-5_pr, -3.e-5_pr]
   model%c2 = [1.e-7_pr, 2.e-7_pr]
   model%k = [0.009_pr, 3.e-5_pr, 2.e-7_pr]
   call check_derivatives(model)
   call check_optional_outputs(model)
   call check_roots(model)
   call check_consistency(model)
   call check_scaling(model)
   call check_pressure_and_calorics(model)

   call model%translation(n, T, C)
   call model%base%residual_helmholtz(n, V+C, T, Ar=base_values(1))
   call model%residual_helmholtz(n, V, T, Ar=translated_values(1))
   call assert(abs(translated_values(1) - base_values(1) - sum(n)*R*T*log(V/(V+C))) < 1.e-12_pr, &
      "Translated Helmholtz includes the ideal-reference correction")

   call model%init(SoaveRedlichKwong([190.564_pr, 305.32_pr], &
      [45.99_pr, 48.72_pr], [0.01142_pr, 0.0995_pr]))
   call assert(index(model%name, "SRK") > 0, "Reinitialization replaces model metadata")
   call check_derivatives(model)
   call check_roots(model)
    call check_non_cubic_and_failures()
    call check_base_requests()
   write(*, *) ""

contains

   logical function close(actual, expected, tolerance)
      real(pr), intent(in) :: actual(:) !! Values to compare.
      real(pr), intent(in) :: expected(:) !! Reference values in the same order.
      real(pr), intent(in) :: tolerance !! Relative tolerance with a small absolute scale floor.

      close = all(abs(actual-expected) <= tolerance*max(1.e-5_pr, abs(expected)))
   end function close

   subroutine evaluate(eos, amounts, volume, temperature, values)
      class(ArModel), intent(in) :: eos !! Model accepting the ordinary residual-energy interface.
      real(pr), intent(in) :: amounts(2) !! Component amounts [mol].
      real(pr), intent(in) :: volume !! Total volume [L].
      real(pr), intent(in) :: temperature !! Temperature [K].
      real(pr), intent(out) :: values(16) !! Scalars, gradients, and flattened Hessian.
      real(pr) :: hessian(2, 2) !! Amount Hessian [bar L/mol^2].

      call eos%residual_helmholtz(amounts, volume, temperature, &
         Ar=values(1), ArV=values(2), ArT=values(3), ArTV=values(4), &
         ArV2=values(5), ArT2=values(6), Arn=values(7:8), &
         ArVn=values(9:10), ArTn=values(11:12), Arn2=hessian)
      values(13:16) = reshape(hessian, [4])
   end subroutine evaluate

   subroutine check_derivatives(eos)
      class(ArModel), intent(in) :: eos !! Translated EOS under test.
      real(pr) :: analytic(16) !! Analytic residual energy and derivatives.
      real(pr) :: numeric(16) !! Central finite-difference reference.
      real(pr) :: hessian(2, 2) !! Numerical amount Hessian [bar L/mol^2].
      real(pr) :: step !! Step-size multiplier.
      integer :: i !! Step-size case.

      call evaluate(eos, n, V, T, analytic)
      do i = 1, 2
         step = 1._pr/real(i, pr)
         call numeric_ar_derivatives(eos, n, V, T, &
            d_n=step*2.e-4_pr, d_v=step*1.e-4_pr, d_t=step*0.05_pr, &
            Ar=numeric(1), ArV=numeric(2), ArT=numeric(3), ArTV=numeric(4), &
            ArV2=numeric(5), ArT2=numeric(6), Arn=numeric(7:8), &
            ArVn=numeric(9:10), ArTn=numeric(11:12), Arn2=hessian)
         numeric(13:16) = reshape(hessian, [4])
         call assert(close(analytic, numeric, 3.e-5_pr), "All derivatives agree with energy finite differences")
      end do
   end subroutine check_derivatives

    subroutine check_optional_outputs(eos)
      class(ArModel), intent(in) :: eos !! Model under test.
      real(pr) :: reference(16) !! All-output reference evaluation.
      real(pr) :: individual(16) !! Outputs obtained one at a time, in the same order.
      real(pr) :: hessian(2, 2) !! Amount Hessian [bar L/mol^2].

      call evaluate(eos, n, V, T, reference)
      call eos%residual_helmholtz(n, V, T, Ar=individual(1))
      call eos%residual_helmholtz(n, V, T, ArV=individual(2))
      call eos%residual_helmholtz(n, V, T, ArT=individual(3))
      call eos%residual_helmholtz(n, V, T, ArTV=individual(4))
      call eos%residual_helmholtz(n, V, T, ArV2=individual(5))
      call eos%residual_helmholtz(n, V, T, ArT2=individual(6))
      call eos%residual_helmholtz(n, V, T, Arn=individual(7:8))
      call eos%residual_helmholtz(n, V, T, ArVn=individual(9:10))
      call eos%residual_helmholtz(n, V, T, ArTn=individual(11:12))
      call eos%residual_helmholtz(n, V, T, Arn2=hessian)
      individual(13:16) = reshape(hessian, [4])
      call assert(close(individual, reference, 1.e-12_pr), "Each individual output agrees with the all-output call")

      ! Common pressure and fugacity calls request several outputs together.
      call eos%residual_helmholtz(n, V, T, ArV=individual(2), ArV2=individual(5))
      call assert(close(individual([2, 5]), reference([2, 5]), 1.e-12_pr), &
         "Pressure derivative outputs can be requested together")
      call eos%residual_helmholtz(n, V, T, ArV=individual(2), Arn=individual(7:8))
      call assert(close(individual([2, 7, 8]), reference([2, 7, 8]), 1.e-12_pr), &
         "Fugacity outputs can be requested together")
      call eos%residual_helmholtz(n, V, T, ArTn=individual(11:12), Arn2=hessian)
      call assert(close(individual(11:12), reference(11:12), 1.e-12_pr) .and. &
         close(reshape(hessian, [4]), reference(13:16), 1.e-12_pr), &
         "Mixed temperature and amount Hessian outputs can be requested together")
      call eos%residual_helmholtz(n, V, T)
    end subroutine check_optional_outputs

    subroutine check_base_requests()
       type(TrackingResidualModel) :: tracked
       type(TestTranslation) :: translated
       real(pr) :: energy
       real(pr) :: volume_gradient
       real(pr) :: volume_hessian
       real(pr) :: pressure
       real(pr) :: pressure_volume_gradient
       real(pr) :: pressure_temperature_gradient
       real(pr) :: pressure_amount_gradient(2)
       real(pr) :: chemical_potential(2)
       real(pr) :: fugacity(2)
       real(pr) :: temperature_hessian
       real(pr) :: amount_volume(2)

       tracked%components = base%components
       call translated%init(tracked)
       call translated%residual_helmholtz(n, V, T, Ar=energy)
       call assert(all(base_requested .eqv. [.true., .false., .false., .false., .false., &
          .false., .false., .false., .false., .false.]), "Energy requests no base derivatives")

       call translated%residual_helmholtz(n, V, T, ArV=volume_gradient, ArV2=volume_hessian)
       call assert(all(base_requested .eqv. [.false., .true., .false., .false., .true., &
          .false., .false., .false., .false., .false.]), "Volume derivatives request only base volume derivatives")

       call translated%residual_helmholtz(n, V, T, ArV=volume_gradient, Arn=chemical_potential)
       call assert(all(base_requested .eqv. [.false., .true., .false., .false., .false., &
          .false., .true., .false., .false., .false.]), "Flash fugacity requests only base ArV and Arn")

       call translated%lnphi_vt(n, V, T, lnPhi=fugacity)
       call assert(all(base_requested .eqv. [.false., .true., .false., .false., .false., &
          .false., .true., .false., .false., .false.]), "Inherited fugacity uses the lean base call")

       call translated%pressure(n, V, T, pressure)
       call assert(all(base_requested .eqv. [.false., .true., .false., .false., .false., &
          .false., .false., .false., .false., .false.]), "Pressure alone requests only the base volume gradient")

       call translated%pressure(n, V, T, pressure, dPdV=pressure_volume_gradient)
       call assert(all(base_requested .eqv. [.false., .true., .false., .false., .true., &
          .false., .false., .false., .false., .false.]), "Pressure volume gradient requests only base ArV2")

       call translated%pressure(n, V, T, pressure, dPdT=pressure_temperature_gradient)
       call assert(all(base_requested .eqv. [.false., .true., .false., .true., .true., &
          .false., .false., .false., .false., .false.]), "Pressure temperature gradient requests its chain-rule inputs")

       call translated%pressure(n, V, T, pressure, dPdn=pressure_amount_gradient)
       call assert(all(base_requested .eqv. [.false., .true., .false., .false., .true., &
          .false., .false., .true., .false., .false.]), "Pressure amount gradient avoids temperature derivatives")

       call translated%residual_helmholtz(n, V, T, ArT2=temperature_hessian, ArVn=amount_volume)
       call assert(all(base_requested .eqv. [.false., .true., .false., .true., .true., &
          .true., .false., .true., .false., .false.]), "Combined temperature and amount derivatives retain dependencies")
    end subroutine check_base_requests

   subroutine check_roots(eos)
      type(TestTranslation), intent(in) :: eos !! Translation with a cubic base model.
      character(len=6), parameter :: roots(3) = ["liquid", "vapor ", "stable"]
      !! Base root labels tested without changing their meaning.
      real(pr) :: shift !! Volume shift [L].
      real(pr) :: base_volume !! Root of the unmodified EOS [L].
      real(pr) :: volume !! Corresponding physical root [L].
      real(pr) :: pressure !! Pressure reconstructed from the translated energy [bar].
      real(pr) :: temperature !! Subcritical test temperature [K].
      real(pr) :: specified_pressure !! Specified pressure [bar].
      integer :: i !! Root-label index.
      integer :: j !! State index: vapor-stable and liquid-stable examples.

      do j = 1, 2
         temperature = 200._pr
         specified_pressure = 1._pr
         if (j == 2) then
            temperature = 150._pr
            specified_pressure = 2._pr
         end if
         call eos%translation(n, temperature, shift)
         do i = 1, 3
            call eos%base%volume(n, specified_pressure, temperature, base_volume, trim(roots(i)))
            call eos%volume(n, specified_pressure, temperature, volume, trim(roots(i)))
            call assert(volume > 0._pr, "Translated root remains physical")
            call assert(abs(volume - (base_volume-shift)) < 1.e-12_pr*base_volume, &
               "Liquid, vapor and stable volume roots delegate to the base solver")
            call eos%pressure(n, volume, temperature, pressure)
            call assert(abs(pressure-specified_pressure) < 1.e-6_pr*specified_pressure, &
               "Delegated volume root satisfies the translated pressure equation")
         end do
         call assert(abs(eos%get_v0(n, specified_pressure, temperature) - &
            (eos%base%get_v0(n, specified_pressure, temperature)-shift)) < 1.e-14_pr, &
            "Liquid-volume initializer is shifted in total-volume units")
      end do
   end subroutine check_roots

   subroutine check_flash(eos)
      type(TestTranslation), intent(in) :: eos !! Constant-per-component translated cubic.
      type(EquilibriumState) :: reference !! Base two-phase equilibrium.
      type(EquilibriumState) :: translated !! Equilibrium from the adapter via the unchanged flash API.
      real(pr) :: shift_liquid !! Molar translation in the liquid phase [L/mol].
      real(pr) :: shift_vapor !! Molar translation in the vapor phase [L/mol].
      real(pr) :: liquid_phi(2) !! Logarithms of liquid fugacity coefficients.
      real(pr) :: vapor_phi(2) !! Logarithms of vapor fugacity coefficients.

      reference = flash(eos%base, n/sum(n), 200._pr, p_spec=10._pr)
      translated = flash(eos, n/sum(n), 200._pr, p_spec=10._pr)
      call assert(reference%kind == "split" .and. translated%kind == "split", &
         "Existing flash algorithm accepts the adapter as an ArModel")
      call assert(close(translated%x, reference%x, 1.e-6_pr) .and. &
         close(translated%y, reference%y, 1.e-6_pr), "Linear volume translation preserves TP phase compositions")
      call assert(abs(translated%beta-reference%beta) < 1.e-6_pr, "Linear translation preserves the TP phase fraction")
      call eos%translation(translated%x, translated%T, shift_liquid)
      call eos%translation(translated%y, translated%T, shift_vapor)
      call assert(abs(translated%Vx - (reference%Vx-shift_liquid)) < 1.e-7_pr .and. &
         abs(translated%Vy - (reference%Vy-shift_vapor)) < 1.e-7_pr, "Flash returns physical translated phase volumes")
      call assert(maxval(abs(n/sum(n) - &
         ((1._pr-translated%beta)*translated%x + translated%beta*translated%y))) < 1.e-8_pr, &
         "Translated flash satisfies component material balances")
      call eos%lnphi_vt(translated%x, translated%Vx, translated%T, lnPhi=liquid_phi)
      call eos%lnphi_vt(translated%y, translated%Vy, translated%T, lnPhi=vapor_phi)
      call assert(maxval(abs(log(translated%x)+liquid_phi-log(translated%y)-vapor_phi)) < 2.e-6_pr, &
         "Translated flash satisfies phase isofugacity")
   end subroutine check_flash

   subroutine check_consistency(eos)
      class(ArModel), intent(in) :: eos !! Translated model used through the ordinary ArModel contract.
      real(pr) :: eq31 !! Gibbs/fugacity identity residual.
      real(pr) :: eq33(2, 2) !! Fugacity amount-derivative symmetry residual.
      real(pr) :: eq34(2) !! Gibbs-Duhem residual.
      real(pr) :: eq36 !! Pressure derivative identity residual.
      real(pr) :: eq37 !! Enthalpy/temperature derivative identity residual.

      call ar_consistency(eos, n, V, T, eq31, eq33, eq34, eq36, eq37)
      call assert(maxval(abs([eq31, reshape(eq33, [4]), eq34, eq36, eq37])) < 1.e-12_pr, &
         "Translated model satisfies residual thermodynamic consistency identities")
   end subroutine check_consistency

   subroutine check_scaling(eos)
      type(TestTranslation), intent(in) :: eos !! Nonlinear extensive translation.
      real(pr) :: values(16) !! Energy and derivatives on the original basis.
      real(pr) :: scaled(16) !! Energy and derivatives on a scaled amount/volume basis.
      real(pr) :: weights(16) !! Expected factors for each derivative under scaling.
      real(pr), parameter :: factor = 3.7_pr !! Arbitrary noninteger scaling of extensive inputs.

      call evaluate(eos, n, V, T, values)
      call evaluate(eos, factor*n, factor*V, T, scaled)
      weights = 1._pr
      weights([1, 3, 6]) = factor
      weights([5, 9, 10, 13, 14, 15, 16]) = 1._pr/factor
      call assert(close(scaled, weights*values, 1.e-11_pr), "Energy and derivatives have the correct extensivity")
      call evaluate(eos, [0._pr, sum(n)], V, T, scaled)
      call assert(all(abs(scaled) < huge(1._pr)), "Zero amount of one component has finite residual derivatives")
   end subroutine check_scaling

   subroutine check_pressure_and_calorics(eos)
      type(TestTranslation), intent(in) :: eos !! Temperature-dependent translation.
      real(pr) :: shift !! Total volume shift [L].
      real(pr) :: shiftT !! Translation temperature derivative [L/K].
      real(pr) :: shiftT2 !! Translation second temperature derivative [L/K^2].
      real(pr) :: base_pressure !! Pressure at the corresponding base state [bar].
      real(pr) :: pressure !! Translated pressure [bar].
      real(pr) :: base_H !! Base residual enthalpy at the same P,T [bar L].
      real(pr) :: translated_H !! Translated residual enthalpy at the same P,T [bar L].
      real(pr) :: base_Cp !! Base residual heat capacity [bar L/K].
      real(pr) :: translated_Cp !! Translated residual heat capacity [bar L/K].
      real(pr) :: base_lnphi(2) !! Base fugacity coefficients in logarithmic form.
      real(pr) :: lnphi(2) !! Translated fugacity coefficients in logarithmic form.
      real(pr) :: shift_n(2) !! Translation amount gradient [L/mol].

      call eos%translation(n, T, shift, CT=shiftT, CT2=shiftT2, Cn=shift_n)
      call eos%base%pressure(n, V+shift, T, base_pressure)
      call eos%pressure(n, V, T, pressure)
      call assert(abs(pressure-base_pressure) < 1.e-12_pr, "Pressure is exactly the base pressure at V+C")
      call eos%base%lnphi_vt(n, V+shift, T, lnPhi=base_lnphi)
      call eos%lnphi_vt(n, V, T, lnPhi=lnphi)
      call assert(close(lnphi, base_lnphi-pressure*shift_n/(R*T), 1.e-11_pr), &
         "Fugacities include the composition-dependent translation")
      call eos%base%enthalpy_residual_vt(n, V+shift, T, Hr=base_H)
      call eos%enthalpy_residual_vt(n, V, T, Hr=translated_H)
      call assert(abs(translated_H - (base_H-pressure*(shift-T*shiftT))) < 1.e-11_pr, &
         "Inherited residual enthalpy includes the temperature-dependent correction")
      call eos%base%Cp_residual_vt(n, V+shift, T, base_Cp)
      call eos%Cp_residual_vt(n, V, T, translated_Cp)
      call assert(abs(translated_Cp - (base_Cp+pressure*T*shiftT2)) < 1.e-12_pr, &
         "Inherited residual Cp includes the second temperature derivative")
   end subroutine check_pressure_and_calorics

   subroutine check_non_cubic_and_failures()
      type(ZeroResidualModel) :: zero !! Non-cubic test base with configurable failure outputs.
      type(TestTranslation) :: translated !! Adapter around the non-cubic base.
      type(TestTranslation) :: nested !! Adapter around an already translated model.
      real(pr) :: shift !! Volume translation [L].
      real(pr) :: volume !! Returned root [L].
      real(pr) :: energy !! Translated residual energy [bar L].
      real(pr) :: initializer !! Translated initial volume [L].
      real(pr) :: invalid !! Quiet NaN used to test failure adaptation.

      invalid = ieee_value(0._pr, ieee_quiet_nan)
      call translated%init(zero)
      call assert(translated%name == "VolumeTranslatedModel", "Unnamed bases receive a default adapter name")
      call translated%translation(n, T, shift)
      call translated%residual_helmholtz(n, V, T, Ar=energy)
      call assert(abs(energy-sum(n)*R*T*log(V/(V+shift))) < 1.e-13_pr, &
         "A non-cubic zero-residual base gives the independent logarithmic reference")
      call check_derivatives(translated)
      call translated%volume(n, 1._pr, T, volume, "stable")
      call assert(abs(volume-(zero%root_volume-shift)) < 1.e-14_pr, &
         "Adapter delegates even when the base has a custom non-cubic volume solver")

      call nested%init(translated)
      call nested%volume(n, 1._pr, T, volume, "stable")
      call assert(abs(volume-(zero%root_volume-2._pr*shift)) < 1.e-14_pr, &
         "Already translated models can be wrapped without slicing their dynamic type")
      call nested%residual_helmholtz(n, V, T, Ar=energy)
      call assert(abs(energy-sum(n)*R*T*log(V/(V+2._pr*shift))) < 1.e-13_pr, &
         "Nested reference corrections compose consistently")

      ! A large negative shift must never turn a failed base root into success.
      zero%root_volume = -1._pr
      zero%initial_volume = -1._pr
      translated%c0 = -10._pr
      call translated%init(zero)
      call translated%volume(n, 1._pr, T, volume, "stable")
      initializer = translated%get_v0(n, 1._pr, T)
      call assert(volume == -1._pr .and. initializer == -1._pr, "Base failure sentinels are not translated")

      zero%root_volume = invalid
      zero%initial_volume = invalid
      call translated%init(zero)
      call translated%volume(n, 1._pr, T, volume, "stable")
      call assert(volume == -1._pr .and. translated%get_v0(n, 1._pr, T) == -1._pr, &
         "Nonfinite base volumes become explicit failures")

      zero%root_volume = 1._pr
      zero%initial_volume = 0.2_pr
      translated%c0 = 10._pr
      call translated%init(zero)
      call translated%volume(n, 1._pr, T, volume, "stable")
      call assert(volume == -1._pr .and. translated%get_v0(n, 1._pr, T) == -1._pr, &
         "Nonpositive translated roots and initializers are rejected")

      translated%c0 = 0.12_pr
      call translated%volume(n, 1._pr, T, volume, "stable")
      call assert(volume > 0._pr .and. translated%get_v0(n, 1._pr, T) == -1._pr, &
         "A rejected initializer does not prevent delegation to a valid base volume root")

      translated%c0 = invalid
      call translated%volume(n, 1._pr, T, volume, "stable")
      call assert(volume == -1._pr .and. translated%get_v0(n, 1._pr, T) == -1._pr, &
         "Nonfinite translations produce explicit volume failures")
   end subroutine check_non_cubic_and_failures

   subroutine reject_invalid_state(scenario)
      character(len=*), intent(in) :: scenario !! Expected error-stop branch to exercise externally.
      type(TestTranslation) :: translated !! Initialized adapter used for invalid-domain tests.
      type(ZeroResidualModel) :: zero !! Valid base used for invalid-domain tests.
      real(pr) :: value !! Requested result, which must not be returned.

      select case (scenario)
       case ("invalid-physical-volume")
         call translated%init(zero)
         call translated%residual_helmholtz(n, 0._pr, T, Ar=value)
       case ("invalid-base-volume")
         call translated%init(zero)
         translated%c0 = -10._pr
         call translated%residual_helmholtz(n, V, T, Ar=value)
       case ("nonfinite-physical-volume")
         call translated%init(zero)
         call translated%residual_helmholtz(n, ieee_value(0._pr, ieee_quiet_nan), T, Ar=value)
       case ("nonfinite-base-volume")
         call translated%init(zero)
         translated%c0 = ieee_value(0._pr, ieee_quiet_nan)
         call translated%residual_helmholtz(n, V, T, Ar=value)
       case default
         error stop "Unknown expected-failure scenario"
      end select
      stop 0
   end subroutine reject_invalid_state

end program test_volume_translated
