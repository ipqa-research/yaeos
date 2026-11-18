module yaeos__models_ar_cubic_implementations
   !! # Implemented cubic equations of state
   !!
   !! # Description
   !! Constructors for the cubic equations implemented in `yaeos`:
   !!
   !! - [[PengRobinson76]]
   !! - [[PengRobinson78]]
   !! - [[SoaveRedlichKwong]]
   !! - [[RKPR]]
   !! - [[PSRK]]
   !! Each constructor returns a configured [[CubicEoS]] with its standard
   !! alpha function, critical-property correlation, and default mixing rule.
   !!
   !! # Examples
   !! ```fortran
   !! model = PengRobinson76(Tc, Pc, omega)
   !! ```
   use yaeos__constants, only: pr, R
   use yaeos__models_ar_genericcubic, only: CubicEoS
   use yaeos__substance, only: Substances
   implicit none

   private

   public :: PengRobinson76
   public :: PengRobinson78
   public :: SoaveRedlichKwong
   public :: RKPR
   public :: PSRK
   public :: refit_rkpr_k

contains

   type(CubicEoS) function PengRobinson76(tc, pc, w, kij, lij) result(model)
      !! # PengRobinson76
      !!
      !! # Description
      !! PengRobinson76.
      !!
      !! Using the critical constants setup the parameters to use the
      !! PengRobinson Equation of State
      !!
      !! - \[\alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2\]
      !! - \[k = 0.37464 + 1.54226 \omega - 0.26993 \omega^2 \]
      !! - \[a_c = 0.45723553  R^2 T_c^2 / P_c\]
      !! - \[b = 0.07779607  R T_c/P_c\]
      !! - \[\delta_1 = 1 + \sqrt{2}\]
      !! - \[\delta_2 = 1 - \sqrt{2}\]
      !!
      !! There is also the optional possibility to include the \(k_{ij}\) and
      !! \(l_{ij}\) matrices. Using by default Classic Van der Waals mixing
      !! rules.
      !!
      !! After setting up the model, it is possible to redefine either the
      !! mixing rule or the alpha function using a different derived type
      !! defined outside the function.
      use yaeos__constants, only: pr, R
      use yaeos__substance, only: Substances
      use yaeos__models_ar_genericcubic, only: CubicEoS
      use yaeos__models_ar_cubic_alphas, only: AlphaSoave
      use yaeos__models_ar_cubic_quadratic_mixing, only: QMR
      real(pr), intent(in) :: tc(:) !! Critical Temperatures [K]
      real(pr), intent(in) :: pc(:) !! Critical Pressures [bar]
      real(pr), intent(in) :: w(:) !! Acentric Factors
      real(pr), optional, intent(in) :: kij(:, :) !! \(k_{ij}\) matrix
      real(pr), optional, intent(in) :: lij(:, :) !! \(l_{ij}\) matrix

      type(Substances) :: composition
      type(QMR) :: mixrule
      type(AlphaSoave) :: alpha !! Soave alpha function
      integer :: nc !! Number of components.
      integer :: i !! Loop/index variable.

      nc = size(tc)

      composition%tc = tc
      composition%pc = pc
      composition%w = w

      alpha%k = 0.37464_pr &
         + 1.54226_pr * composition%w &
         - 0.26993_pr * composition%w**2

      if (present(kij)) then
         mixrule%k = kij
      else
         mixrule%k = reshape([(0, i=1,nc**2)], [nc, nc])
      endif

      if (present(lij)) then
         mixrule%l = lij
      else
         mixrule%l = reshape([(0, i=1,nc**2)], [nc, nc])
      endif

      model%components = composition
      model%ac = 0.45723553_pr * R**2 * composition%tc**2 / composition%pc
      model%b = 0.07779607_pr * R * composition%tc/composition%pc
      model%del1 = [(1 + sqrt(2.0_pr), i=1,nc)]
      model%del2 = [(1 - sqrt(2.0_pr), i=1,nc)]
      model%alpha = alpha
      model%mixrule = mixrule
      model%name = "PR76"
   end function PengRobinson76

   type(CubicEoS) function PengRobinson78(tc, pc, w, kij, lij) result(model)
      !! # PengRobinson78
      !!
      !! # Description
      !! PengRobinson78.
      !!
      !! Using the critical constants setup the parameters to use the
      !! PengRobinson Equation of State
      !!
      !! - \[\alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2\]
      !! - \[k = 0.37464 + 1.54226 \omega - 0.26992 \omega^2  \text{ where } \omega <=0.491\]
      !! - \[k = 0.379642 + 1.48503 \omega - 0.164423 \omega^2  + 0.016666 \omega^3 \text{ where } \omega > 0.491\]
      !! - \[a_c = 0.45723553  R^2 T_c^2 / P_c\]
      !! - \[b = 0.07779607  R T_c/P_c\]
      !! - \[\delta_1 = 1 + \sqrt{2}\]
      !! - \[\delta_2 = 1 - \sqrt{2}\]
      !!
      !! There is also the optional possibility to include the \(k_{ij}\) and
      !! \(l_{ij}\) matrices. Using by default Classic Van der Waals mixing
      !! rules.
      !!
      !! After setting up the model, it is possible to redefine either the
      !! mixing rule or the alpha function using a different derived type
      !! defined outside the function.
      use yaeos__constants, only: pr, R
      use yaeos__substance, only: Substances
      use yaeos__models_ar_genericcubic, only: CubicEoS
      use yaeos__models_ar_cubic_alphas, only: AlphaSoave
      use yaeos__models_ar_cubic_quadratic_mixing, only: QMR
      real(pr), intent(in) :: tc(:) !! Critical Temperatures [K]
      real(pr), intent(in) :: pc(:) !! Critical Pressures [bar]
      real(pr), intent(in) :: w(:) !! Acentric Factors
      real(pr), optional, intent(in) :: kij(:, :) !! \(k_{ij}\) matrix
      real(pr), optional, intent(in) :: lij(:, :) !! \(l_{ij}\) matrix

      type(Substances) :: composition
      type(QMR) :: mixrule
      type(AlphaSoave) :: alpha !! Soave alpha function
      integer :: nc !! Number of components.
      integer :: i !! Loop/index variable.

      nc = size(tc)

      composition%tc = tc
      composition%pc = pc
      composition%w = w

      allocate(alpha%k(nc))
      where (composition%w <=0.491)
         alpha%k = 0.37464 + 1.54226 * composition%w - 0.26992 * composition%w**2
      elsewhere
         alpha%k = 0.379642 + 1.48503 * composition%w - 0.164423 * composition%w**2 + 0.016666 * composition%w**3
      end where

      if (present(kij)) then
         mixrule%k = kij
      else
         mixrule%k = reshape([(0, i=1,nc**2)], [nc, nc])
      endif

      if (present(lij)) then
         mixrule%l = lij
      else
         mixrule%l = reshape([(0, i=1,nc**2)], [nc, nc])
      endif

      model%components = composition
      model%ac = 0.45723553_pr * R**2 * composition%tc**2 / composition%pc
      model%b = 0.07779607_pr * R * composition%tc/composition%pc
      model%del1 = [(1 + sqrt(2.0_pr), i=1,nc)]
      model%del2 = [(1 - sqrt(2.0_pr), i=1,nc)]
      model%alpha = alpha
      model%mixrule = mixrule
      model%name = "PR78"
   end function PengRobinson78

   type(CubicEoS) function SoaveRedlichKwong(tc, pc, w, kij, lij) result(model)
      !! # SoaveRedlichKwong
      !!
      !! # Description
      !! SoaveRedlichKwong.
      !!
      !! Using the critical constants setup the parameters to use the
      !! SoaveRedlichKwong Equation of State
      !!
      !! - \[\alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2\]
      !! - \[k = 0.48 + 1.574 \omega - 0.175 \omega^2 \]
      !! - \[a_c = 0.427480  R^2 T_c^2/P_c\]
      !! - \[b = 0.086640  R T_c/P_c\]
      !! - \[\delta_1 = 1\]
      !! - \[\delta_2 = 0\]
      !!
      !! There is also the optional possibility to include the \(k_{ij}\) and \(l_{ij}\)
      !! matrices. Using by default Classic Van der Waals mixing rules.
      !!
      !! After setting up the model, it is possible to redefine either the
      !! mixing rule or the alpha function using a different derived type
      !! defined outside the function.
      use yaeos__models_ar_genericcubic, only: CubicEoS
      use yaeos__models_ar_cubic_alphas, only: AlphaSoave
      use yaeos__models_ar_cubic_quadratic_mixing, only: QMR
      real(pr), intent(in) :: tc(:) !! Critical temperature [K]
      real(pr), intent(in) :: pc(:) !! Critical pressure [bar]
      real(pr), intent(in) :: w(:) !! Acentric factor
      real(pr), optional, intent(in) :: kij(:, :) !! \(k_{ij}\) matrix
      real(pr), optional, intent(in) :: lij(:, :) !! \(l_{ij}\) matrix

      type(Substances) :: composition
      type(QMR) :: mixrule
      type(AlphaSoave) :: alpha !! Soave alpha function
      integer :: nc !! Number of components.
      integer :: i !! Loop/index variable.

      nc = size(tc)

      composition%tc = tc
      composition%pc = pc
      composition%w = w

      alpha%k = 0.48_pr + 1.574_pr * composition%w - 0.175_pr * composition%w**2

      if (present(kij)) then
         mixrule%k = kij
      else
         mixrule%k = reshape([(0, i=1,nc**2)], [nc, nc])
      endif

      if (present(lij)) then
         mixrule%l = lij
      else
         mixrule%l = reshape([(0, i=1,nc**2)], [nc, nc])
      endif

      model%components = composition
      model%ac = 0.427480_pr * R**2 * composition%tc**2/composition%pc
      model%b = 0.086640_pr * R * composition%tc/composition%pc
      model%del1 = [(1, i=1,nc)]
      model%del2 = [(0, i=1,nc)]
      model%alpha = alpha
      model%mixrule = mixrule
      model%name = "SRK"
   end function SoaveRedlichKwong

   type(CubicEoS) function PSRK(tc, pc, w, molecules, c1, c2, c3) result(model)
      !! # PSRK
      !!
      !! # Description
      !! Predictive Soave-Redlich-Kwong (PSRK).
      !!
      !! SRK equation combined with the modified Huron-Vidal mixing rule at
      !! \(q = -0.64663\) and a PSRK UNIFAC excess Gibbs model. The alpha
      !! function is Soave by default, or Mathias-Copeman when `c1`, `c2`, and
      !! `c3` are provided.
      use yaeos__models_ar_genericcubic, only: CubicEoS
      use yaeos__models_ar_cubic_alphas, only: AlphaMathiasCopeman, AlphaSoave
      use yaeos__models_cubic_mixing_rules_huron_vidal, only: MHV
      use yaeos__models_ge_implementations, only: setup_psrk, UNIFAC
      use yaeos__models_ge_group_contribution_groups, only: Groups
      real(pr), intent(in) :: tc(:) !! Critical temperature [K]
      real(pr), intent(in) :: pc(:) !! Critical pressure [bar]
      real(pr), intent(in) :: w(:) !! Acentric factor
      type(Groups), intent(in) :: molecules(:) !! UNIFAC group assignment
      real(pr), optional, intent(in) :: c1(:)
      real(pr), optional, intent(in) :: c2(:)
      real(pr), optional, intent(in) :: c3(:)
      !! Mathias-Copeman alpha parameters

      type(UNIFAC) :: ge
      type(Substances) :: composition
      type(MHV) :: mixrule
      type(AlphaSoave) :: alpha !! Soave alpha function
      type(AlphaMathiasCopeman) :: alpha_mc
      integer :: nc !! Number of components.
      integer :: i !! Loop/index variable.

      nc = size(tc)

      composition%tc = tc
      composition%pc = pc
      composition%w = w

      ge = setup_psrk(molecules)

      if (present(c1) .and. present(c2) .and. present(c3)) then
         alpha_mc = AlphaMathiasCopeman(c1, c2, c3)
         model%alpha = alpha_mc
      else
         alpha%k = 0.48_pr + 1.574_pr * composition%w - 0.175_pr * composition%w**2
         model%alpha = alpha
      end if

      model%components = composition
      model%ac = 0.427480_pr * R**2 * composition%tc**2/composition%pc
      model%b = 0.086640_pr * R * composition%tc/composition%pc
      model%del1 = [(1, i=1,nc)]
      model%del2 = [(0, i=1,nc)]

      mixrule = MHV(ge=ge, b=model%b, q=-0.64663_pr)
      model%mixrule = mixrule

      model%name = "PSRK"
   end function PSRK

   type(CubicEoS) function RKPR(tc, pc, w, zc, kij, lij, delta_1, k) result(model)
      !! # RKPR
      !!
      !! # Description
      !! RKPR Equation of State
      !!
      !! The RKPR EoS extends the classical formulation of Cubic Equations
      !! of State by freeing the parameter \(\delta_1\). This extra degree
      !! provides extra ways of implementing the equation in comparison
      !! of other Cubic EoS (like PR and SRK) which are limited to definition
      !! of their critical constants.
      !!
      !! Besides that extra parameter, the RKPR includes another \(\alpha\)
      !! function:
      !! \[
      !!  \alpha(T_r) = \left(\frac{3}{2+T_r}\right)^k
      !! \]
      !!
      !! In this implementation we take the simplest form which correlates
      !! the extra parameter to the critical compressibility factor \(Z_c\) and
      !! the \(k\) parameter of the \(\alpha\) function to \(Z_c\) and \(\omega\):
      !!
      !! \[\delta_1 = d_1 + d_2 (d_3 - Z_c^{eos})^{d_4} + d_5 (d_3 - Z_c^{eos})^{d_6}\]
      !! \[k = (A_1 Z_c^{eos} + A_0)\omega^2 + (B_1 Z_c + B_0)\omega + (C_1 Z_c^{eos} + C_0)\]
      !!
      !! where \(Z_c^{eos} = 1.168 Z_c\). When `delta_1` or `k` are given they
      !! replace the correlations. When `k` is not given, each component's
      !! \(k\) is afterwards refitted with `refit_rkpr_k` so that the model
      !! reproduces the acentric factor.
      use yaeos__models_ar_cubic_quadratic_mixing, only: QMR
      use yaeos__models_ar_cubic_alphas, only: AlphaRKPR
      real(pr), intent(in) :: tc(:) !! Critical Temperature [K]
      real(pr), intent(in) :: pc(:) !! Critical Pressure [bar]
      real(pr), intent(in) :: w(:) !! Acentric Factor
      real(pr), intent(in) :: zc(:) !! Critical compressibility
      real(pr), optional, intent(in) :: kij(:, :) !! \(k_{ij}\) matrix
      real(pr), optional, intent(in) :: lij(:, :) !! \(l_{ij}\) matrix
      real(pr), optional, intent(in) :: delta_1(:) !! Pure-component \(\delta_1\) values (overrides the correlation)
      real(pr), optional, intent(in) :: k(:) !! Alpha-function \(k\) values (overrides correlation and refit)

      type(AlphaRKPR) :: alpha !! RKPR alpha function
      type(QMR) :: mixrule
      type(Substances) :: composition

      integer :: i !! Loop/index variable.
      integer :: nc !! Number of components.

      real(pr), parameter :: d1 = 0.428364 !! Constant term in the RKPR delta-one correlation.
      real(pr), parameter :: d2 = 18.496215 !! Amplitude of the first compressibility-dependent term.
      real(pr), parameter :: d3=0.338426 !! Reference compressibility in both power-law terms.
      real(pr), parameter :: d4=0.66 !! Exponent of the first compressibility-dependent term.
      real(pr), parameter :: d5 = 789.723105 !! Amplitude of the second compressibility-dependent term.
      real(pr), parameter :: d6=2.512392 !! Exponent of the second compressibility-dependent term.

      real(pr), parameter :: A1 = -2.4407
      real(pr), parameter :: A0 = 0.0017
      real(pr), parameter :: B1 =7.4513
      real(pr), parameter :: B0 =1.9681
      real(pr), parameter :: C1 =12.504
      real(pr), parameter :: C0 =-2.6238

      real(pr) :: OMa(size(pc))
      real(pr) :: OMb(size(pc))
      real(pr) :: Zc_eos(size(pc))
      real(pr) :: Psat_i
      real(pr) :: diff

      nc = size(tc)

      composition%pc = pc
      composition%tc = tc
      composition%w = w

      Zc_eos = 1.168 * Zc

      if (present(k)) then
         alpha%k = k
      else
         alpha%k = (A1 * Zc_eos + A0)*w**2 + (B1*zc + B0)*w + (C1*Zc_eos + C0)
      end if

      if (present(kij)) then
         mixrule%k = kij
      else
         mixrule%k = reshape([(0, i=1,nc**2)], [nc, nc])
      end if

      if (present(lij)) then
         mixrule%l = lij
      else
         mixrule%l = reshape([(0, i=1,nc**2)], [nc, nc])
      end if

      model%components = composition
      if (present(delta_1)) then
         model%del1 = delta_1
      else
         model%del1 = d1 + d2 * (d3 - zc_eos) ** d4 + d5 * (d3 - zc_eos) ** d6
      end if

      model%del2 = (1._pr - model%del1)/(1._pr + model%del1)
      model%alpha = alpha

      call get_OMa_OMb(model%del1, oma, omb)
      model%ac = OMa * (R*Tc)**2/Pc
      model%b = OMb * (R*Tc)/Pc

      model%mixrule = mixrule
      model%name = "RKPR 2005"

      if (.not. present(k)) then
         do i=1,nc
            call refit_rkpr_k(model, i)
         end do
      end if
   end function RKPR


   subroutine refit_rkpr_k(model, component)
      !! # refit_rkpr_k
      !!
      !! Refit the \(k\) parameter of the RKPR EoS to match the acentric
      !! factor.
      !!
      !! # Description
      !! Uses a secant iteration on \(k\) so that the model's
      !! \(\omega = -1 - \log_{10}(P^{sat}(0.7T_c)/P_c)\) matches the stored
      !! acentric factor of `component`.
      use yaeos__models_ar_cubic_alphas, only: AlphaRKPR
      use yaeos__equilibria_boundaries_pure_saturation, only: &
         PurePsat, pure_saturation_line

      type(CubicEoS), intent(in out) :: model !! The model to be refitted
      integer, intent(in) :: component !! Component index to refit

      real(pr) :: Psat_i
      real(pr) :: Tc
      real(pr) :: Pc
      real(pr) :: w !! Experimental acentric factor of the component
      real(pr) :: Vv !! Volume of vapor [L]
      real(pr) :: VL !! Volume of liquid [L]
      real(pr) :: w_eos
      real(pr) :: f_prev
      real(pr) :: f_curr
      real(pr) :: k_prev
      real(pr) :: k_curr
      real(pr) :: k_next
      real(pr) :: step

      integer :: i !! Loop/index variable.
      integer :: iter !! Iteration count.
      type(AlphaRKPR) :: alpha !! RKPR alpha function
      type(PurePsat) :: Psat

      associate(a => model%alpha)
         select type(a)
          type is (AlphaRKPR)
            alpha = a
         end select
      end associate

      i = component

      Tc = model%components%tc(i)
      Pc = model%components%pc(i)
      w = model%components%w(i)

      ! Solve w_eos(k) - w = 0 using the secant method.
      alpha%k(i) = 5
      k_curr = alpha%k(i)
      k_prev = k_curr + max(abs(k_curr)*1e-3_pr, 1e-4_pr)
      f_prev = huge(1._pr)

      iter = 0
      do while(abs(f_prev) > 1e-3)
         iter = iter + 1
         Psat_i = model%Psat_pure(i, 0.7*Tc, Vl=Vl, Vv=Vv)

         if (Psat_i < 1e-10 .or. isnan(Psat_i)) then
            ! If the saturation pressure did not converge, calculate the
            ! whole line
            Psat = pure_saturation_line(model, i, minp=0._pr, minT=0.6*Tc)
            Psat_i = Psat%get_P(0.7*Tc)

            ! If the saturation pressure still not converges, just give it
            ! a value of 80% of the critical pressure.
            ! TODO: This could be improved by using a better initial guess
            if (Psat_i < 1e-10) Psat_i = 0.1*Pc
         end if

         w_eos = (-1 - log10(Psat_i/Pc))
         f_curr = (w_eos - w)/w

         if (iter == 1) then
            ! The perturbed point is evaluated on the next iteration.
            k_prev = k_curr
            f_prev = f_curr
            k_curr = k_curr + min(abs(k_curr)*1e-3_pr, 1e-4_pr)
         else
            ! if (abs(f_curr - f_prev) <= epsilon(1._pr)) exit
            step = f_curr*(k_curr-k_prev)/(f_curr-f_prev)
            k_next = k_curr - sign(min(abs(step), 0.1*k_curr), step)
            k_prev = k_curr
            f_prev = f_curr
            k_curr = k_next
         end if
         if (iter > 5000) exit

         deallocate(model%alpha)
         alpha%k(i) = k_curr
         model%alpha = alpha
      end do

      if (allocated(model%alpha)) deallocate(model%alpha)
      alpha%k(i) = k_curr
      model%alpha = alpha
   end subroutine refit_rkpr_k

   subroutine get_OMa_OMb(del1, OMa, OMb)
      !! # get_OMa_OMb
      !!
      !! # Description
      !! Universal critical constants for a cubic equation with a free
      !! \(\delta_1\) parameter.
      !!
      !! Returns the \(\Omega_a\) and \(\Omega_b\) parameters that make the
      !! equation satisfy the critical-point conditions for each value of
      !! `del1`.
      real(pr), intent(in) :: del1(:) !! \(\delta_1\) parameter
      real(pr), intent(out) :: OMa(size(del1)) !! \(\Omega_a\)
      real(pr), intent(out) :: OMb(size(del1)) !! \(\Omega_b\)

      real(pr) :: d1(size(del1))
      real(pr) :: y(size(del1)) !! Auxiliary root \(y = 1 + (2(1+\delta_1))^{1/3} + (4/(1+\delta_1))^{1/3}\)

      d1 = (1._pr + del1**2._pr)/(1._pr + del1)
      y = 1._pr + (2._pr*(1._pr + del1))**(1.0_pr/3._pr) + (4._pr/(1._pr + del1))**(1.0_pr/3)
      OMa = (3._pr*y*y + 3._pr*y*d1 + d1**2._pr + d1 - 1.0_pr)/(3._pr*y + d1 - 1.0_pr)**2._pr
      OMb = 1._pr/(3._pr*y + d1 - 1.0_pr)
   end subroutine get_OMa_OMb
end module yaeos__models_ar_cubic_implementations
