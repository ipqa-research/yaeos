module yaeos__models_ar_cubic_implementations
   use yaeos__constants, only: pr, R
   use  yaeos__models_ar_genericcubic, only: CubicEoS
   use yaeos__substance, only: Substances
   implicit none
   !! Implemented Cubic Equations of State.
   !!
   !! - PengRobinson76
   !! - PengRobinson78
   !! - SoaveRedlichKwong
   !! - RKPR

   private

   public :: PengRobinson76
   public :: PengRobinson78
   public :: SoaveRedlichKwong
   public :: RKPR
   public :: PSRK

contains

   type(CubicEoS) function PengRobinson76(tc, pc, w, kij, lij) result(model)
      !! PengRobinson76.
      !!
      !! Using the critical constants setup the parameters to use the
      !! PengRobinson Equation of State
      !!
      !! - \[\alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2\]
      !! - \[k = 0.37464 + 1.54226 * \omega - 0.26993 \omega^2 \]
      !! - \[a_c = 0.45723553  R^2 T_c^2 / P_c\]
      !! - \[b = 0.07779607r  R T_c/P_c\]
      !! - \[\delta_1 = 1 + \sqrt{2}\]
      !! - \[\delta_2 = 1 - \sqrt{2}\]
      !!
      !! There is also the optional posibility to include the \(k_{ij}\) and
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
      type(AlphaSoave) :: alpha
      integer :: nc
      integer :: i

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
      !! PengRobinson78.
      !!
      !! Using the critical constants setup the parameters to use the
      !! PengRobinson Equation of State
      !!
      !! - \[\alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2\]
      !! - \[k = 0.37464 + 1.54226 \omega - 0.26992 \omega^2  \text{ where } \omega <=0.491\]
      !! - \[k = 0.37464 + 1.48503 \omega - 0.16442 \omega^2  + 0.016666 \omega^3 \text{ where } \omega > 0.491\]
      !! - \[a_c = 0.45723553  R^2 T_c^2 / P_c\]
      !! - \[b = 0.07779607r  R T_c/P_c\]
      !! - \[\delta_1 = 1 + \sqrt{2}\]
      !! - \[\delta_2 = 1 - \sqrt{2}\]
      !!
      !! There is also the optional posibility to include the \(k_{ij}\) and
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
      type(AlphaSoave) :: alpha
      integer :: nc
      integer :: i

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
      !! SoaveRedlichKwong.
      !!
      !! Using the critical constants setup the parameters to use the
      !! SoaveRedlichKwong Equation of State
      !!
      !! - \[\alpha(T_r) = (1 + k (1 - \sqrt{T_r}))^2\]
      !! - \[k = 0.48 + 1.574 \omega - 0.175 \omega^2 \]
      !! - \[a_c = 0.427480  R^2 * T_c^2/P_c\]
      !! - \[b = 0.086640  R T_c/P_c\]
      !! - \[\delta_1 = 1\]
      !! - \[\delta_2 = 0\]
      !!
      !! There is also the optional posibility to include the k_{ij} and l_{ij}
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
      type(AlphaSoave) :: alpha
      integer :: nc
      integer :: i

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
      use yaeos__models_ar_genericcubic, only: CubicEoS
      use yaeos__models_ar_cubic_alphas, only: AlphaMathiasCopeman, AlphaSoave
      use yaeos__models_cubic_mixing_rules_huron_vidal, only: MHV
      use yaeos__models_ge_implementations, only: setup_psrk, UNIFAC
      use yaeos__models_ge_group_contribution_groups, only: Groups
      real(pr), intent(in) :: tc(:) !! Critical temperature [K]
      real(pr), intent(in) :: pc(:) !! Critical pressure [bar]
      real(pr), intent(in) :: w(:) !! Acentric factor
      type(Groups), intent(in) :: molecules(:)
      real(pr), optional, intent(in) :: c1(:), c2(:), c3(:)

      type(UNIFAC) :: ge
      type(Substances) :: composition
      type(MHV) :: mixrule
      type(AlphaSoave) :: alpha
      type(AlphaMathiasCopeman) :: alpha_mc
      integer :: nc
      integer :: i

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
      !! RKPR Equation of State
      !!
      !! The RKPR EoS extends the classical formulation of Cubic Equations
      !! of State by freeing the parameter \(\delta_1\). This extra degree
      !! provides extra ways of implementing the equation in comparison
      !! of other Cubic EoS (like PR and SRK) which are limited to definition
      !! of their critical constants.
      !!
      !! Besides that extra parameter, the RKRR includes another \(\alpha\)
      !! function:
      !! \[
      !!  \alpha(T_r) = \left(\frac{3}{2+T_r}\right)^k
      !! \]
      !!
      !! In this implementation we take the simplest form which correlates
      !! the extra parameter to the critical compressibility factor \(Z_c\) and
      !! the \(k\) parameter of the \(\alpha\) function to \(Z_c\) and \(\omega\):
      !!
      !! \[\delta_1 = d_1 + d_2 (d_3 - Z_c)^d_4 + d_5 (d_3 - Z_c) ^ d_6\]
      !! \[k = (A_1  Z_c + A_0)\omega^2 + (B_1 Z_c + B_0)\omega + (C_1 Z_c + C_0)\]
      use yaeos__models_ar_cubic_quadratic_mixing, only: QMR_RKPR
      use yaeos__models_ar_cubic_alphas, only: AlphaRKPR
      real(pr), intent(in) :: tc(:) !! Critical Temperature [K]
      real(pr), intent(in) :: pc(:) !! Critical Pressure [bar]
      real(pr), intent(in) :: w(:) !! Acentric Factor
      real(pr), intent(in) :: zc(:) !! Critical compressibility
      real(pr), optional, intent(in) :: kij(:, :) !! k_{ij} matrix
      real(pr), optional, intent(in) :: lij(:, :) !! l_{ij} matrix
      real(pr), optional, intent(in) :: delta_1(:)
      real(pr), optional, intent(in) :: k(:)

      type(AlphaRKPR) :: alpha
      type(QMR_RKPR) :: mixrule
      type(Substances) :: composition

      integer :: i, nc

      real(pr), parameter :: d1 = 0.428364, &
         d2 = 18.496215, &
         d3=0.338426, &
         d4=0.66, &
         d5 = 789.723105, &
         d6=2.512392

      real(pr), parameter :: A1 = -2.4407
      real(pr), parameter :: A0 = 0.0017
      real(pr), parameter :: B1 =7.4513
      real(pr), parameter :: B0 =1.9681
      real(pr), parameter :: C1 =12.504
      real(pr), parameter :: C0 =-2.6238

      real(pr) :: OMa(size(pc)), OMb(size(pc))
      real(pr) :: Zc_eos(size(pc))
      real(pr) :: Psat_i, diff

      nc = size(tc)

      composition%pc = pc
      composition%tc = tc
      composition%w = w

      Zc_eos = 1.168 * Zc

      if (present(k)) then
         alpha%k = k
      else
         alpha%k = (A1 * zc + A0)*w**2 + (B1*zc + B0)*w + (C1*Zc + C0)
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
         model%del1 = d1 + d2 * (d3 - zc) ** d4 + d5 * (d3 - zc) ** d6
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
            diff = 1
            do while (abs(diff) > 1e-6)
               Psat_i = model%Psat_pure(i, 0.7*Tc(i))
               diff = (w(i) - (-1 - log10(Psat_i/Pc(i))))
               alpha%k(i) = alpha%k(i) + 0.1*diff
               model%alpha = alpha
            end do
         end do
      end if
   end function RKPR

   subroutine get_OMa_OMb(del1, OMa, OMb)
      real(pr), intent(in) :: del1(:)
      real(pr), intent(out) :: OMa(size(del1))
      real(pr), intent(out) :: OMb(size(del1))

      real(pr) :: d1(size(del1)), y(size(del1))

      d1 = (1._pr + del1**2._pr)/(1._pr + del1)
      y = 1._pr + (2._pr*(1._pr + del1))**(1.0_pr/3._pr) + (4._pr/(1._pr + del1))**(1.0_pr/3)
      OMa = (3._pr*y*y + 3._pr*y*d1 + d1**2._pr + d1 - 1.0_pr)/(3._pr*y + d1 - 1.0_pr)**2._pr
      OMb = 1._pr/(3._pr*y + d1 - 1.0_pr)
   end subroutine get_OMa_OMb
end module yaeos__models_ar_cubic_implementations
