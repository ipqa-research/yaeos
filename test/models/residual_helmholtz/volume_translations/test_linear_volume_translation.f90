program test_linear_volume_translation
   use yaeos, only: pr, R, ArModel, PengRobinson76, LinearVolumeTranslation, ConstantVolumeTranslation
   use yaeos__consistency, only: numeric_ar_derivatives, ar_consistency, individual_ar_calls
   use auxiliar_functions, only: allclose, rel_error
   use testing_aux, only: assert, test_title
   implicit none

   write(*, *) test_title("Linear volume translation")
   call test_translation()
   call test_helmholtz_consistency()
   call test_translated_properties()
   call test_property_derivatives()
   call test_constant_limit()
   call test_excess_properties()
   write(*, *) ""

contains

   subroutine test_translation()
      class(ArModel), allocatable :: base
      type(LinearVolumeTranslation) :: linear
      real(pr) :: n(4)
      real(pr) :: a(4)
      real(pr) :: b(4)
      real(pr) :: T
      real(pr) :: C
      real(pr) :: CT
      real(pr) :: CT2
      real(pr) :: Cn(4)
      real(pr) :: CTn(4)
      real(pr) :: Cn2(4,4)

      base = PengRobinson76(&
         [190.564_pr, 305.32_pr, 369.83_pr, 425.12_pr], &
         [45.99_pr, 48.72_pr, 42.48_pr, 37.96_pr], &
         [0.0115_pr, 0.0995_pr, 0.1523_pr, 0.2002_pr] &
         )
      n = [1.5_pr, 0.2_pr, 0.7_pr, 2.3_pr]
      a = [0.010_pr, 0.020_pr, -0.005_pr, 0.030_pr]
      b = [0.00001_pr, -0.00002_pr, 0.00003_pr, -0.00001_pr]
      T = 600._pr
      linear = LinearVolumeTranslation(base, a, b)

      ! The total shift is extensive; optional outputs also work individually.
      call linear%translation(n, T, C, CT, CT2, Cn, CTn, Cn2)

      call assert(abs(C-dot_product(n,a+b*T)) < 1.e-14_pr, "C is the total shift")
      call assert(abs(CT-dot_product(n,b)) < 1.e-14_pr, "CT")
      call assert(CT2 == 0._pr .and. all(Cn2 == 0._pr), "Second derivatives vanish")
      call assert(allclose(Cn,a+b*T,1.e-12_pr), "Cn")
      call assert(allclose(CTn,b,1.e-12_pr), "CTn")

      call linear%translation(2*n, T, C)
      call assert(abs(C-2*dot_product(n,a+b*T)) < 1.e-14_pr, "Extensive C")

      call linear%translation(n, T, C, CT=CT)
      call assert(abs(CT-dot_product(n,b)) < 1.e-14_pr, "Individual CT")

      call linear%translation(n, T, C, CT2=CT2)
      call assert(CT2 == 0._pr, "Individual CT2")

      call linear%translation(n, T, C, Cn=Cn)
      call assert(allclose(Cn,a+b*T,1.e-12_pr), "Individual Cn")

      call linear%translation(n, T, C, CTn=CTn)
      call assert(allclose(CTn,b,1.e-12_pr), "Individual CTn")

      call linear%translation(n, T, C, Cn2=Cn2)
      call assert(all(Cn2 == 0._pr), "Individual Cn2")
   end subroutine test_translation

   subroutine test_helmholtz_consistency()
      class(ArModel), allocatable :: base
      type(LinearVolumeTranslation) :: linear
      real(pr) :: n(4)
      real(pr) :: a(4)
      real(pr) :: b(4)
      real(pr) :: T
      real(pr) :: V
      class(ArModel), allocatable :: model
      real(pr) :: Ar
      real(pr) :: ArV
      real(pr) :: ArT
      real(pr) :: ArTV
      real(pr) :: ArV2
      real(pr) :: ArT2
      real(pr) :: Arn(4)
      real(pr) :: ArVn(4)
      real(pr) :: ArTn(4)
      real(pr) :: Arn2(4,4)
      real(pr) :: ArTV_combined
      real(pr) :: ArVn_combined(4)
      real(pr) :: Ar_num
      real(pr) :: ArV_num
      real(pr) :: ArT_num
      real(pr) :: ArTV_num
      real(pr) :: ArV2_num
      real(pr) :: ArT2_num
      real(pr) :: Arn_num(4)
      real(pr) :: ArVn_num(4)
      real(pr) :: ArTn_num(4)
      real(pr) :: Arn2_num(4,4)
      real(pr) :: eq31
      real(pr) :: eq33(4,4)
      real(pr) :: eq34(4)
      real(pr) :: eq36
      real(pr) :: eq37
      real(pr) :: step
      integer :: i
      logical :: passed

      base = PengRobinson76(&
         [190.564_pr, 305.32_pr, 369.83_pr, 425.12_pr], &
         [45.99_pr, 48.72_pr, 42.48_pr, 37.96_pr], &
         [0.0115_pr, 0.0995_pr, 0.1523_pr, 0.2002_pr] &
         )
      n = [1.5_pr, 0.2_pr, 0.7_pr, 2.3_pr]
      a = [0.010_pr, 0.020_pr, -0.005_pr, 0.030_pr]
      b = [0.00001_pr, -0.00002_pr, 0.00003_pr, -0.00001_pr]
      T = 600._pr
      V = 0.8_pr
      linear = LinearVolumeTranslation(base, a, b)

      model = linear

      call individual_ar_calls(model, n, V, T, passed)
      call assert(passed, "Individual Helmholtz outputs")

      call model%residual_helmholtz(n,V,T,Ar,ArV,ArT,ArTV,ArV2,ArT2,Arn,ArVn,ArTn,Arn2)
      call model%residual_helmholtz(n,V,T,ArTV=ArTV_combined,ArVn=ArVn_combined)

      call assert(abs(ArTV-ArTV_combined)<1.e-12_pr, "Combined ArTV")
      call assert(allclose(ArVn,ArVn_combined,1.e-12_pr), "Combined ArVn")

      do i = 1, 2
         step = 1._pr / real(i, pr)
         call numeric_ar_derivatives(&
            model, n, V, T, d_n=step*8.e-4_pr, d_v=step*1.e-4_pr, d_t=step*0.1_pr, &
            Ar=Ar_num, ArV=ArV_num, ArT=ArT_num, ArTV=ArTV_num, ArV2=ArV2_num, &
            ArT2=ArT2_num, Arn=Arn_num, ArVn=ArVn_num, ArTn=ArTn_num, Arn2=Arn2_num &
            )

         call assert(rel_error(Ar, Ar_num) < 1.e-12_pr, "Ar numerical")
         call assert(rel_error(ArV, ArV_num) < 5.e-8_pr, "ArV numerical")
         call assert(rel_error(ArT, ArT_num) < 3.e-8_pr, "ArT numerical")
         call assert(rel_error(ArTV, ArTV_num) < 3.e-7_pr, "ArTV numerical")
         call assert(rel_error(ArV2, ArV2_num) < 5.e-6_pr, "ArV2 numerical")
         call assert(rel_error(ArT2, ArT2_num) < 2.e-6_pr, "ArT2 numerical")
         call assert(allclose(Arn, Arn_num, 5.e-8_pr), "Arn numerical")
         call assert(allclose(ArVn, ArVn_num, 4.e-6_pr), "ArVn numerical")
         call assert(allclose(ArTn, ArTn_num, 2.e-7_pr), "ArTn numerical")
         call assert(maxval(rel_error(Arn2, Arn2_num)) < 5.e-6_pr, "Arn2 numerical")
      end do

      call test_second_derivatives(model, n, V, T, ArV2, ArT2, ArTV, ArVn, ArTn, Arn2)

      call ar_consistency(model,n,V,T,eq31,eq33,eq34,eq36,eq37)

      call assert(abs(eq31)<1.e-12_pr .and. maxval(abs(eq33))<1.e-12_pr, "Gibbs identities")
      call assert(maxval(abs(eq34))<1.e-12_pr, "Gibbs-Duhem")
      call assert(abs(eq36)<1.e-12_pr .and. abs(eq37)<1.e-12_pr, "Pressure and temperature identities")
   end subroutine test_helmholtz_consistency

   subroutine test_second_derivatives(model, n, V, T, ArV2, ArT2, ArTV, ArVn, ArTn, Arn2)
      ! Differencing first derivatives avoids subtracting nearly equal energies
      ! twice. The separate energy-based checks above validate those first derivatives.
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: V
      real(pr), intent(in) :: T
      real(pr), intent(in) :: ArV2
      real(pr), intent(in) :: ArT2
      real(pr), intent(in) :: ArTV
      real(pr), intent(in) :: ArVn(:)
      real(pr), intent(in) :: ArTn(:)
      real(pr), intent(in) :: Arn2(:,:)
      real(pr) :: ArV_plus
      real(pr) :: ArV_minus
      real(pr) :: ArT_plus
      real(pr) :: ArT_minus
      real(pr) :: Arn_plus(size(n))
      real(pr) :: Arn_minus(size(n))
      real(pr) :: n_work(size(n))
      real(pr) :: step
      integer :: i
      integer :: k

      do k = 1, 2
         step = 1.e-5_pr / real(k, pr)
         call model%residual_helmholtz(n, V+step, T, ArV=ArV_plus, ArT=ArT_plus, Arn=Arn_plus)
         call model%residual_helmholtz(n, V-step, T, ArV=ArV_minus, ArT=ArT_minus, Arn=Arn_minus)
         call assert(rel_error(ArV2, (ArV_plus-ArV_minus)/(2*step)) < 1.e-7_pr, "ArV2 from ArV")
         call assert(rel_error(ArTV, (ArT_plus-ArT_minus)/(2*step)) < 1.e-7_pr, "ArTV from ArT")
         call assert(allclose(ArVn, (Arn_plus-Arn_minus)/(2*step), 1.e-7_pr), "ArVn from Arn")

         step = 0.01_pr / real(k, pr)
         call model%residual_helmholtz(n, V, T+step, ArV=ArV_plus, ArT=ArT_plus, Arn=Arn_plus)
         call model%residual_helmholtz(n, V, T-step, ArV=ArV_minus, ArT=ArT_minus, Arn=Arn_minus)
         call assert(rel_error(ArT2, (ArT_plus-ArT_minus)/(2*step)) < 1.e-8_pr, "ArT2 from ArT")
         call assert(rel_error(ArTV, (ArV_plus-ArV_minus)/(2*step)) < 1.e-8_pr, "ArTV from ArV")
         call assert(allclose(ArTn, (Arn_plus-Arn_minus)/(2*step), 1.e-8_pr), "ArTn from Arn")

         step = 1.e-4_pr / real(k, pr)
         do i = 1, size(n)
            n_work = n
            n_work(i) = n(i) + step
            call model%residual_helmholtz(n_work, V, T, Arn=Arn_plus)
            n_work(i) = n(i) - step
            call model%residual_helmholtz(n_work, V, T, Arn=Arn_minus)
            call assert(allclose(Arn2(:,i), (Arn_plus-Arn_minus)/(2*step), 1.e-7_pr), "Arn2 from Arn")
         end do
      end do
   end subroutine test_second_derivatives

   subroutine test_translated_properties()
      class(ArModel), allocatable :: base
      type(LinearVolumeTranslation) :: linear
      real(pr) :: n(4)
      real(pr) :: a(4)
      real(pr) :: b(4)
      real(pr) :: T
      real(pr) :: V
      real(pr) :: W
      real(pr) :: P
      real(pr) :: Pbase
      real(pr) :: C
      real(pr) :: CT
      real(pr) :: Ar
      real(pr) :: Arbase
      real(pr) :: Gr
      real(pr) :: Grbase
      real(pr) :: Sr
      real(pr) :: Srbase
      real(pr) :: Hr
      real(pr) :: Hrbase
      real(pr) :: Cp
      real(pr) :: Cpbase
      real(pr) :: lnphi(4)
      real(pr) :: lnphibase(4)
      real(pr) :: volume_plus
      real(pr) :: volume_minus
      real(pr) :: enthalpy_plus
      real(pr) :: enthalpy_minus
      real(pr), parameter :: step = 0.01_pr
      integer :: i

      base = PengRobinson76(&
         [190.564_pr, 305.32_pr, 369.83_pr, 425.12_pr], &
         [45.99_pr, 48.72_pr, 42.48_pr, 37.96_pr], &
         [0.0115_pr, 0.0995_pr, 0.1523_pr, 0.2002_pr] &
         )
      n = [1.5_pr, 0.2_pr, 0.7_pr, 2.3_pr]
      a = [0.010_pr, 0.020_pr, -0.005_pr, 0.030_pr]
      b = [0.00004_pr, -0.00002_pr, 0.00003_pr, 0.00002_pr]
      V = 0.8_pr
      linear = LinearVolumeTranslation(base, a, b)

      ! Compare corresponding states at two temperatures, with CT /= 0.
      do i = 1, 2
         T = 400._pr + 200._pr*real(i - 1, pr)
         C = dot_product(n, a + b*T)
         CT = dot_product(n, b)
         W = V + C

         call base%pressure(n, W, T, Pbase)
         call linear%pressure(n, V, T, P)
         call assert(rel_error(P, Pbase) < 1.e-12_pr, "Pressure at corresponding volumes")

         ! The logarithm changes the ideal-gas reference from W to V.
         call base%residual_helmholtz(n, W, T, Ar=Arbase)
         call linear%residual_helmholtz(n, V, T, Ar=Ar)
         call assert(abs(Ar - Arbase - sum(n)*R*T*log(V/W)) < 1.e-10_pr, &
            "Helmholtz includes the ideal-reference correction")

         call base%lnphi_vt(n, W, T, lnPhi=lnphibase)
         call linear%lnphi_vt(n, V, T, lnPhi=lnphi)
         call assert(allclose(lnphi, lnphibase - P*(a + b*T)/(R*T), 1.e-12_pr), &
            "Fugacity coefficient correction")

         ! Residual G and S use ideal references at their respective volumes.
         call base%gibbs_residual_vt(n, W, T, Gr=Grbase)
         call linear%gibbs_residual_vt(n, V, T, Gr=Gr)
         call assert(abs(Gr - (Grbase - P*C + sum(n)*R*T*log(V/W))) < 1.e-10_pr, &
            "Residual Gibbs includes the volume reference correction")

         call base%entropy_residual_vt(n, W, T, Sr=Srbase)
         call linear%entropy_residual_vt(n, V, T, Sr=Sr)
         call assert(abs(Sr - (Srbase + P*CT - sum(n)*R*log(V/W))) < 1.e-10_pr, &
            "Residual entropy includes the temperature and reference corrections")

         call base%enthalpy_residual_vt(n, W, T, Hr=Hrbase)
         call linear%enthalpy_residual_vt(n, V, T, Hr=Hr)
         call assert(abs(Hr - (Hrbase - P*(C - T*CT))) < 1.e-10_pr, &
            "Enthalpy changes by -P*(C - T*CT)")

         call base%Cp_residual_vt(n, W, T, Cpbase)
         call linear%Cp_residual_vt(n, V, T, Cp)
         call assert(abs(Cp - Cpbase) < 1.e-10_pr, "Linear shift preserves Cp at fixed P")
      end do

      ! Cp is also the temperature derivative of enthalpy at fixed pressure.
      call linear%volume(n, P, T + step, volume_plus, "vapor")
      call linear%volume(n, P, T - step, volume_minus, "vapor")
      call assert(volume_plus > 0._pr .and. volume_minus > 0._pr, "Cp test volumes are positive")

      call linear%enthalpy_residual_vt(n, volume_plus, T + step, Hr=enthalpy_plus)
      call linear%enthalpy_residual_vt(n, volume_minus, T - step, Hr=enthalpy_minus)
      call assert(abs(Cp - (enthalpy_plus - enthalpy_minus)/(2*step)) < 1.e-8_pr, &
         "Cp agrees with the fixed-pressure enthalpy derivative")

      call base%volume(n, 1._pr, T, W, "vapor")
      call linear%volume(n, 1._pr, T, V, "vapor")
      call assert(W > 0._pr .and. V > 0._pr, "Positive vapor volumes")
      call assert(abs(V - W + C) < 1.e-12_pr, "Translated vapor volume")
   end subroutine test_translated_properties

   subroutine test_property_derivatives()
      ! Check property derivatives independently, at their specified constraints.
      class(ArModel), allocatable :: base
      type(LinearVolumeTranslation) :: model
      real(pr) :: n(2)
      real(pr) :: T
      real(pr) :: V
      real(pr) :: P
      real(pr) :: derivative
      real(pr) :: plus
      real(pr) :: minus
      real(pr) :: volume_plus
      real(pr) :: volume_minus
      real(pr) :: lnphi_p(2)
      real(pr) :: lnphi_plus(2)
      real(pr) :: lnphi_minus(2)
      real(pr), parameter :: dv = 1.e-5_pr
      real(pr), parameter :: dt = 0.01_pr
      real(pr), parameter :: dp = 0.001_pr

      base = PengRobinson76(&
         [190.564_pr, 305.32_pr], [45.99_pr, 48.72_pr], [0.0115_pr, 0.0995_pr] &
         )
      n = [0.8_pr, 1.5_pr]
      T = 400._pr
      V = 0.8_pr
      model = LinearVolumeTranslation(base, [0.010_pr, 0.020_pr], [0.00004_pr, -0.00002_pr])

      ! Enthalpy's volume derivative holds T and n fixed.
      call model%enthalpy_residual_vt(n, V, T, HrV=derivative)
      call model%enthalpy_residual_vt(n, V + dv, T, Hr=plus)
      call model%enthalpy_residual_vt(n, V - dv, T, Hr=minus)
      call assert(abs(derivative - (plus - minus)/(2*dv)) < 1.e-8_pr*max(1._pr, abs(derivative)), &
         "Enthalpy volume derivative")

      ! Entropy's temperature derivative holds V and n fixed.
      call model%entropy_residual_vt(n, V, T, SrT=derivative)
      call model%entropy_residual_vt(n, V, T + dt, Sr=plus)
      call model%entropy_residual_vt(n, V, T - dt, Sr=minus)
      call assert(abs(derivative - (plus - minus)/(2*dt)) < 1.e-9_pr*max(1._pr, abs(derivative)), &
         "Entropy temperature derivative")

      ! Fugacity pressure derivatives hold T and n fixed, so V must be solved again.
      call model%pressure(n, V, T, P)
      call model%lnphi_vt(n, V, T, dlnPhidP=lnphi_p)
      call model%volume(n, P + dp, T, volume_plus, "vapor")
      call model%volume(n, P - dp, T, volume_minus, "vapor")
      call assert(volume_plus > 0._pr .and. volume_minus > 0._pr, "Fugacity test volumes are positive")

      call model%lnphi_vt(n, volume_plus, T, lnPhi=lnphi_plus)
      call model%lnphi_vt(n, volume_minus, T, lnPhi=lnphi_minus)
      call assert(all(abs(lnphi_p - (lnphi_plus - lnphi_minus)/(2*dp)) < 1.e-9_pr), &
         "Fugacity pressure derivatives")
   end subroutine test_property_derivatives

   subroutine test_constant_limit()
      class(ArModel), allocatable :: base
      type(LinearVolumeTranslation) :: linear
      real(pr) :: n(4)
      real(pr) :: a(4)
      real(pr) :: b(4)
      real(pr) :: T
      real(pr) :: V
      type(ConstantVolumeTranslation) :: constant
      real(pr) :: ArT
      real(pr) :: ArT2
      real(pr) :: ArTn(4)
      real(pr) :: ArT_num
      real(pr) :: ArT2_num
      real(pr) :: ArTn_num(4)

      base = PengRobinson76(&
         [190.564_pr, 305.32_pr, 369.83_pr, 425.12_pr], &
         [45.99_pr, 48.72_pr, 42.48_pr, 37.96_pr], &
         [0.0115_pr, 0.0995_pr, 0.1523_pr, 0.2002_pr] &
         )
      n = [1.5_pr, 0.2_pr, 0.7_pr, 2.3_pr]
      a = [0.010_pr, 0.020_pr, -0.005_pr, 0.030_pr]
      b = [0.00001_pr, -0.00002_pr, 0.00003_pr, -0.00001_pr]
      T = 600._pr
      V = 0.8_pr
      ! Zero slopes recover the constant translation, including temperature derivatives.
      linear = LinearVolumeTranslation(base, a, 0._pr*b)
      constant = ConstantVolumeTranslation(base, a)

      call linear%residual_helmholtz(n, V, T, ArT=ArT, ArT2=ArT2, ArTn=ArTn)
      call constant%residual_helmholtz(n, V, T, ArT=ArT_num, ArT2=ArT2_num, ArTn=ArTn_num)

      call assert(abs(ArT-ArT_num)<1.e-12_pr .and. abs(ArT2-ArT2_num)<1.e-12_pr, "Constant limit")
      call assert(maxval(abs(ArTn-ArTn_num))<1.e-12_pr, "Constant limit amount derivatives")
   end subroutine test_constant_limit

   subroutine test_excess_properties()
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

      ! Store the same properties for the base and translated models.
      type :: excess_results
         real(pr) :: lngamma(3)
         real(pr) :: lngamma_alone(3)
         real(pr) :: lngamma_p(3)
         real(pr) :: lngamma_t(3)
         real(pr) :: lngamma_n(3,3)
         real(pr) :: gibbs
         real(pr) :: gibbs_alone
         real(pr) :: gibbs_p
         real(pr) :: gibbs_t
         real(pr) :: gibbs_n(3)
         real(pr) :: enthalpy
         real(pr) :: volume
         real(pr) :: entropy
         real(pr) :: helmholtz
         real(pr) :: internal_energy
      end type excess_results

      class(ArModel), allocatable :: base
      class(ArModel), allocatable :: eos
      type(LinearVolumeTranslation) :: translated
      type(excess_results) :: result(2)
      real(pr) :: n(3)
      real(pr) :: npure(3)
      real(pr) :: kij(3,3)
      real(pr) :: lij(3,3)
      real(pr) :: p
      real(pr) :: t
      real(pr) :: volume
      real(pr), parameter :: tol = 1.e-10_pr
      integer :: state
      integer :: root_index
      integer :: model_index
      integer :: i
      character(len=6) :: root

      ! A ternary PR76 mixture with nonzero binary interactions.
      kij = reshape(&
         [0._pr, 0.06_pr, 0.10_pr, 0.06_pr, 0._pr, 0.03_pr, &
         0.10_pr, 0.03_pr, 0._pr], [3,3])

      lij = 0._pr

      base = PengRobinson76(&
         [190.564_pr, 305.32_pr, 369.83_pr], &
         [45.99_pr, 48.72_pr, 42.48_pr], &
         [0.0115_pr, 0.0995_pr, 0.1523_pr], kij, lij)

      translated = LinearVolumeTranslation(base, [0.002_pr, 0.004_pr, -0.001_pr], &
         [0.00001_pr, -0.000005_pr, 0.000015_pr])

      ! Compare two states and each root choice at the same P, T and amounts.
      do state = 1, 2
         if (state == 1) then
            n = [0.8_pr, 1.5_pr, 0.7_pr]
            p = 10._pr
            t = 240._pr
         else
            n = [1.2_pr, 0.4_pr, 2.1_pr]
            p = 40._pr
            t = 350._pr
         end if

         do root_index = 1, 3
            select case (root_index)
             case (1)
               root = "liquid"
             case (2)
               root = "vapor"
             case (3)
               root = "stable"
            end select

            ! Result 1 belongs to the base EOS; result 2 to its translation.
            do model_index = 1, 2
               if (model_index == 1) then
                  eos = base
               else
                  eos = translated
               end if

               ! The mixture root and all stable pure reference roots must be physical.
               call eos%volume(n, p, t, volume, trim(root))
               call assert(ieee_is_finite(volume) .and. volume > 0._pr, "Excess mixture volume")

               do i = 1, size(n)
                  npure = 0._pr
                  npure(i) = 1._pr
                  call eos%volume(npure, p, t, volume, "stable")
                  call assert(ieee_is_finite(volume) .and. volume > 0._pr, "Excess pure volume")
               end do

               ! Request ln(gamma) and G^E both with and without derivatives,
               ! since those calls follow different paths in ArModel.
               call eos%ln_activity_coefficient(n, p, t, trim(root), &
                  result(model_index)%lngamma, result(model_index)%lngamma_p, &
                  result(model_index)%lngamma_t, result(model_index)%lngamma_n)

               call eos%ln_activity_coefficient(n, p, t, trim(root), lngamma=result(model_index)%lngamma_alone)

               call eos%gibbs_excess(n, p, t, trim(root), result(model_index)%gibbs, &
                  result(model_index)%gibbs_p, result(model_index)%gibbs_t, result(model_index)%gibbs_n)

               call eos%gibbs_excess(n, p, t, trim(root), Ge=result(model_index)%gibbs_alone)

               ! Remaining extensive excess properties at the same state.
               call eos%enthalpy_excess(n, p, t, trim(root), result(model_index)%enthalpy)
               call eos%volume_excess(n, p, t, trim(root), result(model_index)%volume)
               call eos%entropy_excess(n, p, t, trim(root), result(model_index)%entropy)
               call eos%helmholtz_excess(n, p, t, trim(root), result(model_index)%helmholtz)
               call eos%internal_energy_excess(n, p, t, trim(root), result(model_index)%internal_energy)
            end do

            associate(reference => result(1), shifted => result(2))
               ! Avoid a trivial comparison in which every activity coefficient is one.
               call assert(maxval(abs(reference%lngamma)) > 1.e-6_pr, "Nontrivial activity coefficients")

               ! Compare every component using an absolute tolerance.
               call assert(all(abs(reference%lngamma - shifted%lngamma) < &
                  tol), "ln(gamma)")
               call assert(all(abs(reference%lngamma - shifted%lngamma_alone) < &
                  tol), "ln(gamma) alone")
               call assert(all(abs(reference%lngamma_p - shifted%lngamma_p) < &
                  tol), "dln(gamma)/dP")
               call assert(all(abs(reference%lngamma_t - shifted%lngamma_t) < &
                  tol), "dln(gamma)/dT")
               call assert(all(abs(reference%lngamma_n - shifted%lngamma_n) < &
                  tol), "dln(gamma)/dn")

               ! Gibbs excess and its derivatives.
               call assert(abs(reference%gibbs - shifted%gibbs) < &
                  tol, "G^E")
               call assert(abs(reference%gibbs - shifted%gibbs_alone) < &
                  tol, "G^E alone")
               call assert(abs(reference%gibbs_p - shifted%gibbs_p) < &
                  tol, "dG^E/dP")
               call assert(abs(reference%gibbs_t - shifted%gibbs_t) < &
                  tol, "dG^E/dT")
               call assert(all(abs(reference%gibbs_n - shifted%gibbs_n) < &
                  tol), "dG^E/dn")

               ! Other excess properties must also remain unchanged.
               call assert(abs(reference%enthalpy - shifted%enthalpy) < &
                  tol, "H^E")
               call assert(abs(reference%volume - shifted%volume) < &
                  tol, "V^E")
               call assert(abs(reference%entropy - shifted%entropy) < &
                  tol, "S^E")
               call assert(abs(reference%helmholtz - shifted%helmholtz) < &
                  tol, "A^E")
               call assert(abs(reference%internal_energy - shifted%internal_energy) < &
                  tol, "U^E")
            end associate
         end do
      end do
   end subroutine test_excess_properties
end program test_linear_volume_translation
