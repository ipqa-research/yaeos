program test_constant_volume_translation
   use yaeos, only: pr, R, ArModel, PengRobinson76, ConstantVolumeTranslation
   use yaeos__consistency, only: numeric_ar_derivatives, ar_consistency, individual_ar_calls
   use auxiliar_functions, only: allclose, rel_error
   use testing_aux, only: assert, test_title
   implicit none

   write(*, *) test_title("Constant volume translation")

   call test_translation()
   call test_helmholtz_consistency()
   call test_translated_properties()
   call test_zero_translation()

   write(*, *) ""

contains

   subroutine test_translation()
      ! C is an extensive volume, even when the amounts do not sum to one.
      class(ArModel), allocatable :: base
      type(ConstantVolumeTranslation) :: model
      real(pr) :: n(4)
      real(pr) :: c(4)
      real(pr) :: C_total
      real(pr) :: CT
      real(pr) :: CT2
      real(pr) :: Cn(4)
      real(pr) :: CTn(4)
      real(pr) :: Cn2(4, 4)

      base = PengRobinson76(&
         [190.564_pr, 305.32_pr, 369.83_pr, 425.12_pr], &
         [45.99_pr, 48.72_pr, 42.48_pr, 37.96_pr], &
         [0.0115_pr, 0.0995_pr, 0.1523_pr, 0.2002_pr] &
         )

      c = [0.010_pr, 0.020_pr, -0.005_pr, 0.030_pr] ! L/mol
      n = [1.5_pr, 0.2_pr, 0.7_pr, 2.3_pr]             ! mol
      model = ConstantVolumeTranslation(base, c)

      call model%translation(n, 300._pr, C_total, CT, CT2, Cn, CTn, Cn2)

      call assert(abs(C_total - 0.0845_pr) < 1.e-14_pr, "C = sum(n_i*c_i)")
      call assert(allclose(c, Cn, 1.e-12_pr), "The amount derivatives of C are c_i")
      call assert(CT == 0._pr .and. CT2 == 0._pr, "Temperature derivatives of C are zero")
      call assert(all(CTn == 0._pr), "Temperature/amount derivatives of C are zero")
      call assert(all(Cn2 == 0._pr), "The amount Hessian of C is zero")

      ! Changing temperature does not change C; scaling amounts scales C.
      call model%translation(2._pr*n, 500._pr, C_total)
      call assert(abs(C_total - 2._pr*0.0845_pr) < 1.e-14_pr, "C scales with the total amount")

      ! Each optional translation derivative can also be requested on its own.
      call model%translation(n, 300._pr, C_total, CT=CT)
      call assert(CT == 0._pr, "Individual CT")

      call model%translation(n, 300._pr, C_total, CT2=CT2)
      call assert(CT2 == 0._pr, "Individual CT2")

      call model%translation(n, 300._pr, C_total, Cn=Cn)
      call assert(allclose(c, Cn, 1.e-12_pr), "Individual Cn")

      call model%translation(n, 300._pr, C_total, CTn=CTn)
      call assert(all(CTn == 0._pr), "Individual CTn")

      call model%translation(n, 300._pr, C_total, Cn2=Cn2)
      call assert(all(Cn2 == 0._pr), "Individual Cn2")
   end subroutine test_translation

   subroutine test_helmholtz_consistency()
      ! Four components and non-normalized amounts exercise composition derivatives.
      class(ArModel), allocatable :: base
      class(ArModel), allocatable :: model
      real(pr) :: n(4)
      real(pr) :: c(4)
      real(pr) :: V
      real(pr) :: T

      real(pr) :: Ar
      real(pr) :: ArV
      real(pr) :: ArT
      real(pr) :: ArTV
      real(pr) :: ArV2
      real(pr) :: ArT2
      real(pr) :: Arn(4)
      real(pr) :: ArVn(4)
      real(pr) :: ArTn(4)
      real(pr) :: Arn2(4, 4)
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
      real(pr) :: Arn2_num(4, 4)

      real(pr) :: eq31
      real(pr) :: eq33(4, 4)
      real(pr) :: eq34(4)
      real(pr) :: eq36
      real(pr) :: eq37
      real(pr) :: step
      integer :: i
      logical :: passed

      n = [1.5_pr, 0.2_pr, 0.7_pr, 2.3_pr]
      c = [0.010_pr, 0.020_pr, -0.005_pr, 0.030_pr]
      T = 600._pr
      V = 0.8_pr

      base = PengRobinson76(&
         [190.564_pr, 305.32_pr, 369.83_pr, 425.12_pr], &
         [45.99_pr, 48.72_pr, 42.48_pr, 37.96_pr], &
         [0.0115_pr, 0.0995_pr, 0.1523_pr, 0.2002_pr] &
         )
      model = ConstantVolumeTranslation(base, c)

      ! Each optional Helmholtz output must work when requested on its own.
      call individual_ar_calls(model, n, V, T, passed)
      call assert(passed, "Individual Helmholtz calls agree with the full call")

      ! Compare every analytic derivative with differences of the energy.
      call model%residual_helmholtz(&
         n, V, T, Ar=Ar, ArV=ArV, ArT=ArT, ArTV=ArTV, ArV2=ArV2, &
         ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
         )

      ! Mixed volume derivatives are also requested together by property calls.
      call model%residual_helmholtz(n, V, T, ArTV=ArTV_combined, ArVn=ArVn_combined)
      call assert(rel_error(ArTV, ArTV_combined) < 1.e-12_pr, "Combined ArTV agrees with the full call")
      call assert(allclose(ArVn, ArVn_combined, 1.e-12_pr), "Combined ArVn agrees with the full call")

      ! Composition steps are larger because second differences lose precision.
      ! Two step sizes guard against a fortuitously accurate finite difference.
      do i = 1, 2
         step = 1._pr / real(i, pr)
         call numeric_ar_derivatives(&
            model, n, V, T, d_n=step*8.e-4_pr, d_v=step*1.e-4_pr, d_t=step*0.1_pr, &
            Ar=Ar_num, ArV=ArV_num, ArT=ArT_num, ArTV=ArTV_num, ArV2=ArV2_num, &
            ArT2=ArT2_num, Arn=Arn_num, ArVn=ArVn_num, ArTn=ArTn_num, Arn2=Arn2_num &
            )

         call assert(rel_error(Ar, Ar_num) < 1.e-6_pr, "Ar vs numerical value")
         call assert(rel_error(ArV, ArV_num) < 1.e-7_pr, "ArV vs numerical derivative")
         call assert(rel_error(ArT, ArT_num) < 1.e-7_pr, "ArT vs numerical derivative")
         call assert(rel_error(ArTV, ArTV_num) < 1.e-6_pr, "ArTV vs numerical derivative")
         call assert(rel_error(ArV2, ArV2_num) < 2.e-5_pr, "ArV2 vs numerical derivative")
         call assert(rel_error(ArT2, ArT2_num) < 1.e-5_pr, "ArT2 vs numerical derivative")
         call assert(allclose(Arn, Arn_num, 1.e-7_pr), "Arn vs numerical derivatives")
         call assert(allclose(ArVn, ArVn_num, 3.e-5_pr), "ArVn vs numerical derivatives")
         call assert(allclose(ArTn, ArTn_num, 1.e-6_pr), "ArTn vs numerical derivatives")
         call assert(maxval(rel_error(Arn2, Arn2_num)) < 1.e-5_pr, "Arn2 vs numerical derivatives")
      end do

      ! Michelsen and Mollerup's thermodynamic consistency identities.
      call ar_consistency(model, n, V, T, eq31, eq33, eq34, eq36, eq37)

      call assert(abs(eq31) < 1.e-12_pr, "Consistency: Gibbs energy and fugacities")
      call assert(maxval(abs(eq33)) < 1.e-12_pr, "Consistency: symmetric fugacity derivatives")
      call assert(maxval(abs(eq34)) < 1.e-12_pr, "Consistency: Gibbs-Duhem relation")
      call assert(abs(eq36) < 1.e-12_pr, "Consistency: pressure derivative")
      call assert(abs(eq37) < 1.e-12_pr, "Consistency: temperature derivative and enthalpy")
   end subroutine test_helmholtz_consistency

   subroutine test_translated_properties()
      ! Compare corresponding states of a binary EOS with and without a shift.
      class(ArModel), allocatable :: base
      type(ConstantVolumeTranslation) :: model
      real(pr) :: n(2)
      real(pr) :: c(2)
      real(pr) :: C_total
      real(pr) :: V
      real(pr) :: W
      real(pr) :: T
      real(pr) :: P
      real(pr) :: P_base
      real(pr) :: Ar
      real(pr) :: Ar_base
      real(pr) :: lnphi(2)
      real(pr) :: lnphi_base(2)
      real(pr) :: Hr
      real(pr) :: Hr_base
      real(pr) :: Gr
      real(pr) :: Gr_base
      real(pr) :: Sr
      real(pr) :: Sr_base
      real(pr) :: Cp
      real(pr) :: Cp_base
      character(len=6) :: root
      integer :: i

      base = PengRobinson76(&
         [190.564_pr, 305.32_pr], [45.99_pr, 48.72_pr], [0.01142_pr, 0.0995_pr] &
         )
      n = [0.8_pr, 1.5_pr]
      c = [0.012_pr, 0.021_pr]
      model = ConstantVolumeTranslation(base, c)

      C_total = dot_product(n, c)
      V = 0.8_pr
      W = V + C_total
      T = 330._pr

      ! The shifted EOS at V has the same pressure as the base EOS at V+C.
      call base%pressure(n, W, T, P_base)
      call model%pressure(n, V, T, P)
      call assert(rel_error(P_base, P) < 1.e-12_pr, "Translated pressure equals base pressure at V+C")

      call base%residual_helmholtz(n, W, T, Ar=Ar_base)
      call model%residual_helmholtz(n, V, T, Ar=Ar)
      call assert(abs(Ar - Ar_base - sum(n)*R*T*log(V/W)) < 1.e-12_pr, &
         "Helmholtz includes the ideal-reference correction")

      ! For a linear shift, ln(phi_i) changes by -P*c_i/(RT).
      call base%lnphi_vt(n, W, T, lnPhi=lnphi_base)
      call model%lnphi_vt(n, V, T, lnPhi=lnphi)
      call assert(allclose(lnphi_base-P*c/(R*T), lnphi, 1.e-12_pr), "Fugacity coefficient correction")

      ! These residual properties use the ideal reference at each state's volume.
      call base%gibbs_residual_vt(n, W, T, Gr=Gr_base)
      call model%gibbs_residual_vt(n, V, T, Gr=Gr)
      call assert(abs(Gr - (Gr_base - P*C_total + sum(n)*R*T*log(V/W))) < 1.e-10_pr, &
         "Residual Gibbs includes the volume reference correction")

      call base%entropy_residual_vt(n, W, T, Sr=Sr_base)
      call model%entropy_residual_vt(n, V, T, Sr=Sr)
      call assert(abs(Sr - (Sr_base - sum(n)*R*log(V/W))) < 1.e-10_pr, &
         "Residual entropy includes the volume reference correction")

      call base%enthalpy_residual_vt(n, W, T, Hr=Hr_base)
      call model%enthalpy_residual_vt(n, V, T, Hr=Hr)
      call assert(rel_error(Hr_base-P*C_total, Hr) < 1.e-12_pr, "Residual enthalpy changes by -P*C")

      call base%Cp_residual_vt(n, W, T, Cp_base)
      call model%Cp_residual_vt(n, V, T, Cp)
      call assert(rel_error(Cp_base, Cp) < 1.e-12_pr, "Constant shift preserves residual Cp")

      ! The inherited volume method uses the base solver for each root type.
      do i = 1, 3
         select case (i)
          case (1)
            root = "liquid"
          case (2)
            root = "vapor"
          case (3)
            root = "stable"
         end select

         call base%volume(n, 1._pr, 200._pr, W, trim(root))
         call model%volume(n, 1._pr, 200._pr, V, trim(root))
         call assert(W > 0._pr .and. V > 0._pr, trim(root) // " volumes are positive")
         call assert(rel_error(W-C_total, V) < 1.e-12_pr, trim(root) // " volume is shifted by C")
         call model%pressure(n, V, 200._pr, P)
         call assert(abs(P-1._pr) < 1.e-6_pr, trim(root) // " volume reproduces the specified pressure")
      end do

      W = base%get_v0(n, 1._pr, 200._pr)
      V = model%get_v0(n, 1._pr, 200._pr)
      call assert(rel_error(W-C_total, V) < 1.e-12_pr, "Initial volume is shifted by C")
   end subroutine test_translated_properties

   subroutine test_zero_translation()
      ! With c_i=0 the translated model must reduce to the original EOS.
      class(ArModel), allocatable :: base
      type(ConstantVolumeTranslation) :: model
      real(pr) :: n(2)
      real(pr) :: Ar
      real(pr) :: Ar_base
      real(pr) :: Arn(2)
      real(pr) :: Arn_base(2)
      real(pr) :: V
      real(pr) :: V_base

      base = PengRobinson76(&
         [190.564_pr, 305.32_pr], [45.99_pr, 48.72_pr], [0.01142_pr, 0.0995_pr] &
         )
      model = ConstantVolumeTranslation(base, [0._pr, 0._pr])
      n = [0.8_pr, 1.5_pr]

      call base%residual_helmholtz(n, 0.8_pr, 330._pr, Ar=Ar_base, Arn=Arn_base)
      call model%residual_helmholtz(n, 0.8_pr, 330._pr, Ar=Ar, Arn=Arn)
      call assert(rel_error(Ar_base, Ar) < 1.e-12_pr, "Zero translation preserves Helmholtz energy")
      call assert(allclose(Arn_base, Arn, 1.e-12_pr), "Zero translation preserves amount derivatives")

      call base%volume(n, 1._pr, 200._pr, V_base, "vapor")
      call model%volume(n, 1._pr, 200._pr, V, "vapor")
      call assert(V_base > 0._pr .and. V > 0._pr, "Zero-shift vapor volumes are positive")
      call assert(rel_error(V_base, V) < 1.e-12_pr, "Zero translation preserves volume")
   end subroutine test_zero_translation

end program test_constant_volume_translation
