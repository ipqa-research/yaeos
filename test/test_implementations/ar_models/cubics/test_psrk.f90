module test_psrk
   use yaeos, only: pr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use auxiliar_functions, only: allclose, rel_error
   implicit none

contains
   subroutine collect_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("Test PSRK Ar consistency mix", test_psrk_cons_mixture), &
         new_unittest("Test PSRK Ar consistency pure", test_psrk_cons_pure), &
         new_unittest("Test PSRK paper 2005 1", test_psrk_paper2005_1), &
         new_unittest("Test PSRK paper 2005 2", test_psrk_paper2005_2), &
         new_unittest("Test PSRK paper 2005 3", test_psrk_paper2005_3), &
         new_unittest("Test PSRK paper 2005 4", test_psrk_paper2005_4) &
         ]
   end subroutine collect_suite

   subroutine test_psrk_cons_mixture(error)
      use yaeos, only: pr, PSRK, ArModel, Groups
      use yaeos__consistency, only: numeric_ar_derivatives, ar_consistency
      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: model
      real(pr) :: tc(2), pc(2), w(2), C(2, 3)
      type(Groups) :: molecules(2)

      real(pr) :: n(2), t, v

      real(pr) :: Ar, ArV, ArT, Arn(size(n)), ArV2, ArT2, ArTV
      real(pr) :: ArVn(size(n)), ArTn(size(n)), Arn2(size(n), size(n))

      real(pr) :: Ar_num, ArV_num, ArT_num, Arn_num(size(n)), ArV2_num, ArT2_num
      real(pr) :: ArTV_num, ArVn_num(size(n)), ArTn_num(size(n))
      real(pr) :: Arn2_num(size(n), size(n))
      real(pr) :: eq31, eq33(size(n), size(n)), eq34(size(n)), eq36, eq37

      n = [60, 40]
      tc = [304.21_pr, 553.8_pr]
      pc = [7.383e6_pr, 4.080358e6_pr] / 1e5
      w = [0.223621_pr, 0.213_pr]
      c(1, :) = [0.8255_pr, 0.16755_pr, -1.7039_pr]
      c(2, :)= [0.84082_pr, -0.39847_pr, 0.94148_pr]

      molecules(1)%groups_ids = [117]
      molecules(1)%number_of_groups = [2]
      molecules(2)%groups_ids = [2]
      molecules(2)%number_of_groups = [7]

      t = 600_pr
      v = 10.0_pr

      ! ========================================================================
      ! Model without kij and lij
      ! ------------------------------------------------------------------------
      model = PSRK(tc, pc, w, molecules, c1=C(:, 1), c2=C(:, 2), c3=C(:, 3))

      call model%residual_helmholtz(&
         n, v, t, Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn, &
         ArTV=ArTV, ArV2=ArV2, ArT2=ArT2, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2)

      call numeric_ar_derivatives(&
         model, n, v, t, d_n = 0.01_pr, d_v = 0.0001_pr, d_t = 0.01_pr, &
         Ar=Ar_num, ArV=ArV_num, ArT=ArT_num, ArTV=ArTV_num, ArV2=ArV2_num, &
         ArT2=ArT2_num, Arn=Arn_num, ArVn=ArVn_num, ArTn=ArTn_num, &
         Arn2=Arn2_num &
         )

      call ar_consistency(&
         model, n, v, t, eq31=eq31, eq33=eq33, eq34=eq34, eq36=eq36, eq37=eq37 &
         )

      ! Numeric derivatives
      call check(error, rel_error(Ar, Ar_num) < 1e-10)
      call check(error, rel_error(ArV, ArV_num) < 1e-5)
      call check(error, rel_error(ArT, ArT_num) < 1e-5)
      call check(error, allclose(Arn, Arn_num, 1e-5_pr))
      call check(error, rel_error(ArV2, ArV2_num) < 1e-5)
      call check(error, rel_error(ArT2, ArT2_num) < 1e-6)
      call check(error, rel_error(ArTV, ArTV_num) < 1e-6)
      call check(error, allclose(ArVn, ArVn_num, 1e-6_pr))
      call check(error, allclose(ArTn, ArTn_num, 1e-6_pr))
      call check(error, maxval(rel_error(Arn2, Arn2_num)) < 1e-5)

      ! Consistency tests
      call check(error, abs(eq31) <= 1e-12)
      call check(error, maxval(abs(eq33)) < 1e-13)
      call check(error, maxval(abs(eq34)) < 1e-13)
      call check(error, abs(eq36) <= 1e-13)
      call check(error, abs(eq37) <= 1e-13)
   end subroutine test_psrk_cons_mixture

   subroutine test_psrk_cons_pure(error)
      use yaeos, only: pr, PSRK, ArModel, Groups
      use yaeos__consistency, only: numeric_ar_derivatives, ar_consistency
      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: model
      type(Groups) :: molecules(1)
      real(pr) :: tc(1), pc(1), w(1), c(1, 3)

      real(pr) :: n(1), t, v

      real(pr) :: Ar, ArV, ArT, Arn(size(n)), ArV2, ArT2, ArTV
      real(pr) :: ArVn(size(n)), ArTn(size(n)), Arn2(size(n), size(n))

      real(pr) :: Ar_num, ArV_num, ArT_num, Arn_num(size(n)), ArV2_num, ArT2_num
      real(pr) :: ArTV_num, ArVn_num(size(n)), ArTn_num(size(n))
      real(pr) :: Arn2_num(size(n), size(n))

      real(pr) :: eq31, eq33(size(n), size(n)), eq34(size(n)), eq36, eq37

      n = [5]
      tc = [304.21_pr]
      pc = [7.383e6_pr] / 1e5
      w = [0.223621_pr]
      c(1, :) = [0.8255_pr, 0.16755_pr, -1.7039_pr]

      molecules(1)%groups_ids = [117]
      molecules(1)%number_of_groups = [2]

      t = 600_pr
      v = 15._pr

      model = PSRK(tc, pc, w, molecules, c1=c(:, 1), c2=c(:, 2), c3=c(:, 3))

      call model%residual_helmholtz(&
         n, v, t, Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn, &
         ArTV=ArTV, ArV2=ArV2, ArT2=ArT2, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2)

      call numeric_ar_derivatives(&
         model, n, v, t, d_n = 0.001_pr, d_v = 0.001_pr, d_t = 0.001_pr, &
         Ar=Ar_num, ArV=ArV_num, ArT=ArT_num, ArTV=ArTV_num, ArV2=ArV2_num, &
         ArT2=ArT2_num, Arn=Arn_num, ArVn=ArVn_num, ArTn=ArTn_num, &
         Arn2=Arn2_num &
         )

      call ar_consistency(&
         model, n, v, t, eq31=eq31, eq33=eq33, eq34=eq34, eq36=eq36, eq37=eq37 &
         )

      ! Numeric derivatives
      call check(error, rel_error(Ar, Ar_num) < 1e-10)
      call check(error, rel_error(ArV, ArV_num) < 1e-6)
      call check(error, rel_error(ArT, ArT_num) < 1e-6)
      call check(error, allclose(Arn, Arn_num, 1e-6_pr))
      call check(error, rel_error(ArV2, ArV2_num) < 1e-4)
      call check(error, rel_error(ArT2, ArT2_num) < 1e-4)
      call check(error, rel_error(ArTV, ArTV_num) < 1e-6)
      call check(error, allclose(ArVn, ArVn_num, 1e-5_pr))
      call check(error, allclose(ArTn, ArTn_num, 1e-5_pr))
      call check(error, maxval(rel_error(Arn2, Arn2_num)) < 1e-5)

      ! Consistency tests
      call check(error, abs(eq31) <= 1e-13)
      call check(error, maxval(abs(eq33)) < 1e-13)
      call check(error, maxval(abs(eq34)) < 1e-13)
      call check(error, abs(eq36) <= 1e-13)
      call check(error, abs(eq37) <= 1e-13)
   end subroutine test_psrk_cons_pure

   subroutine test_psrk_paper2005_1(error)
      use yaeos, only: pr, PSRK, ArModel, Groups, saturation_pressure
      use yaeos, only: EquilibriumState
      type(error_type), allocatable, intent(out) :: error
      real(pr) :: x1_exp(14), y1_exp(16), pbub_exp(14), pdew_exp(16), t
      integer :: i, j
      type(EquilibriumState) :: bubble, dew

      class(ArModel), allocatable :: model
      real(pr) :: tc(2), pc(2), w(2), C(2,3)
      type(Groups) :: molecules(2)

      tc = [456.8_pr, 647.3_pr]
      pc = [53.9049_pr, 220.48321_pr]
      w = [0.407_pr, 0.344_pr]
      C(1,:) = [1.62716_pr, -5.00731_pr, 10.34189_pr]
      C(2,:) = [1.0783_pr, -0.58321_pr, 0.54619_pr]

      x1_exp = [0.00754_pr, 0.04774_pr, 0.07538_pr, 0.12312_pr, &
         0.17839_pr, 0.26382_pr, 0.33668_pr, 0.4196_pr, 0.5402_pr, &
         0.63317_pr, 0.74623_pr, 0.86935_pr, 0.94221_pr, 0.97739_pr]

      pbub_exp = [0.06384039999999999_pr, 0.2374065_pr, &
         0.32718200000000003_pr, 0.4349127_pr, 0.5067332_pr, &
         0.5645884999999999_pr, 0.5865337_pr, 0.5985037_pr, 0.6104738_pr, &
         0.6244389_pr, 0.6503741_pr, 0.6962594_pr, 0.7341646000000001_pr, &
         0.7541146999999999_pr]

      y1_exp = [0.01759_pr, 0.19849_pr, 0.3593_pr, 0.5_pr, &
         0.57789_pr, 0.68844_pr, 0.76382_pr, 0.83417_pr, 0.88442_pr, &
         0.91457_pr, 0.94472_pr, 0.95729_pr, 0.96482_pr, 0.97236_pr, &
         0.9799_pr, 0.98995_pr]

      pdew_exp = [0.021945100000000002_pr, 0.035910199999999996_pr, &
         0.0478803_pr, 0.0578554_pr, 0.06384039999999999_pr, &
         0.07182039999999999_pr, 0.09177060000000001_pr, &
         0.12169579999999999_pr, 0.1715711_pr, 0.2354115_pr, &
         0.3431421_pr, 0.42892769999999997_pr, 0.5087282_pr, &
         0.6024938_pr, 0.6663342_pr, 0.7261845_pr]

      molecules(1)%groups_ids = [145]
      molecules(1)%number_of_groups = [1]
      molecules(2)%groups_ids = [16]
      molecules(2)%number_of_groups = [1]

      model = PSRK(tc, pc, w, molecules, c1=C(:, 1), c2=C(:, 2), c3=C(:, 3))
      t = 291.15_pr

      do i=1,14
         bubble = saturation_pressure(model, [x1_exp(i), 1-x1_exp(i)], t, kind="bubble", p0=pbub_exp(i))
         call check(error, abs(bubble%p-pbub_exp(i)) < 4e-3)
      end do

      do j = 1,16
         dew = saturation_pressure(model, [y1_exp(j), 1-y1_exp(j)], t, kind="dew", p0=pdew_exp(j))
         call check(error, abs(dew%p-pdew_exp(j)) < 2e-2)
      end do
   end subroutine test_psrk_paper2005_1

   subroutine test_psrk_paper2005_2(error)
      use yaeos, only: pr, PSRK, ArModel, Groups, saturation_pressure
      use yaeos, only: EquilibriumState
      type(error_type), allocatable, intent(out) :: error
      real(pr) :: x1_exp(12), y1_exp(16), pbub_exp(12), pdew_exp(16), t
      integer :: i, j
      type(EquilibriumState) :: bubble, dew

      class(ArModel), allocatable :: model
      real(pr) :: tc(2), pc(2), w(2), C(2,3)
      type(Groups) :: molecules(2)

      tc = [324.6_pr, 417.0_pr]
      pc = [83.0865_pr, 77.007_pr]
      w = [0.12_pr, 0.073_pr]
      C(1,:) = [0.66635_pr, 0.35497_pr, -1.3766_pr]
      C(2,:) = [0.55192_pr, 0.01934_pr, 0.59414_pr]

      x1_exp = [0.01942_pr, 0.07864_pr, 0.14284_pr, 0.20951_pr, 0.28609_pr, &
         0.36516_pr, 0.44671_pr, 0.53322_pr, 0.62219_pr, 0.73341_pr, 0.84957_pr, 0.96077_pr]

      pbub_exp = [4.6423_pr, 6.8492_pr, 8.8347_pr, 10.746199999999998_pr, &
         12.5833_pr, 14.2728_pr, 15.814699999999998_pr, 17.43_pr, 19.045_pr, &
         21.027_pr, 23.229799999999997_pr, 25.506600000000002_pr]

      y1_exp = [0.06402_pr, 0.1432_pr, 0.23227_pr, 0.30896_pr, 0.39059_pr, &
         0.45983_pr, 0.54885_pr, 0.6131_pr, 0.67487_pr, 0.72425_pr, 0.77361_pr, &
         0.818_pr, 0.86483_pr, 0.91409_pr, 0.9535_pr, 0.9806_pr]

      pdew_exp = [4.049300000000001_pr, 4.412_pr, 4.9215_pr, 5.4318_pr, &
         6.089099999999999_pr, 6.8210999999999995_pr, 7.9939_pr, 9.1686_pr, &
         10.564499999999999_pr, 12.0351_pr, 13.874200000000002_pr, &
         15.861099999999999_pr, 18.2901_pr, 21.2348_pr, 23.5907_pr, 25.136599999999998_pr]

      molecules(1)%groups_ids = [130]
      molecules(1)%number_of_groups = [1]
      molecules(2)%groups_ids = [143]
      molecules(2)%number_of_groups = [1]

      model = PSRK(tc, pc, w, molecules, c1=C(:, 1), c2=C(:, 2), c3=C(:, 3))
      t = 273.15_pr

      do i=1,12
         bubble = saturation_pressure(model, [x1_exp(i), 1-x1_exp(i)], t, kind="bubble", p0=pbub_exp(i))
         call check(error, rel_error(bubble%p,pbub_exp(i)) < 1.1e-2)
      end do

      do j = 1,16
         dew = saturation_pressure(model, [y1_exp(j), 1-y1_exp(j)], t, kind="dew", p0=pdew_exp(j))
         call check(error, rel_error(dew%p,pdew_exp(j)) < 7e-3)
      end do
   end subroutine test_psrk_paper2005_2

   subroutine test_psrk_paper2005_3(error)
      use yaeos, only: pr, PSRK, ArModel, Groups, saturation_pressure
      use yaeos, only: EquilibriumState
      type(error_type), allocatable, intent(out) :: error
      real(pr) :: x1_exp(13), y1_exp(13), pbub_exp(13), pdew_exp(13), t
      integer :: i, j
      type(EquilibriumState) :: bubble, dew

      class(ArModel), allocatable :: model
      real(pr) :: tc(2), pc(2), w(2), C(2,3)
      type(Groups) :: molecules(2)

      tc = [324.6_pr, 417.0_pr]
      pc = [83.0865_pr, 77.007_pr]
      w = [0.12_pr, 0.073_pr]
      C(1,:) = [0.66635_pr, 0.35497_pr, -1.3766_pr]
      C(2,:) = [0.55192_pr, 0.01934_pr, 0.59414_pr]

      x1_exp = [0.02719_pr, 0.10883_pr, 0.17316_pr, 0.24987_pr, 0.31174_pr, &
         0.36618_pr, 0.42805_pr, 0.50477_pr, 0.58149_pr, 0.67553_pr, 0.81659_pr, &
         0.89826_pr, 0.96508_pr]

      pbub_exp = [0.514_pr, 0.9501999999999999_pr, 1.2403_pr, 1.5295_pr, &
         1.7460999999999998_pr, 1.8895_pr, 2.0322999999999998_pr, 2.1741_pr, &
         2.3896_pr, 2.62_pr, 2.8882000000000003_pr, 3.1033_pr, 3.2458_pr]

      y1_exp = [0.02968_pr, 0.10641_pr, 0.18561_pr, 0.26977_pr, 0.35888_pr, &
         0.41581_pr, 0.52719_pr, 0.61876_pr, 0.69548_pr, 0.79446_pr, 0.86621_pr, &
         0.92062_pr, 0.97006_pr]

      pdew_exp = [0.2927_pr, 0.287_pr, 0.3549_pr, 0.3487_pr, 0.4158_pr, &
         0.4116_pr, 0.4771_pr, 0.6178_pr, 0.7595000000000001_pr, 1.1208_pr, &
         1.5578_pr, 2.1434_pr, 2.8769_pr]

      molecules(1)%groups_ids = [130]
      molecules(1)%number_of_groups = [1]
      molecules(2)%groups_ids = [143]
      molecules(2)%number_of_groups = [1]

      model = PSRK(tc, pc, w, molecules, c1=C(:, 1), c2=C(:, 2), c3=C(:, 3))
      t = 213.15_pr

      do i=1,13
         bubble = saturation_pressure(model, [x1_exp(i), 1-x1_exp(i)], t, kind="bubble", p0=pbub_exp(i))
         call check(error, rel_error(bubble%p,pbub_exp(i)) < 1e-1)
      end do

      do j = 1,13
         dew = saturation_pressure(model, [y1_exp(j), 1-y1_exp(j)], t, kind="dew", p0=pdew_exp(j))
         call check(error, rel_error(dew%p,pdew_exp(j)) < 0.2)
      end do
   end subroutine test_psrk_paper2005_3

   subroutine test_psrk_paper2005_4(error)
      use yaeos, only: pr, PSRK, ArModel, Groups, saturation_pressure
      use yaeos, only: EquilibriumState
      type(error_type), allocatable, intent(out) :: error
      real(pr) :: x1_exp(10), y1_exp(12), pbub_exp(10), pdew_exp(12), t
      integer :: i, j
      type(EquilibriumState) :: bubble, dew

      class(ArModel), allocatable :: model
      real(pr) :: tc(2), pc(2), w(2), C(2,3)
      type(Groups) :: molecules(2)

      tc = [440.0_pr, 523.2_pr]
      pc = [91.1925_pr, 38.300850000000004_pr]
      w = [0.318_pr, 0.363_pr]
      C(1,:) = [0.96273_pr, 0.10351_pr, -1.6102_pr]
      C(2,:) = [1.0408_pr, -0.17686_pr, 0.49506_pr]

      x1_exp = [0.03202_pr, 0.12562_pr, 0.24877_pr, 0.34975_pr, 0.46798_pr, &
         0.58867_pr, 0.68473_pr, 0.80049_pr, 0.94089_pr, 0.9803_pr]

      pbub_exp = [0.014705900000000001_pr, 0.0392157_pr, 0.0784314_pr, &
         0.11764709999999999_pr, 0.17279409999999998_pr, 0.2389706_pr, &
         0.2965686_pr, 0.372549_pr, 0.4607843_pr, 0.48284309999999997_pr]

      y1_exp = [0.03695_pr, 0.13547_pr, 0.26601_pr, 0.45074_pr, 0.61576_pr, &
         0.74877_pr, 0.85714_pr, 0.94335_pr, 0.97291_pr, 0.98768_pr, 0.99261_pr, 1.0_pr]

      pdew_exp = [0.0073529_pr, 0.0085784_pr, 0.009803899999999999_pr, &
         0.013480399999999998_pr, 0.019607799999999998_pr, 0.0281863_pr, &
         0.0465686_pr, 0.09313729999999999_pr, 0.1446078_pr, &
         0.21446079999999998_pr, 0.2830882_pr, 0.47058819999999996_pr]

      molecules(1)%groups_ids = [149]
      molecules(1)%number_of_groups = [1]
      molecules(2)%groups_ids = [1, 2, 21]
      molecules(2)%number_of_groups = [1, 1, 1]

      model = PSRK(tc, pc, w, molecules, c1=C(:, 1), c2=C(:, 2), c3=C(:, 3))
      t = 251.95_pr

      do i=1,10
         bubble = saturation_pressure(model, [x1_exp(i), 1-x1_exp(i)], t, kind="bubble", p0=pbub_exp(i))
         call check(error, rel_error(bubble%p,pbub_exp(i)) < 2e-2)
      end do

      do j = 1,12
         dew = saturation_pressure(model, [y1_exp(j), 1-y1_exp(j)], t, kind="dew", p0=pdew_exp(j))
         call check(error, rel_error(dew%p,pdew_exp(j)) < 1e-1)
      end do
   end subroutine test_psrk_paper2005_4
end module test_psrk
