module test_pr76
   use yaeos, only: pr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use auxiliar_functions, only: allclose, rel_error
   implicit none

contains
   subroutine collect_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("Test PR76 Ar consistency mix", test_pr76_cons_mixture), &
         new_unittest("Test PR76 Ar consistency pure", test_pr76_cons_pure), &
         new_unittest("Test PR76 Z", test_pr76_compressibility_factor), &
         new_unittest("Test PR76 CO2 isotherm", test_pr76_co2_volume), &
         new_unittest("Test PR76 Fugacities Elliot", test_pr76_fugacities), &
         new_unittest("Test PR76 Txy methanol-benzene Elliot", test_pr76_txy_methanol_benzene), &
         new_unittest("Test PR76 Pxy nitrogen-methane Elliot", test_pr76_pxy_nitrogen_methane), &
         new_unittest("Test PR76 ethane-heptane envelope Elliot", test_pr76_envelope_ethane_heptane) &
         ]
   end subroutine collect_suite

   subroutine test_pr76_cons_mixture(error)
      use yaeos, only: pr, PengRobinson76, ArModel
      use yaeos__consistency, only: numeric_ar_derivatives, ar_consistency
      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: model, model_kij
      real(pr) :: tc(4), pc(4), w(4)

      real(pr) :: n(4), t, v

      real(pr) :: Ar, ArV, ArT, Arn(size(n)), ArV2, ArT2, ArTV
      real(pr) :: ArVn(size(n)), ArTn(size(n)), Arn2(size(n), size(n))

      real(pr) :: Ar_num, ArV_num, ArT_num, Arn_num(size(n)), ArV2_num, ArT2_num
      real(pr) :: ArTV_num, ArVn_num(size(n)), ArTn_num(size(n))
      real(pr) :: Arn2_num(size(n), size(n))
      real(pr) :: eq31, eq33(size(n), size(n)), eq34(size(n)), eq36, eq37

      real(pr) :: kij(size(n), size(n)), lij(size(n), size(n))

      n = [1.5, 0.2, 0.7, 2.3]
      tc = [190.564, 425.12, 300.11, 320.25]
      pc = [45.99, 37.96, 39.23, 40.21]
      w = [0.0115478, 0.200164, 0.3624, 0.298]

      t = 600_pr
      v = 0.5_pr

      ! ========================================================================
      ! Model without kij and lij
      ! ------------------------------------------------------------------------
      model = PengRobinson76(tc, pc, w)

      call model%residual_helmholtz(&
         n, v, t, Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn, &
         ArTV=ArTV, ArV2=ArV2, ArT2=ArT2, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2)

      call numeric_ar_derivatives(&
         model, n, v, t, d_n = 0.0001_pr, d_v = 0.00001_pr, d_t = 0.001_pr, &
         Ar=Ar_num, ArV=ArV_num, ArT=ArT_num, ArTV=ArTV_num, ArV2=ArV2_num, &
         ArT2=ArT2_num, Arn=Arn_num, ArVn=ArVn_num, ArTn=ArTn_num, &
         Arn2=Arn2_num &
         )

      call ar_consistency(&
         model, n, v, t, eq31=eq31, eq33=eq33, eq34=eq34, eq36=eq36, eq37=eq37 &
         )

      ! Numeric derivatives
      call check(error, rel_error(Ar, Ar_num) < 1e-6)
      call check(error, rel_error(ArV, ArV_num) < 1e-6)
      call check(error, rel_error(ArT, ArT_num) < 1e-6)
      call check(error, allclose(Arn, Arn_num, 1e-6_pr))
      call check(error, rel_error(ArV2, ArV2_num) < 1e-6)
      call check(error, rel_error(ArT2, ArT2_num) < 1e-4)
      call check(error, rel_error(ArTV, ArTV_num) < 1e-6)
      call check(error, allclose(ArVn, ArVn_num, 1e-6_pr))
      call check(error, allclose(ArTn, ArTn_num, 1e-6_pr))
      call check(error, maxval(rel_error(Arn2, Arn2_num)) < 1e-5)

      ! Consistency tests
      call check(error, abs(eq31) <= 1e-15)
      call check(error, maxval(abs(eq33)) < 1e-15)
      call check(error, maxval(abs(eq34)) < 1e-15)
      call check(error, abs(eq36) <= 1e-15)
      call check(error, abs(eq37) <= 1e-15)

      ! ========================================================================
      ! Model with kij and lij
      ! ------------------------------------------------------------------------
      kij = reshape([&
         0.0_pr, 0.1_pr, 0.2_pr, 0.1_pr, &
         0.1_pr, 0.0_pr, 0.3_pr, 0.25_pr, &
         0.2_pr, 0.3_pr, 0.0_pr, 0.18_pr, &
         0.1_pr, 0.25_pr, 0.18_pr, 0.0_pr], [size(n), size(n)])

      lij = reshape([&
         0.0_pr, 0.001_pr, 0.002_pr, 0.001_pr, &
         0.001_pr, 0.0_pr, 0.003_pr, 0.0025_pr, &
         0.002_pr, 0.003_pr, 0.0_pr, 0.0018_pr, &
         0.001_pr, 0.0025_pr, 0.0018_pr, 0.0_pr], [size(n), size(n)])


      model_kij = PengRobinson76(tc, pc, w, kij, lij)

      call model_kij%residual_helmholtz(&
         n, v, t, Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn, &
         ArTV=ArTV, ArV2=ArV2, ArT2=ArT2, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2)

      call numeric_ar_derivatives(&
         model_kij, n, v, t, d_n = 0.0001_pr, d_v = 0.0001_pr, d_t = 0.01_pr, &
         Ar=Ar_num, ArV=ArV_num, ArT=ArT_num, ArTV=ArTV_num, ArV2=ArV2_num, &
         ArT2=ArT2_num, Arn=Arn_num, ArVn=ArVn_num, ArTn=ArTn_num, &
         Arn2=Arn2_num &
         )

      call ar_consistency(&
         model_kij, n, v, t, eq31=eq31, eq33=eq33, &
         eq34=eq34, eq36=eq36, eq37=eq37 &
         )

      ! Numeric derivatives
      call check(error, rel_error(Ar, Ar_num) < 1e-6)
      call check(error, rel_error(ArV, ArV_num) < 1e-6)
      call check(error, rel_error(ArT, ArT_num) < 1e-6)
      call check(error, allclose(Arn, Arn_num, 1e-6_pr))
      call check(error, rel_error(ArV2, ArV2_num) < 1e-6)
      call check(error, rel_error(ArT2, ArT2_num) < 1e-5)
      call check(error, rel_error(ArTV, ArTV_num) < 1e-6)
      call check(error, allclose(ArVn, ArVn_num, 1e-6_pr))
      call check(error, allclose(ArTn, ArTn_num, 1e-6_pr))
      call check(error, maxval(rel_error(Arn2, Arn2_num)) < 1e-6)

      ! Consistency tests
      call check(error, abs(eq31) <= 1e-14)
      call check(error, maxval(abs(eq33)) < 1e-15)
      call check(error, maxval(abs(eq34)) < 1e-14)
      call check(error, abs(eq36) <= 1e-15)
      call check(error, abs(eq37) <= 1e-15)
   end subroutine test_pr76_cons_mixture

   subroutine test_pr76_cons_pure(error)
      use yaeos, only: pr, PengRobinson76, ArModel
      use yaeos__consistency, only: numeric_ar_derivatives, ar_consistency
      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: model
      real(pr) :: tc(1), pc(1), w(1)

      real(pr) :: n(1), t, v

      real(pr) :: Ar, ArV, ArT, Arn(size(n)), ArV2, ArT2, ArTV
      real(pr) :: ArVn(size(n)), ArTn(size(n)), Arn2(size(n), size(n))

      real(pr) :: Ar_num, ArV_num, ArT_num, Arn_num(size(n)), ArV2_num, ArT2_num
      real(pr) :: ArTV_num, ArVn_num(size(n)), ArTn_num(size(n))
      real(pr) :: Arn2_num(size(n), size(n))
      real(pr) :: eq31, eq33(size(n), size(n)), eq34(size(n)), eq36, eq37

      n = [5.0]
      tc = [190.564]
      pc = [45.99]
      w = [0.0115478]

      model = PengRobinson76(tc, pc, w)

      t = 600_pr
      v = 0.5_pr

      call model%residual_helmholtz(&
         n, v, t, Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn, &
         ArTV=ArTV, ArV2=ArV2, ArT2=ArT2, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2)

      call numeric_ar_derivatives(&
         model, n, v, t, d_n = 0.0001_pr, d_v = 0.0001_pr, d_t = 0.01_pr, &
         Ar=Ar_num, ArV=ArV_num, ArT=ArT_num, ArTV=ArTV_num, ArV2=ArV2_num, &
         ArT2=ArT2_num, Arn=Arn_num, ArVn=ArVn_num, ArTn=ArTn_num, &
         Arn2=Arn2_num &
         )

      call ar_consistency(&
         model, n, v, t, eq31=eq31, eq33=eq33, eq34=eq34, eq36=eq36, eq37=eq37 &
         )

      ! Numeric derivatives
      call check(error, rel_error(Ar, Ar_num) < 1e-5)
      call check(error, rel_error(ArV, ArV_num) < 1e-5)
      call check(error, rel_error(ArT, ArT_num) < 1e-5)
      call check(error, allclose(Arn, Arn_num, 1e-5_pr))
      call check(error, rel_error(ArV2, ArV2_num) < 1e-5)
      call check(error, rel_error(ArT2, ArT2_num) < 1e-5)
      call check(error, rel_error(ArTV, ArTV_num) < 1e-5)
      call check(error, allclose(ArVn, ArVn_num, 1e-5_pr))
      call check(error, allclose(ArTn, ArTn_num, 1e-5_pr))
      call check(error, maxval(rel_error(Arn2, Arn2_num)) < 1e-5)

      ! Consistency tests
      call check(error, abs(eq31) <= 1e-14)
      call check(error, maxval(abs(eq33)) < 1e-15)
      call check(error, maxval(abs(eq34)) < 1e-15)
      call check(error, abs(eq36) <= 1e-15)
      call check(error, abs(eq37) <= 1e-15)
   end subroutine test_pr76_cons_pure

   subroutine test_pr76_compressibility_factor(error)
      ! From original paper.
      use yaeos, only: pr, R, PengRobinson76, ArModel
      use yaeos__thermoprops, only: volume

      type(error_type), allocatable, intent(out) :: error
      integer :: i

      class(ArModel), allocatable :: model

      real(pr) :: T, P, v, Zcomp, z(2)
      real(pr) :: P1(5), Z1(5)
      real(pr) :: P2(4), Z2(4)
      real(pr) :: P3(5), Z3(5)
      real(pr) :: P4(4), Z4(4)
      real(pr) :: P5(5), Z5(5)
      real(pr) :: P6(5), Z6(5)

      model = PengRobinson76(&
         [425.12_pr, 304.21_pr], &
         [37.96_pr, 73.83_pr], &
         [0.200164_pr, 0.223621_pr], &
         kij=reshape([0.0_pr, 0.130_pr, 0.130_pr, 0.0_pr], shape=[2,2]) &
         )

      ! =======================================================================
      ! Composition 0.9
      ! -----------------------------------------------------------------------
      T = 310.928_pr ! 100 F
      z = [0.9_pr, 0.1_pr]
      P1 = [41.3685_pr, 68.9476_pr, 137.8951_pr, 206.8427_pr, 275.7903_pr]
      Z1 = [0.151_pr, 0.248_pr, 0.482_pr, 0.707_pr, 0.926_pr]

      do i=1,5
         call volume(model, z, P1(i), T, V=v, root_type="stable")
         Zcomp = P1(i) * v / (R * T)

         call check(error, abs(Zcomp - Z1(i)) < 1e-3)
      end do

      T = 410.928_pr ! 280 F
      z = [0.9_pr, 0.1_pr]
      P2 = [68.9476_pr, 137.8951_pr, 206.8427_pr, 275.7903_pr]
      Z2 = [0.289_pr, 0.482_pr, 0.665_pr, 0.840_pr]

      do i=1,4
         call volume(model, z, P2(i), T, V=v, root_type="stable")
         Zcomp = P2(i) * v / (R * T)

         call check(error, abs(Zcomp - Z2(i)) < 1e-3)
      end do

      T = 510.928_pr ! 460 F
      z = [0.9_pr, 0.1_pr]
      P3 = [41.3685_pr, 68.9476_pr, 137.8951_pr, 206.8427_pr, 275.7903_pr]
      Z3 = [0.804_pr, 0.696_pr, 0.643_pr, 0.744_pr, 0.869_pr]

      do i=1,5
         call volume(model, z, P3(i), T, V=v, root_type="stable")
         Zcomp = P3(i) * v / (R * T)

         call check(error, abs(Zcomp - Z3(i)) < 1e-3)
      end do

      ! =======================================================================
      ! Composition 0.5
      ! -----------------------------------------------------------------------
      T = 310.928_pr ! 100 F
      z = [0.5_pr, 0.5_pr]
      P4 = [68.9476_pr, 137.8951_pr, 206.8427_pr, 275.7903_pr]
      Z4 = [0.215_pr, 0.404_pr, 0.580_pr, 0.750_pr]

      do i=1,4
         call volume(model, z, P4(i), T, V=v, root_type="stable")
         Zcomp = P4(i) * v / (R * T)

         call check(error, abs(Zcomp - Z4(i)) < 1e-3)
      end do

      T = 410.928_pr ! 280 F
      z = [0.5_pr, 0.5_pr]
      P5 = [41.3685_pr, 68.9476_pr, 137.8951_pr, 206.8427_pr, 275.7903_pr]
      Z5 = [0.782_pr, 0.638_pr, 0.545_pr, 0.645_pr, 0.765_pr]

      do i=1,5
         call volume(model, z, P5(i), T, V=v, root_type="stable")
         Zcomp = P5(i) * v / (R * T)

         call check(error, abs(Zcomp - Z5(i)) < 1e-3)
      end do

      T = 510.928_pr ! 460 F
      z = [0.5_pr, 0.5_pr]
      P6 = [41.3685_pr, 68.9476_pr, 137.8951_pr, 206.8427_pr, 275.7903_pr]
      Z6 = [0.920_pr, 0.870_pr, 0.796_pr, 0.806_pr, 0.877_pr]

      do i=1,5
         call volume(model, z, P6(i), T, V=v, root_type="stable")
         Zcomp = P6(i) * v / (R * T)

         call check(error, abs(Zcomp - Z6(i)) < 2e-2)
      end do
   end subroutine test_pr76_compressibility_factor

   subroutine test_pr76_co2_volume(error)
      ! From Elliot's book.
      use yaeos, only : pr, PengRobinson76, ArModel
      use yaeos__thermoprops, only: volume
      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: model

      real(pr) :: mw, V, n(1)
      integer :: i

      model = PengRobinson76([304.2_pr], [73.82_pr], [0.228_pr])
      n = [1.0_pr]
      mw = 44.01 ! g / mol

      call volume(model, n, 8.0_pr, 310.0_pr, V=V, root_type="stable")
      call check(error, abs(V / mw * 1000 - 70.37) < 0.06)

      call volume(model, n, 75.0_pr, 310.0_pr, V=V, root_type="stable")
      call check(error, abs(V / mw * 1000 - 3.84) < 0.01)
   end subroutine test_pr76_co2_volume

   subroutine test_pr76_fugacities(error)
      ! K values of N2-CH4 (0.5, 0.5) mixture from Elliot's book.
      use yaeos, only: pr, R, PengRobinson76, ArModel
      use yaeos__thermoprops, only: fugacity_tp, volume
      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: model

      real(pr) :: T, P, z_v(2), z_l(2), v_v, lnphip_l(2), lnphip_v(2)

      T = 100 ! K
      P = 4.119 ! bar
      z_v = [0.958_pr, 1.0_pr - 0.958_pr]
      z_l = [0.5_pr, 0.5_pr]

      model = PengRobinson76(&
         [126.1_pr, 190.6_pr], &
         [33.94_pr, 46.04_pr], &
         [0.040_pr, 0.011_pr] &
         )

      call volume(model, z_v, P, T, root_type="vapor", V=v_v)
      call fugacity_tp(model, z_v, T, P, root_type="vapor", lnphip = lnphip_v)
      call fugacity_tp(model, z_l, T, P, root_type="liquid", lnphip = lnphip_l)

      ! Elliot Z value of vapor
      call check(error, abs(P * v_v / R / T - 0.9059) <  1e-4)

      ! Elliot vapor fugacities
      call check(error, abs(exp(lnphip_v(1) - log(P)) - 0.9162) < 1e-4)
      call check(error, abs(exp(lnphip_v(2) - log(P)) - 0.8473) < 1e-4)

      ! Elliot liquid fugacities
      call check(error, abs(exp(lnphip_l(1) - log(P)) - 1.791) < 1e-3)
      call check(error, abs(exp(lnphip_l(2) - log(P)) - 0.0937) < 1e-4)
   end subroutine test_pr76_fugacities

   subroutine test_pr76_txy_methanol_benzene(error)
      ! Txy methanol-benzene from Elliot's book using saturation_temperature
      ! function.
      use yaeos, only: pr, PengRobinson76, ArModel, EquilibriaState
      use yaeos, only: saturation_temperature
      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: model
      type(EquilibriaState) :: sat_t

      real(pr) :: z(2), error_sum
      real(pr) :: x_bubble_nokij(18), t_bubble_nokij(18)
      real(pr) :: y_dew_nokij(18), t_dew_nokij(18)
      real(pr) :: x_bubble_kij(20), t_bubble_kij(20)
      real(pr) :: y_dew_kij(20), t_dew_kij(20)

      integer :: i

      ! =======================================================================
      ! No kij simulation points
      ! -----------------------------------------------------------------------
      x_bubble_nokij = [&
         0.010, 0.036, 0.071, 0.118, 0.163, 0.226, &
         0.282, 0.341, 0.410, 0.481, 0.546, 0.626, &
         0.720, 0.796, 0.864, 0.903, 0.948, 0.984  &
         ]

      t_bubble_nokij = [&
         352.584, 351.494, 350.013, 348.455, 347.052, 345.26, &
         343.857, 342.766, 341.675, 340.584, 339.727, 338.87, &
         338.247, 337.857, 337.701, 337.701, 337.701, 337.779 &
         ]

      y_dew_nokij = [&
         0.010, 0.032, 0.058, 0.097, 0.137, 0.182, &
         0.231, 0.270, 0.370, 0.442, 0.492, 0.533, &
         0.583, 0.649, 0.713, 0.809, 0.864, 0.967  &
         ]

      t_dew_nokij = [&
         352.896, 352.429, 351.883, 351.182, 350.325, 349.390, &
         348.377, 347.597, 345.571, 344.013, 343.000, 342.221, &
         341.286, 340.117, 339.104, 338.013, 337.701, 337.623  &
         ]

      ! =======================================================================
      ! kij simulation points
      ! -----------------------------------------------------------------------
      x_bubble_kij = [&
         0.005, 0.017, 0.042, 0.073, 0.098, 0.135, 0.161, 0.222, &
         0.279, 0.348, 0.396, 0.436, 0.491, 0.559, 0.655, 0.688, &
         0.764, 0.847, 0.904, 0.980&
         ]

      t_bubble_kij = [&
         352.412, 350.876, 347.884, 345.054, 343.194, 340.606, &
         338.908, 336.402, 334.461, 333.248, 332.844, 332.520, &
         331.954, 331.469, 331.388, 331.307, 331.388, 331.873, &
         332.439, 335.836 &
         ]

      y_dew_kij = [&
         0.029, 0.082, 0.146, 0.201, 0.236, 0.286, 0.345, 0.392, &
         0.462, 0.529, 0.568, 0.607, 0.652, 0.688, 0.744, 0.793, &
         0.845, 0.880, 0.927, 0.977 &
         ]

      t_dew_kij = [&
         352.251, 350.795, 348.854, 347.156, 346.024, 344.488, &
         342.466, 340.849, 338.342, 335.916, 334.542, 333.248, &
         332.035, 331.388, 332.278, 333.329, 334.461, 335.189, &
         336.240, 337.372 &
         ]

      ! methanol - benzene no kij
      model = PengRobinson76(&
         [512.5_pr, 562.05_pr], &
         [80.84_pr, 48.95_pr], &
         [0.565831_pr, 0.2103_pr] &
         )

      ! bubble no kij
      error_sum = 0.0_pr

      do i = 1, 18
         z = [x_bubble_nokij(i), 1.0_pr - x_bubble_nokij(i)]
         sat_t = saturation_temperature(&
            model, z, p=1.01325_pr, kind="bubble", t0=t_bubble_nokij(i) &
            )

         error_sum = error_sum + (sat_t%T - t_bubble_nokij(i))**2
      end do
      call check(error, sqrt(error_sum) / 18 < 2e-2)

      ! dew no kij
      error_sum = 0.0_pr

      do i = 1, 18
         z = [y_dew_nokij(i), 1.0_pr - y_dew_nokij(i)]
         sat_t = saturation_temperature(&
            model, z, p=1.01325_pr, kind="dew", t0=t_dew_nokij(i)&
            )

         error_sum = error_sum + (sat_t%T - t_bubble_nokij(i))**2
      end do
      call check(error, sqrt(error_sum) / 18 < 0.7)

      ! methanol - benzene kij
      model = PengRobinson76(&
         [512.5_pr, 562.05_pr], &
         [80.84_pr, 48.95_pr], &
         [0.565831_pr, 0.2103_pr], &
         kij=reshape([0.0_pr, 0.084_pr, 0.084_pr, 0.0_pr], shape=[2, 2]) &
         )

      ! bubble kij
      error_sum = 0.0_pr

      do i = 1, 20
         z = [x_bubble_kij(i), 1.0_pr - x_bubble_kij(i)]
         sat_t = saturation_temperature(&
            model, z, p=1.01325_pr, kind="bubble", t0=t_bubble_kij(i) &
            )

         error_sum = error_sum + (sat_t%T - t_bubble_kij(i))**2
      end do
      call check(error, sqrt(error_sum) / 20 < 6e-2)

      ! Dew kij
      error_sum = 0.0_pr

      do i = 1, 20
         z = [y_dew_kij(i), 1.0_pr - y_dew_kij(i)]
         sat_t = saturation_temperature(&
            model, z, p=1.01325_pr, kind="dew", t0=t_dew_kij(i) &
            )

         error_sum = error_sum + (sat_t%T - t_bubble_kij(i))**2
      end do
      call check(error, sqrt(error_sum) / 20 < 0.6)
   end subroutine test_pr76_txy_methanol_benzene

   subroutine test_pr76_pxy_nitrogen_methane(error)
      ! Pxy nitrogen-methane from Elliot's book using saturation_pressure
      ! function.
      use yaeos, only: pr, PengRobinson76, ArModel, EquilibriaState
      use yaeos, only: saturation_pressure

      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: model
      type(EquilibriaState) :: psat
      real(pr) :: T, y_dew(9), p_dew(9), x_bubble(10), p_bubble(10), z(2)
      real(pr) :: error_sum

      integer :: i

      y_dew = [0.05545, 0.1349, 0.2237, 0.3420, 0.4344, 0.5287, 0.6118, 0.6765, 0.7006]
      p_dew = [10.97, 12.71, 13.96, 16.45, 19.20, 23.43, 29.67, 35.65, 40.89]

      x_bubble = [0.06839, 0.1257, 0.1848, 0.2791, 0.2292, 0.3475, 0.4122, 0.4621, 0.5379, 0.6100]
      p_bubble = [14.71, 18.20, 21.44, 26.68, 24.18, 30.17, 33.91, 36.90, 40.89, 44.63]

      T = 150.0_pr

      model = PengRobinson76(&
         [126.2_pr, 190.564_pr], &
         [34.0_pr, 45.99_pr], &
         [0.0377215_pr, 0.0115478_pr] &
         )

      error_sum = 0.0_pr
      do i = 1,9
         z = [y_dew(i), 1.0_pr - y_dew(i)]
         psat = saturation_pressure(model, z, T, kind="dew", p0=p_dew(i))

         error_sum = error_sum + (psat%P - p_dew(i))**2
      end do
      call check(error, sqrt(error_sum) / 9 < 0.15)

      error_sum = 0.0_pr
      do i = 1,10
         z = [x_bubble(i), 1.0_pr - x_bubble(i)]
         psat = saturation_pressure(model, z, T, kind="bubble", p0=p_bubble(i))

         error_sum = error_sum + (psat%P - p_bubble(i))**2
      end do
      call check(error, sqrt(error_sum) / 10 < 0.07)
   end subroutine test_pr76_pxy_nitrogen_methane

   subroutine test_pr76_envelope_ethane_heptane(error)
      use yaeos, only: pr, ArModel, EquilibriaState, PTEnvel2, PengRobinson76
      use yaeos, only: saturation_pressure, pt_envelope_2ph

      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: model
      type(EquilibriaState) :: bubble
      type(PTEnvel2) :: envelope

      real(pr) :: z(2)
      real(pr), allocatable :: bubble_points_t(:), dew_points_t(:)
      real(pr), allocatable :: bubble_points_p(:), dew_points_p(:)
      real(pr) :: t_bub(21), p_bub(21), t_dew(14), p_dew(14)
      real(pr) :: p_env, T, y1, y2, y3, x1, x2, x3, term1, term2, term3
      real(pr) :: error_sum
      ! real(pr) :: crit_point_t, crit_point_p

      integer :: i, j, index

      model = PengRobinson76( &
         [305.32_pr, 540.2_pr], &
         [48.72_pr, 27.4_pr], &
         [0.099493_pr, 0.349469_pr] &
         )

      z = [0.1_pr, 0.9_pr]

      bubble = saturation_pressure(model, z, T=230._pr, kind="bubble", p0=5.0_pr)

      envelope = pt_envelope_2ph(model, z, bubble)

      allocate(bubble_points_t(0))
      allocate(bubble_points_p(0))
      allocate(dew_points_t(0))
      allocate(dew_points_p(0))

      do i = 1, size(envelope%points)
         if (envelope%points(i)%kind .eq. "bubble") then
            bubble_points_t = [bubble_points_t, envelope%points(i)%T]
            bubble_points_p = [bubble_points_p, envelope%points(i)%P]
         elseif (envelope%points(i)%kind .eq. "dew") then
            dew_points_t = [dew_points_t, envelope%points(i)%T]
            dew_points_p = [dew_points_p, envelope%points(i)%P]
         end if
      end do

      ! =======================================================================
      ! Critical point check
      ! -----------------------------------------------------------------------
      call check(error, abs(envelope%cps(1)%T - 532.9_pr) < 1)
      call check(error, abs(envelope%cps(1)%P - 32.95_pr) < 0.2)

      ! =======================================================================
      ! Bubble points check
      ! -----------------------------------------------------------------------
      t_bub = [&
         241.0, 252.4, 262.0, 272.3, 283.1, 295.6, 309.8, &
         326.9, 344.0, 361.1, 374.1, 395.2, 408.3, 419.1, &
         433.9, 451.0, 469.2, 484.0, 501.1, 513.6, 524.4  &
         ]
      p_bub = [&
         0.8238, 1.236, 1.442, 1.854, 2.471, 3.089, 3.913, &
         4.943, 6.178, 7.414, 8.650, 10.71, 12.15, 13.59,  &
         15.45, 17.92, 21.01, 23.89, 27.39, 30.07, 32.33   &
         ]

      error_sum = 0.0_pr
      do i = 1, size(t_bub)
         do j = 1, size(bubble_points_t)
            if (bubble_points_t(j) > t_bub(i)) then
               ! lagrange polynomics
               x1 = bubble_points_t(j-1)
               x2 = bubble_points_t(j)
               x3 = bubble_points_t(j+1)

               y1 = bubble_points_p(j-1)
               y2 = bubble_points_p(j)
               y3 = bubble_points_p(j+1)
               exit
            end if
         end do
         T = t_bub(i)
         term1 = y1 * (T - x2) * (T - x3) / ((x1 - x2) * (x1 - x3))
         term2 = y2 * (T - x1) * (T - x3) / ((x2 - x1) * (x2 - x3))
         term3 = y3 * (T - x1) * (T - x2) / ((x3 - x1) * (x3 - x2))

         p_env = term1 + term2 + term3

         error_sum = error_sum + (p_env - p_bub(i))**2
      end do

      call check(error, sqrt(error_sum) / size(t_bub) < 0.03)

      ! =======================================================================
      ! Dew points check
      ! -----------------------------------------------------------------------
      t_dew = [&
         368.5, 387.8, 406.6, 420.8, 437.3, 450.4, 464.1, &
         474.9, 489.1, 497.6, 507.9, 515.9, 524.4, 531.2  &
         ]

      p_dew= [&
         0.8238, 1.854, 2.883, 3.913, 5.561, 7.208, 9.474, &
         11.53, 14.83, 17.09, 20.18, 23.07, 26.57, 30.48   &
         ]

      error_sum = 0.0_pr
      do i = 1, size(t_dew)
         do j = 1, size(dew_points_t)
            if (dew_points_t(j) < t_dew(i)) then
               x1 = dew_points_t(j-1)
               x2 = dew_points_t(j)
               x3 = dew_points_t(j+1)

               y1 = dew_points_p(j-1)
               y2 = dew_points_p(j)
               y3 = dew_points_p(j+1)
               exit
            end if
         end do
         T = t_dew(i)
         term1 = y1 * (T - x2) * (T - x3) / ((x1 - x2) * (x1 - x3))
         term2 = y2 * (T - x1) * (T - x3) / ((x2 - x1) * (x2 - x3))
         term3 = y3 * (T - x1) * (T - x2) / ((x3 - x1) * (x3 - x2))

         p_env = term1 + term2 + term3

         error_sum = error_sum + (p_env - p_dew(i))**2
      end do

      call check(error, sqrt(error_sum) / size(t_dew) < 0.04)
   end subroutine test_pr76_envelope_ethane_heptane
end module test_pr76
