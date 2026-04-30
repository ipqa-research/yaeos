! program simple_gca_test
!    !! Simple validation test for GCA-EoS
!    !! Tests individual contributions separately
!    use yaeos__constants, only: pr, R
!    use yaeos__models_ar_gc_gca
!    implicit none
! 
!    type(GCAEOS) :: eos
!    integer, parameter :: nc = 1
!    real(pr) :: T, V, n(1)
!    real(pr) :: Ar, ArV, ArT
!    integer :: i
! 
!    print *, ""
!    print *, "=========================================="
!    print *, "  GCA-EoS Simple Validation Test"
!    print *, "=========================================="
!    print *, ""
! 
!    ! ========================================================================
!    ! Test 1: Pure component without association
!    ! ========================================================================
!    print *, "Test 1: Pure component (no association)"
!    print *, "------------------------------------------"
! 
!    ! Initialize simple system
!    eos = GCAEOS(nc=1, ng=1, nga=0, nst=0)
! 
!    ! Set properties (simplified n-hexane-like)
!    eos%tc(1) = 507.82_pr
!    eos%pc(1) = 30.34_pr
!    eos%omega(1) = 0.2975_pr
!    eos%dc(1) = 5.95_pr
! 
!    ! Single group representing the whole molecule
!    eos%ny(1, 1) = 1
!    eos%q(1) = 4.5_pr
!    eos%tstr(1) = 600.0_pr
!    eos%gstr(1) = 150000.0_pr
!    eos%g1(1) = 0.8_pr
!    eos%g2(1) = -0.1_pr
!    eos%tspl(1) = 800.0_pr
!    eos%epx(1) = 1.0_pr
!    eos%kstr(1, 1) = 1.0_pr
!    eos%kp(1, 1) = 0.0_pr
!    eos%alpha(1, 1) = 0.0_pr
! 
!    ! Conditions
!    T = 298.15_pr
!    V = 1.0_pr
!    n(1) = 1.0_pr
! 
!    print '(A,F7.2,A)', "  T = ", T, " K"
!    print '(A,F7.3,A)', "  V = ", V, " L"
!    print '(A,F7.3,A)', "  n = ", n(1), " mol"
!    print *, ""
! 
!    ! Update parameters and calculate
!    call eos%residual_helmholtz(n, V, T, Ar=Ar, ArV=ArV, ArT=ArT)
! 
!    print *, "Results:"
!    print '(A,ES15.6,A)', "  Ar  = ", Ar, " bar·L"
!    print '(A,ES15.6,A)', "  P   = ", -ArV, " bar"
!    print '(A,ES15.6,A)', "  Z   = ", -ArV * V / (R * T * n(1))
! 
!    ! Check if results are reasonable
!    if (Ar < 0.0_pr .and. ArV < 0.0_pr) then
!       print *, "  ✓ Signs are correct (Ar < 0, P > 0)"
!    else
!       print *, "  ✗ WARNING: Unexpected signs!"
!    end if
!    print *, ""
! 
!    ! ========================================================================
!    ! Test 2: Temperature sensitivity
!    ! ========================================================================
!    print *, "Test 2: Temperature dependence"
!    print *, "------------------------------------------"
!    print *, "  T (K)      P (bar)        Z"
!    print *, "  ------   ----------   ----------"
! 
!    V = 1.0_pr
!    n(1) = 1.0_pr
! 
!    do i = 1, 5
!       T = 250.0_pr + real(i - 1, pr) * 50.0_pr
!       call eos%residual_helmholtz(n, V, T, ArV=ArV)
!       print '(F8.2,F13.4,F13.6)', T, -ArV, -ArV * V / (R * T * n(1))
!    end do
!    print *, ""
! 
!    ! ========================================================================
!    ! Test 3: Volume/Pressure dependence
!    ! ========================================================================
!    print *, "Test 3: Volume/Pressure dependence"
!    print *, "------------------------------------------"
!    print *, "  V (L)      P (bar)        Z"
!    print *, "  ------   ----------   ----------"
! 
!    T = 298.15_pr
!    n(1) = 1.0_pr
! 
!    do i = 1, 5
!       V = 0.1_pr + real(i - 1, pr) * 0.5_pr
!       call eos%residual_helmholtz(n, V, T, ArV=ArV)
!       print '(F7.3,F13.4,F13.6)', V, -ArV, -ArV * V / (R * T * n(1))
!    end do
!    print *, ""
! 
!    ! ========================================================================
!    ! Test 4: With association
!    ! ========================================================================
!    print *, "Test 4: Component with association (alcohol)"
!    print *, "------------------------------------------"
! 
!    ! Reinitialize with association
!    eos = GCAEOS(nc=1, ng=1, nga=1, nst=2)
! 
!    ! Properties
!    eos%tc(1) = 513.92_pr     ! Ethanol-like
!    eos%pc(1) = 61.48_pr
!    eos%omega(1) = 0.6436_pr
!    eos%dc(1) = 3.45_pr
! 
!    ! Group
!    eos%ny(1, 1) = 1
!    eos%q(1) = 2.5_pr
!    eos%tstr(1) = 500.0_pr
!    eos%gstr(1) = 180000.0_pr
!    eos%g1(1) = 1.0_pr
!    eos%g2(1) = -0.12_pr
!    eos%tspl(1) = 700.0_pr
!    eos%epx(1) = 1.2_pr
!    eos%kstr(1, 1) = 1.0_pr
!    eos%kp(1, 1) = 0.0_pr
!    eos%alpha(1, 1) = 0.0_pr
! 
!    ! Association sites
!    eos%sigma(1, 1) = 1  ! Donor
!    eos%sigma(2, 1) = 1  ! Acceptor
!    eos%eps_r(1, 2) = 2500.0_pr
!    eos%eps_r(2, 1) = 2500.0_pr
!    eos%kappa(1, 2) = 0.02_pr
!    eos%kappa(2, 1) = 0.02_pr
!    eos%eps_r(1, 1) = 0.0_pr
!    eos%eps_r(2, 2) = 0.0_pr
!    eos%kappa(1, 1) = 0.0_pr
!    eos%kappa(2, 2) = 0.0_pr
! 
!    T = 298.15_pr
!    V = 1.0_pr
!    n(1) = 1.0_pr
! 
!    print '(A,F7.2,A)', "  T = ", T, " K"
!    print '(A,F7.3,A)', "  V = ", V, " L"
!    print *, ""
! 
!    call eos%residual_helmholtz(n, V, T, Ar=Ar, ArV=ArV)
! 
!    print *, "Results:"
!    print '(A,ES15.6,A)', "  Ar  = ", Ar, " bar·L"
!    print '(A,ES15.6,A)', "  P   = ", -ArV, " bar"
!    print '(A,ES15.6)', "  Z   = ", -ArV * V / (R * T * n(1))
!    print *, ""
! 
!    print *, "Comparison with/without association:"
!    print *, "(Association should increase pressure)"
!    print *, ""
! 
!    ! ========================================================================
!    ! Test 5: Derivative consistency check
!    ! ========================================================================
!    print *, "Test 5: Numerical derivative check"
!    print *, "------------------------------------------"
! 
!    eos = GCAEOS(nc=3, ng=5, nga=3, nst=8)
! 
!    eos%Tc = [562.0, 536.8, 645.6]
!    eos%Pc = [48.3, 51.1, 41.8]
!    eos%omega = [0.193, 0.251, 0.301]
!    eos%dc = [4.3823, 4.0330, 5.0639]
! 
!    eos%tstr = [600.0, 600.0, 512.6, 600.0, 600.0]
!    eos%gstr = [316910.0, 356080.0, 531330.0, 723210.0, 614185.0]
!    eos%g1 = [-0.9274, -0.8755, -0.3201, -0.6060, -0.5459]
!    eos%g2 = [0.0, 0.0, -0.0168, 0.0, -0.2]
!    eos%tspl = [1246.97, 1285.32, 2076.34, 1590.10, 1497.98]
!    eos%epx = [7.71, 7.50, 5.25, 6.42, 6.25]
!    eos%ny(1, :) = [0, 0, 0, 6, 0]
!    eos%ny(2, :) = [1, 1, 1, 0, 0]
!    eos%ny(3, :) = [0, 0, 0, 5, 1]
! 
!    eos%kstr = reshape([ &
!       1.0, 0.995, 0.980, 0.985, 0.990, &
!       0.995, 1.0, 0.985, 0.990, 0.995, &
!       0.980, 0.985, 1.0, 0.975, 0.980, &
!       0.985, 0.990, 0.975, 1.0, 0.995, &
!       0.990, 0.995, 0.980, 0.995, 1.0], &
!       shape=[5,5])
!    eos%kp = 0
!    eos%kp(1, :) = [0., 0., -0.09, 0.0944, 0.0]
!    eos%kp(2, :) = [0., 0., 0.005, 0.0944, 0.0]
!    eos%kp(3, :) = [-0.09, 0.005, 0.0, 0.0944, 0.0]
!    eos%kp(4, :) = [0.0944, 0.0944, 0.0, 0.0, 0.0]
!    eos%kp(5, :) = [0.0, 0.0, 0.0, 0.0, 0.0]
! 
!    eos%alpha = 0
!    eos%alpha(1, :) = [0.0, 0.0, 0.0, 0.3915, 0.0]
!    eos%alpha(2, :) = [0.0, 0.0, 0.0, 0.3915, 0.0]
!    eos%alpha(3, :) = [0.0, 0.0, 0.0, 0.0, 0.0]
!    eos%alpha(4, :) = [0.3915, 0.3915, 0.0, 0.0, 0.0]
!    eos%alpha(5, :) = [0.0, 0.0, 0.0, 0.0, 0.0]
! 
! 
!  Number of associating groups =  3
! 
!  Associating group configuration 
!                                   Compound No
!    No.  Name       No of sites    1    2    3
!  --------------------------------------------
!     3   -OH_1º          2         0    1    0
!    16   o-AOH           2         0    0    1
!    20   C=C             1         1    0    1
! 
! 
! 
! contains
! 
!    subroutine test_derivatives(eos, n, V, T)
!       use yaeos__consistency_armodel, only: numeric_ar_derivatives
!       !! Test derivative consistency using finite differences
!       type(GCAEOS), intent(inout) :: eos
!       real(pr), intent(in) :: n(:), V, T
! 
! 
! 
!       real(pr) :: Ar1, Ar2, ArV_analytic, ArV_numeric
!       real(pr) :: ArT_numeric, ArV2_numeric, ArT2_numeric, Arn_numeric(nc)
!       real(pr) :: ArTV_numeric, ArVn_numeric(nc), ArTn_numeric(nc), Arn2_numeric(nc, nc)
! 
!       real(pr) :: ArT_analytic, ArV2_analytic, ArT2_analytic, Arn_analytic(nc)
!       real(pr) :: ArTV_analytic, ArVn_analytic(nc), ArTn_analytic(nc), Arn2_analytic(nc, nc)
! 
!       real(pr) :: error
!       real(pr), parameter :: tol = 1.0e-6_pr
! 
!       call numeric_ar_derivatives(&
!          eos, n, V, T, 1e-5_pr, 1e-5_pr, 1e-5_pr, Ar, &
!          ArV_numeric, ArT_numeric, Arn_numeric, &
!          ArV2_numeric, ArT2_numeric, ArTV_numeric, &
!          ArVn_numeric, ArTn_numeric, Arn2_numeric)
! 
!       call eos%residual_helmholtz(n, V, T, Ar=Ar1, ArV=ArV_analytic, ArT=ArT_analytic, &
!          ArV2=ArV2_analytic, ArT2=ArT2_analytic, ArTV=ArTV_analytic, &
!          ArVn=ArVn_analytic, ArTn=ArTn_analytic, Arn2=Arn2_analytic)
! 
! 
! 
!       print *, "Testing analytical vs numerical derivatives:"
!       print *, ""
! 
!       ! Test dAr/dV
!       error = abs(ArV_analytic - ArV_numeric)
! 
!       print '(A,ES15.6)', "  ArV (analytical): ", ArV_analytic
!       print '(A,ES15.6)', "  ArV (numerical):  ", ArV_numeric
!       print '(A,ES12.4)', "  Relative error:   ", error
! 
!       if (error < tol) then
!          print *, "  ✓ Volume derivative OK"
!       else
!          print *, "  ✗ WARNING: Large error in volume derivative"
!       end if
!       print *, ""
! 
!       ! Test dAr/dT
!       call eos%residual_helmholtz(n, V, T, Ar=Ar1, ArT=ArT_analytic)
!       error = abs((ArT_analytic - ArT_numeric))
! 
!       print '(A,ES15.6)', "  ArT (analytical): ", ArT_analytic
!       print '(A,ES15.6)', "  ArT (numerical):  ", ArT_numeric
!       print '(A,ES12.4)', "  Relative error:   ", error
! 
!       if (error < tol) then
!          print *, "  ✓ Temperature derivative OK"
!       else
!          print *, "  ✗ WARNING: Large error in temperature derivative"
!       end if
! 
!    end subroutine test_derivatives
! 
! end program simple_gca_test
! 