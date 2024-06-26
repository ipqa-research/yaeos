program tester
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type

    use test_legacy, only: suite_legacy => collect_suite
    use test_cubic_alphas, only: suite_alphas => collect_suite
    use test_cubic_implementations, only: suite_implementations => collect_suite
    use test_cubic_mixrules, only: suite_cubic_mixrules => collect_suite
    use test_autodiff_api, only: suite_autodiff_hd => collect_suite
    use test_thermoprops, only: suite_thermoprops => collect_suite
    use test_flash, only: suite_flash => collect_suite
    use test_saturation, only: suite_saturation => collect_suite
    use test_math, only: suite_math => collect_suite

    ! =========================================================================
    ! Implemented ArModels testings
    ! -------------------------------------------------------------------------
    use test_pr76, only: suite_pr76 => collect_suite
    use test_pr78, only: suite_pr78 => collect_suite
    use test_srk, only: suite_srk => collect_suite
    use test_rkpr, only: suite_rkpr => collect_suite

    ! =========================================================================
    ! Implemented GeModels testings
    ! -------------------------------------------------------------------------
    use test_unifac, only: suite_unifac => collect_suite
    use test_unifac_parameters, only: suite_unifac_parameters => collect_suite
    use test_tape_nrtl, only: suite_nrtl => collect_suite

    use stdlib_ansi, only: fg_color_green, fg_color_red, operator(//), style_reset

    implicit none

    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("legacy", suite_legacy), &
        new_testsuite("Alphas", suite_alphas), &
        new_testsuite("Cubic EoS", suite_implementations), &
        new_testsuite("Cubic MixRules", suite_cubic_mixrules), &
        new_testsuite("Autodiff APIs", suite_autodiff_hd), &
        new_testsuite("Thermoprops", suite_thermoprops), &
        new_testsuite("Flash", suite_flash), &
        new_testsuite("Saturation Points", suite_saturation), &
        new_testsuite("Math module", suite_math), &
        ! =====================================================================
        ! Armodel particular tests
        ! ---------------------------------------------------------------------
        ! Cubic models
        new_testsuite("PengRobinson76", suite_pr76), &
        new_testsuite("PengRobinson78", suite_pr78), &
        new_testsuite("SoaveRedlichKwong", suite_srk), &
        new_testsuite("RKPR", suite_rkpr), &
        ! =====================================================================
        ! Ge particular tests
        ! ---------------------------------------------------------------------
        new_testsuite("UNIFAC", suite_unifac), &
        new_testsuite("UNIFACParameters", suite_unifac_parameters), &
        new_testsuite("NRTL", suite_nrtl) &
        ]

    do is = 1, size(testsuites)
        write (error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    else
        write (*, *) fg_color_green//"All tests run!"//style_reset
    end if

end program tester
