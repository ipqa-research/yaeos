program tester
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type

    use test_legacy, only: suite_legacy => collect_suite
    use test_cubic_alphas, only: suite_alphas => collect_suite
    use test_thermoprops, only: suite_thermoprops => collect_suite
    use test_flash, only: suite_flash => collect_suite

    use stdlib_ansi, only : fg_color_green, fg_color_red, operator(//), style_reset
    
    implicit none

    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("legacy", suite_legacy ),  &
        new_testsuite("Alphas", suite_alphas), &
        new_testsuite("Thermoprops", suite_thermoprops), &
        new_testsuite("Flash", suite_flash) &
    ]

    do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    else 
        write(*, *) fg_color_green // "All tests run!" // style_reset
    end if

end program tester