module test_flash
    use testdrive, only: new_unittest, unittest_type, error_type, check
    implicit none

contains
    subroutine collect_suite(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("FlashPT", test_flash_pt) &
            ]
    end subroutine collect_suite

    subroutine test_flash_pt(error)
        use yaeos, only: pr, EquilibriaState, flash, PengRobinson76, ArModel, fugacity_tp
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nc = 2

        real(pr) :: x(nc) = [0.32424471950363210, 0.67575531029866709]
        real(pr) :: y(nc) = [0.91683466155334536, 8.3165368249135715E-002]
        real(pr) :: vx = 8.4918883298198036E-002
        real(pr) :: vy = 0.32922132295944545

        class(ArModel), allocatable :: model
        type(EquilibriaState) :: flash_result

        real(pr) :: tc(nc), pc(nc), w(nc)
        real(pr) :: n(nc), p, t, k0(nc)
        integer :: iters

        n = [0.4, 0.6]
        tc = [190.564, 425.12]
        pc = [45.99, 37.96]
        w = [0.0115478, 0.200164]
        model = PengRobinson76(tc, pc, w)

        P = 60
        t = 294
        k0 = (PC/P)*exp(5.373*(1 + w)*(1 - TC/T))

        flash_result = flash(model, n, t=t, p_spec=p, k0=k0, iters=iters)

        call check(error, maxval(abs(flash_result%x - x)) < 1e-5)
        call check(error, maxval(abs(flash_result%y - y)) < 1e-5)
        call check(error, (abs(flash_result%Vx - Vx)) < 1e-5)
        call check(error, (abs(flash_result%Vy - Vy)) < 1e-5)
    end subroutine test_flash_pt
end module test_flash
