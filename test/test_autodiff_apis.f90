module test_autodiff_api
    use yaeos_constants, only: pr
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use auxiliar_functions, only: allclose
    implicit none

    real(pr) :: absolute_tolerance = 1e-6_pr

contains

    subroutine collect_suite(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("hdPR76", test_pr76_hd), &
            new_unittest("tapePR76", test_pr76_tape) &
            ]
    end subroutine collect_suite

    subroutine test_pr76_hd(error)
        use yaeos_constants, only: pr
        use fixtures_models, only: binary_PR76_hd
        use yaeos, only: ArModel
        type(error_type), allocatable, intent(out) :: error

        class(ArModel), allocatable :: eos
        integer, parameter :: n = 2
        real(pr) :: z(n), V, T
        real(pr) :: Ar, ArV, ArV2, ArT, ArTV, ArT2
        real(pr) :: Arn(n), ArVn(n), ArTn(n), Arn2(n, n)

        real(pr) :: Ar_val, ArV_val, ArV2_val, ArT_val, ArTV_val, ArT2_val
        real(pr) :: Arn_val(n), ArVn_val(n), ArTn_val(n), Arn2_val(n, n)

        Ar_val = -9.5079213387597061
        ArV_val = 8.8348105702414230
        ArT_val = 2.5288760006412853E-002
        ArT2_val = -8.1263714911056052E-005
        ArV2_val = -16.452712169871607
        ArTV_val = -2.4354181554918298E-002
        Arn_val = [-14.760083989416412, -19.878152533126190]
        ArVn_val = [12.970846906902654, 17.944940224423746]
        ArTn_val = [4.7299709855544367E-002, 5.0647183777961201E-002]
        Arn2_val(1, :) = [-11.697767407192650, -13.516452437750393]
        Arn2_val(2, :) = [-13.516452437750393, -19.842863669307611]

        eos = binary_PR76_hd()
        z = [0.3, 0.7]
        v = 1
        T = 150
        call eos%residual_helmholtz( &
            z, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
            ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
            )

        call check(error, allclose([Ar], [Ar_val], absolute_tolerance))
        call check(error, allclose([ArV], [ArV_val], absolute_tolerance))
        call check(error, allclose([ArT], [ArT_val], absolute_tolerance))
        call check(error, allclose([ArTV], [ArTV_val], absolute_tolerance))
        call check(error, allclose([ArV2], [ArV2_val], absolute_tolerance))
        call check(error, allclose([ArT2], [ArT2_val], absolute_tolerance))

        call check(error, allclose([ArVn], [ArVn_val], absolute_tolerance))
        call check(error, allclose([ArTn], [ArTn_val], absolute_tolerance))
        call check(error, allclose([Arn2], [Arn2_val], absolute_tolerance))
    end subroutine
    
    subroutine test_pr76_tape(error)
        use yaeos_constants, only: pr
        use fixtures_models, only: binary_PR76_tape
        use yaeos, only: ArModel
        type(error_type), allocatable, intent(out) :: error

        class(ArModel), allocatable :: eos
        integer, parameter :: n = 2
        real(pr) :: z(n), V, T
        real(pr) :: Ar, ArV, ArV2, ArT, ArTV, ArT2
        real(pr) :: Arn(n), ArVn(n), ArTn(n), Arn2(n, n)

        real(pr) :: Ar_val, ArV_val, ArV2_val, ArT_val, ArTV_val, ArT2_val
        real(pr) :: Arn_val(n), ArVn_val(n), ArTn_val(n), Arn2_val(n, n)

        Ar_val = -9.5079213387597061
        ArV_val = 8.8348105702414230
        ArT_val = 2.5288760006412853E-002
        ArT2_val = -8.1263714911056052E-005
        ArV2_val = -16.452712169871607
        ArTV_val = -2.4354181554918298E-002
        Arn_val = [-14.760083989416412, -19.878152533126190]
        ArVn_val = [12.970846906902654, 17.944940224423746]
        ArTn_val = [4.7299709855544367E-002, 5.0647183777961201E-002]
        Arn2_val(1, :) = [-11.697767407192650, -13.516452437750393]
        Arn2_val(2, :) = [-13.516452437750393, -19.842863669307611]

        eos = binary_PR76_tape()
        z = [0.3, 0.7]
        v = 1
        T = 150
        call eos%residual_helmholtz( &
            z, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
            ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn &
            )
        call check(error, allclose([Ar], [Ar_val], absolute_tolerance))
        call check(error, allclose([ArV], [ArV_val], absolute_tolerance))
        call check(error, allclose([ArT], [ArT_val], absolute_tolerance))
        call check(error, allclose([ArTV], [ArTV_val], absolute_tolerance))
        call check(error, allclose([ArV2], [ArV2_val], absolute_tolerance))
        call check(error, allclose([ArT2], [ArT2_val], absolute_tolerance))
        call check(error, allclose([ArVn], [ArVn_val], absolute_tolerance))
        call check(error, allclose([ArTn], [ArTn_val], absolute_tolerance))


        call eos%residual_helmholtz( &
            z, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
            ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
            )

        call check(error, allclose([Ar], [Ar_val], absolute_tolerance))
        call check(error, allclose([ArV], [ArV_val], absolute_tolerance))
        call check(error, allclose([ArT], [ArT_val], absolute_tolerance))
        call check(error, allclose([ArTV], [ArTV_val], absolute_tolerance))
        call check(error, allclose([ArV2], [ArV2_val], absolute_tolerance))
        call check(error, allclose([ArT2], [ArT2_val], absolute_tolerance))

        call check(error, allclose([ArVn], [ArVn_val], absolute_tolerance))
        call check(error, allclose([ArTn], [ArTn_val], absolute_tolerance))
        call check(error, allclose([Arn2], [Arn2_val], absolute_tolerance))
    end subroutine
end module
