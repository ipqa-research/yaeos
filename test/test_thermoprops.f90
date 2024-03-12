module test_thermoprops
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use auxiliar_functions, only: rel_error
    implicit none

contains
    subroutine collect_suite(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("Fugacity VT", test_fugacity_VT), &
            new_unittest("Fugacity tp", test_fugacity_tp) &
            ]
    end subroutine collect_suite

    subroutine test_fugacity_VT(error)
        use fixtures_models, only: binary_PR76
        use yaeos, only: pr, CubicEoS, fugacity_vt
        type(error_type), allocatable, intent(out) :: error
        type(CubicEoS) :: eos

        real(pr) :: lnfug(2), dlnphidp(2), dlnphidt(2), dlnphidn(2, 2)

        real(pr), allocatable :: z(:)
        real(pr) ::  v, t, p

        real(pr) :: lnfug_val(2), dlnphidp_val(2), dlnphidt_val(2)

        lnfug_val = [2.0759140949373416, -2.2851989270402058]
        dlnphidp_val = [-0.99059224575177762, -0.99388122357848807]
        dlnphidt_val = [3.0263769083149254E-002, 7.6204871541712640E-002]

        eos = binary_PR76()

        z = [0.3_pr, 0.7_pr]
        v = 8.8780451065729321E-002_pr
        t = 150

        call fugacity_vt(eos, &
            z, V, T, P, lnfug, dlnPhidP, dlnphidT, dlnPhidn &
            )

        call check( &
            error, maxval(abs(lnfug - lnfug_val)) < 1e-5 &
            )
        call check( &
            error, maxval(abs(dlnphidp - dlnphidp_val)) < 1e-5 &
            )
        call check( &
            error, maxval(abs(dlnphidt - dlnphidt_val)) < 1e-5 &
            )
    end subroutine test_fugacity_VT

    subroutine test_fugacity_TP(error)
        use fixtures_models, only: binary_PR76
        use yaeos, only: pr, CubicEoS, fugacity_tp
        type(error_type), allocatable, intent(out) :: error
        type(CubicEoS) :: eos

        real(pr) :: lnfug(2), dlnphidp(2), dlnphidt(2), dlnphidn(2, 2)

        real(pr), allocatable :: z(:)
        real(pr) ::  v, t, p

        real(pr) :: lnfug_val(2), dlnphidp_val(2), dlnphidt_val(2)

        character(len=:), allocatable :: root_type

        lnfug_val = [2.0759140949373416, -2.2851989270402058]
        dlnphidp_val = [-0.99059224575177762, -0.99388122357848807]
        dlnphidt_val = [3.0263769083149254E-002, 7.6204871541712640E-002]

        eos = binary_PR76()
        z = [0.3, 0.7]

        p = 1
        t = 150

        root_type = "liquid"

        call fugacity_tp(eos, &
            z, T, P, V, root_type, lnfug, dlnPhidP, dlnphidT, dlnPhidn &
            )

        call check( &
            error, maxval(rel_error(lnfug_val, lnfug)) < 1e-4 &
            )
        call check( &
            error, maxval(rel_error(dlnphidp_val, dlnphidp)) < 1e-4 &
            )
        call check( &
            error, maxval(rel_error(dlnphidt_val, dlnphidt)) < 1e-4 &
            )
    end subroutine test_fugacity_TP
end module test_thermoprops
