module fixture_set_pr76
        use yaeos_constants, only: pr
        use pengrobinson76, only: setup_pr76
        use yaeos_cubic_eos, only: ac, b, wmod => w, k
        use yaeos_thermo_properties, only: ln_phi
        implicit none

contains

    subroutine set_binary_system
        integer, parameter :: n = 2
        integer :: i

        real(pr) :: tc(n), pc(n), w(n), z(n), kij(n, n), lij(n, n)
        real(pr) :: t, v, lnphi(n)

        z  = [0.3, 0.7]
        tc = [190, 310]
        pc = [14, 30]
        w  = [0.001, 0.03]
        kij = 0
        lij = 0
        call setup_pr76(n, tc, pc, w, kij, lij)
    end subroutine

    subroutine set_binary_system_legacy
    use legacy_ar_models, only: setup, fact => PR76_factory, &
                                ar_srkpr, ArVnder, &
                                kij, lij, ac, b, wmod => w, k
    integer, parameter :: n = 2

    real(pr) :: tc(n), pc(n), w(n), z(n) 
    real(pr) :: t, p, v
    real(pr) :: ar, arv, artv, arv2, arn(n), arvn(n), artn(n), arn2(n,n), et, st

    integer :: i

    z = [0.3, 0.7]
    tc = [190, 310]
    pc = [14, 30]
    w = [0.001, 0.03]
    
    call setup(n, 1, 0, 1)
    kij = 0
    lij = 0
    call fact(z, tc_in=tc, pc_in=pc, w_in=w)
    end subroutine
end module
