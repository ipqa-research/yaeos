program main
    use yaeos_constants, only: pr
    use pengrobinson76, only: setup_pr76
    use yaeos_cubic_eos, only: ac, b, wmod => w, k
    use yaeos_ar_models, only: residual_helmholtz
    use yaeos_thermo_properties, only: get_volume, pressure
    implicit none

    integer, parameter :: n = 2, npoints=10000
    real(pr), parameter :: v0=0.06_pr, vf=20_pr
    real(pr) :: et, st
    
    integer :: i

    real(pr) :: tc(n), pc(n), w(n), z(n), kij(n, n), lij(n, n)
    real(pr) :: t, v, v_vap, p , v_liq, dv, vs(npoints)

    z = [0.3, 0.7]
    tc = [190, 310]
    pc = [14, 30]
    w = [0.001, 0.03]

    kij = 0
    lij = 0

    call setup_pr76(n, tc, pc, w, kij, lij)

    t = 150._pr

    v = v0
    dv = (vf - v0)/real(npoints,pr)
    do i=1,npoints
        vs(i) = v0 + i*dv
    end do

    call cpu_time(st)
    do i=1,npoints
        call pressure(z, 1.0_pr, t, p)
        ! print *, vs(i), p, v_vap, v_liq
    end do
    call cpu_time(et)
    print *, (et-st) * 1e6 / npoints
end program main
