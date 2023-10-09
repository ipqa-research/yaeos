program main
    use yaeos_constants, only: pr, R
    use yaeos_thermo_properties, only: ln_phi
    use fixture_set_pr76, only: set_binary_system, set_binary_system_legacy
    use legacy_ar_models, only: zTVTERMO
    implicit none

    integer, parameter :: n = 2
    real(pr) :: z(n), t, v, yaeos_lnphi(n)

    real(pr) :: p, dpv, philog(n), dlphip(n), dlphit(n), fugn(n, n)

    call set_binary_system
    call set_binary_system_legacy

    z = [0.3, 0.7]
    v = 1.0
    t = 150

    call ln_phi(z, v, t, yaeos_lnphi)
    call zTVTERMO(2, 4, t, z, v, p, dpv, philog, dlphip, dlphit, fugn)

    philog = philog - log(p)

    print *, yaeos_lnphi
    print *, philog
    if (any(abs(yaeos_lnphi - philog) > 0.0000001_pr)) call exit(1)
end program main
