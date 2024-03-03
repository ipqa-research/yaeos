program main
    use yaeos_constants, only: pr
    use yaeos_substance, only: Substances
    use yaeos_models, only: PengRobinson76, ArModel
    use yaeos_thermoprops, only: fugacity_tp
    use yaeos_phase_equilibria_stability, only: tangent_plane_distance
    type(Substances) :: compos
    class(ArModel), allocatable :: model

    real(pr) :: tc(2), pc(2), af(2), lnfug_z(2), lnfug_w(2)
    real(pr) :: w(2), z(2), t, p
    real(pr) :: tpd

    integer :: i

    z = [0.4, 0.6]
    tc = [190.564, 425.12]
    pc = [45.99, 37.96]
    af = [0.0115478, 0.200164]
    model = PengRobinson76(tc, pc, af)

    P = 60
    t = 294

    call fugacity_tp(model, z, t, p, root_type="stable", lnfug=lnfug_z)

    lnfug_z = lnfug_z - log(P)
    do i=1,99
        w(1) = real(i, pr)/100
        w(2) = 1 - w(1)
        call fugacity_tp(model, w, t, p, lnfug=lnfug_w, root_type="stable")

        call tangent_plane_distance(z, lnfug_z, w, lnfug_w-log(P), tpd)
        print *, i, tpd
    end do

end program
