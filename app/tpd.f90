program main
    use yaeos_constants, only: pr
    use yaeos_substance, only: Substances
    use yaeos_models, only: PengRobinson76, ArModel
    use yaeos_thermoprops, only: fugacity_tp
    use yaeos_phase_equilibria_stability, only: tangent_plane_distance
    type(Substances) :: compos
    class(ArModel), allocatable :: model
    integer, parameter :: n = 2

    real(pr) :: tc(n), pc(n), af(n), lnfug_z(n), lnfug_w(n)
    real(pr) :: w(n), z(n), t, p
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

    opt: block
        use powellopt, only: bobyqa
        real(pr) :: wl(2), wu(2)
        integer :: npt

        npt = (n+2 + (n+1)*(n+2)/2) / 2
        wl = 1.e-5_pr
        wu = 1
        call bobyqa(n, npt, w, wl, wu, 0.01_pr, 1.e-5_pr, 0, 100, func)
        print *, w/sum(w)
    end block opt

contains

    subroutine func(n, x, f)
        integer, intent(in) :: n
        real(pr), intent(in) :: x(:)
        real(pr), intent(out) :: f

        real(pr), dimension(size(x)) :: w, lnfug_w

        w = x/sum(x)
        call fugacity_tp(model, w, t, p, lnfug=lnfug_w, root_type="stable")
        call tangent_plane_distance(z, lnfug_z, w, lnfug_w-log(P), f)
    end subroutine
end program
