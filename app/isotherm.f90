program main
    use yaeos, only: pr, ArModel, PengRobinson76, volume
    implicit none
    class(ArModel), allocatable :: eos

    integer, parameter :: nc = 2
    real(pr), dimension(nc) :: n, tc, pc, w
    real(pr), dimension(nc, nc) :: kij, lij

    real(pr) :: v
    real(pr) :: p0, pf, dp, p
    real(pr) :: t0, tf, dt, t
    integer :: i, j, n_p_points=500, n_t_points=5

    n = [0.3_pr, 0.7_pr]
    tc = [190._pr, 310._pr]
    pc = [14._pr, 30._pr]
    w = [0.001_pr, 0.03_pr]
    kij = reshape([0., 0.1, 0.1, 0.], [nc,nc])
    lij = kij / 2

    eos = PengRobinson76(tc, pc, w, kij, lij)

    p0 = 100
    pf = 0.1
    dp = (pf-p0)/n_p_points

    t0 = 150
    tf = 350
    dt = (tf-t0)/n_t_points

    do j=0, n_t_points - 1
        t = t0 + dt * j
        print *, "# ", t
        do i=0,n_p_points-1
            p = p0 + dp * i
            call volume(eos, n, p, t, v, root_type="stable")
            print *, v, p
        end do
        print *, ""
        print *, ""
    end do
end program main
