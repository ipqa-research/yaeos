program main
    use constants
    use mixing_rules
    use hyperdual_mod
    use cubic_eos, only: ClassicVdW, a_res
    use models, only: residual_helmholtz
    use pengrobinson76, only: PR76, setup_pr76
    ! use system, only: nc, kijs => kij, lijs => lij, PR76_factory, setup, bij,&
    !     m_ac => ac, m_b => b, m_k => k

    implicit none

    integer, parameter :: n=2
    real(pr) :: kij(n, n), lij(n, n)
    real(pr) :: z(n), p, v, t, phi_dual(n), phi_leg(n)
    real(pr) :: ares, dar(n+2), dar2(n+2, n+2)
    real(pr) :: tc(n), pc(n), w(n)

    real(pr) ::  aresl, arv, artv, arv2, arn(n), arvn(n), artn(n), arn2(n,n)

    real(pr) :: et, st, dt, at

    type(ClassicVdW), target :: mixrule
    type(PR76) :: model

    integer :: i, j, k

    kij(:, 1) = [0, 2]
    kij(:, 2) = [2, 0]
    lij = 0* 1.5*kij
    
    z = [0.3_pr, 0.7_pr]
    p = 2.0_pr
    v = 100.0_pr
    t = 500.0_pr
    
    tc = [190, 304]
    pc = [45, 74]
    w = [0.19, 0.1]

    call setup_MixingRule(mixrule, kij, lij)

    model%size = n
    model%names = ["A", "B"]

    call setup_pr76(model, n, tc, pc, w)

    ! nc = n
    ! call setup(nc, 2, 0, 0)
    ! kijs = kij
    ! lijs = lij

    ! call PR76_factory(z, tc_in=tc, pc_in=pc, w_in=w)

    ! m_ac = model%ac
    ! m_b = model%b
    ! m_k = model%k
    
    ! associate(b => model%b)
    ! do i=1,n
    !     do j=1, n
    !         bij(i, j) = (b(i) + b(j))/2 * (1 - lij(i, j))
    !         bij(i, i) = b(i)
    !     end do
    ! end do
    ! end associate

    model%mixrule => mixrule
    model%residual_helmholtz => a_res

    call setup_MixingRule(mixrule, kij, lij)
    
    ! call ArVnder(nc, 2, 1, z, v, t, aresl, arv, artv, arv2, arn, arvn, artn, arn2)
    call residual_helmholtz(model, z, v, t, ares, dar, dar2)

    print *, model%ac
    print *, model%b
    print *, model%k

    print *, "=================================================================="
    print *, "Ar: ", ares, aresl
    print *, "=================================================================="

    print *, "First order derivs: "
    do i=1,n
        print *, "dAr/dn", i, dar(i), arn(i)
    end do
    print *, "dAr/dv", dar(n+1), arv
    print *, " "
    
    print *, "=================================================================="

    print *, "Second order derivs: "
    print *, "dAr2/dv2", dar2(n+1, n+1), arv2
    print *, "dAr2/dtv", dar2(n+2, n+1), artv
    
    print *, "------------------------------------------------------------------"

    do i=1, n
        print *, "dAr2/dn2", i, dar2(i, i), arn2(i, i)
    end do

    print *, "------------------------------------------------------------------"

    do i=1, n
        print *, "dAr2/dvn", i, dar2(i, n+1), arvn(i)
    end do
    
    print *, "------------------------------------------------------------------"

    do i=1, n
        print *, "dAr2/dtn", i, dar2(i, n+2), artn(i)
    end do

    open(1, file="testfile")
    write(1, *) "v arv arn_1 arn_2 arn^2_{11} arn^2_{12} arn^2_{21} arn^2_{22} ares "

    call cpu_time(st)
    do i=1,100000
        v = real(i, kind=pr)/10
        ! call ArVnder(nc, 2, 1, z, v, t, aresl, arv, artv, arv2, arn, arvn, artn, arn2)
        call residual_helmholtz(model, z, v, t, ares, dar, dar2)

        ! write(1, "(F6.1, A)", advance="no") v, " "
        ! write(1, "(E10.4, A)", advance="no") arv - dar(n+1), " "
       
        do j=1,n
            ! write(1, "(E10.4, A)", advance="no") arn(j) - dar(j), " "
        end do

        do j=1,n
            do k=1,n
                ! write(1, "(E10.4, A)", advance="no") arn2(j, k) - dar2(j, k), " "
            end do
        end do

        ! write(1, "(E10.4, A)") ares - aresl
    end do
    close(1)
    call cpu_time(et)

    print *, et-st

end program main
