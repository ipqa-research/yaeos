module bench
    use yaeos, only: pr, R, Substances, AlphaSoave, CubicEoS, fugacity_vt, QMR
    implicit none

    real(pr), allocatable :: z(:), tc(:), pc(:), w(:), kij(:,:), lij(:,:)

contains

    subroutine set_vars(n)
        integer, intent(in) :: n
        integer :: i
        z = [(i, i=1,n)]
        z = z/sum(z)
        tc = z * 2
        pc = z * 100
        w  = z/10
        kij = reshape([(i,i=1,n**2)], [n,n])
        lij = kij/2
    end subroutine

    subroutine yaeos_run(n, dn)
        integer :: n
        logical :: dn
        type(CubicEoS) :: model
        type(Substances) :: composition
        type(AlphaSoave) :: alpha
        type(QMR) :: mixrule
        real(pr) :: lnfug(n), dlnphidp(n), dlnphidt(n), dlnphidn(n,n)

        integer :: i

        real(pr) :: v, t, p

        call set_vars(n)
        
        composition%tc = tc
        composition%pc = pc
        composition%w = w
        
        alpha%k = 0.37464_pr + 1.54226_pr * w - 0.26993_pr * w**2
        
        mixrule%k = kij
        mixrule%l = lij

        v = 1.0_pr
        t = maxval(composition%tc)

        model%ac = 0.45723553_pr * R**2 * composition%tc**2 / composition%pc
        model%b = 0.07779607_pr * R * composition%tc/composition%pc
        model%del1 = [1 + sqrt(2.0_pr)]
        model%del2 = [1 - sqrt(2.0_pr)]
        model%alpha = alpha
        model%components = composition
        model%mixrule = mixrule

        if (dn) then
            call fugacity_vt(model, z, V, T, P, lnfug, dlnPhidP, dlnphidT, dlnphidn)
        else
            call fugacity_vt(model, z, V, T, P, lnfug, dlnPhidP, dlnphidT)
        end if
    end subroutine
end module


program main
    use bench
    implicit none

    integer, parameter :: nevals=1e3
    integer :: i, n
    real(8) :: time, std, mean
    real(8) :: et, st


    do n=25,30
        time = 0
        std = 0
        mean = 0
        do i=1,nevals
            call cpu_time(st)
                call yaeos_run(n, .true.)
            call cpu_time(et)

            time = (et-st)*1e6
            mean = mean + time
            std = std + time**2
        end do

        mean = mean/nevals
        std = sqrt(1._pr/(nevals-1) * (std - nevals*mean**2))

        print *, n, mean, std
    end do
end program
