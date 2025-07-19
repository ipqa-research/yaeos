module bench
    use yaeos, only: pr, R, Substances, AlphaSoave, CubicEoS, &
        QMR, PengRobinson76, ArModel
    use hyperdual_pr76, only: PR76, setup_adiff_pr76 => setup
    use autodiff_tapenade_pr76_demo, only: setup_tapen_pr76 => setup_model
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
        kij = reshape([(0, i=1,n**2)], [n,n])
        lij = reshape([(0, i=1,n**2)], [n,n])
    end subroutine

    subroutine yaeos__run(n, dn, f_p, model_name)
        integer :: n
        logical :: dn
        logical :: f_p
        character(len=*), intent(in) :: model_name
        class(ArModel), allocatable :: model

        real(pr) :: Ar, ArV, ArT, ArVT, ArV2, ArT2
        real(pr) :: Arn(n), Arn2(n, n), ArTn(n), ArVn(n)

        real(pr) :: lnphi(n), lnphin(n,n)


        real(pr) :: v, t, p

        call set_vars(n)

        select case(model_name)
        case ("Analytic PR76")
            model = PengRobinson76(tc, pc, w, kij, lij)
        case ("Adiff PR76")
            model = setup_adiff_pr76(tc, pc, w, kij, lij)
        case ("Tape PR76")
            model = setup_tapen_pr76(tc, pc, w, kij, lij)
        end select

        v = 10.0_pr
        t = 250._pr
        p = 15

        if (dn) then
            call model%lnphi_vt(z, V, T, P=P, lnphi=lnPhi, dlnPhidn=lnphin)
        else
            call model%lnphi_vt(z, V, T, P=P, lnphi=lnPhi)
        end if

    end subroutine

    subroutine run_bench(nmax, all_derivs, f_p, eos)
        integer, parameter :: nevals=1e5
        integer :: nmax
        logical :: all_derivs
        logical :: f_p
        character(len=*), intent(in) :: eos
        real(pr) :: time, std, mean
        real(pr) :: et, st
        integer :: i, n
        time = 0
        std = 0
        mean = 0

        print "(//)"
        print *,  "#", eos, "all derivs: ", all_derivs
        do n=1,nmax
            do i=1,nevals
                call cpu_time(st)
                    call yaeos__run(n, all_derivs, f_p, eos)
                call cpu_time(et)

                time = (et-st)*1e6
                mean = mean + time
                std = std + time**2
            end do
            mean = mean/nevals
            std = sqrt(1._pr/(nevals-1) * abs(std - nevals*mean**2))

            print *, n, mean, std
        end do
    end subroutine
    
    subroutine main()
        integer :: n=30
        logical :: allderivs=.false.
        logical :: fug_p = .false.

        call run_bench(n, allderivs, fug_p, "Analytic PR76")
        call run_bench(n, allderivs, fug_p, "Tape PR76")
        call run_bench(n, allderivs, fug_p, "Adiff PR76")

        allderivs = .true.
        call run_bench(n, allderivs, fug_p, "Analytic PR76")
        call run_bench(n, allderivs, fug_p, "Tape PR76")
        call run_bench(n, allderivs, fug_p, "Adiff PR76")
    end subroutine

end module
