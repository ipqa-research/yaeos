module bench
    use yaeos, only: pr, R, Substances, AlphaSoave, CubicEoS, fugacity_vt, QMR, PengRobinson76, ArModel
    use hyperdual_pr76, only: PR76, setup_adiff_pr76 => setup
    use TapeRobinson, only: setup_taperobinson => setup_model, tape_model => model, ArModelTapenade
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

    subroutine yaeos_run(n, dn, model_name)
        integer :: n
        logical :: dn
        character(len=*), intent(in) :: model_name
        class(ArModel), allocatable :: model
        real(pr) :: lnfug(n), dlnphidp(n), dlnphidt(n), dlnphidn(n,n)

        integer :: i

        real(pr) :: v, t, p

        call set_vars(n)

        select case(model_name)
        case ("Analytic PR76")
            model = PengRobinson76(tc, pc, w, kij, lij)
        case ("Adiff PR76")
            model = setup_adiff_pr76(tc, pc, w, kij, lij)
        case ("Tape PR76")
            call setup_taperobinson(tc, pc, w, kij, lij)
            model = tape_model
        end select

        v = 1.0_pr
        t = 150._pr

        if (dn) then
            call fugacity_vt(model, z, V, T, P, lnfug, dlnPhidP, dlnphidT, dlnphidn)
        else
            call fugacity_vt(model, z, V, T, lnphip=lnfug)
        end if
    end subroutine

    subroutine run_bench(nmax, all_derivs, eos)
        integer, parameter :: nevals=1e3
        integer :: nmax
        logical :: all_derivs
        character(len=*), intent(in) :: eos
        real(8) :: time, std, mean
        real(8) :: et, st
        integer :: i, n
        time = 0
        std = 0
        mean = 0

        print *, "running: ", eos, "all derivs: ", all_derivs

        do n=1,nmax
            do i=1,nevals
                call cpu_time(st)
                    call yaeos_run(n, all_derivs, eos)
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
        integer :: n=20
        logical :: allderivs=.false.

        call run_bench(n, allderivs, "Analytic PR76")
        call run_bench(n, allderivs, "Tape PR76")
        call run_bench(n, allderivs, "Adiff PR76")

        allderivs = .true.
        
        call run_bench(n, allderivs, "Analytic PR76")
        call run_bench(n, allderivs, "Tape PR76")
        call run_bench(n, allderivs, "Adiff PR76")
    end subroutine

end module
