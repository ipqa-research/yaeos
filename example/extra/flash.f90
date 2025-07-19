module flashing
    use yaeos, only: pr, EquilibriumState, flash, PengRobinson76, ArModel
    implicit none

contains

    subroutine main()
        class(ArModel), allocatable :: model
        type(EquilibriumState) :: flash_result

        real(pr) :: tc(2), pc(2), w(2)

        real(pr) :: n(2), p, v, t, k0(2)
        integer :: iters, i, j

        print *, "FLASH EXAMPLE:"
        
        n = [0.4, 0.6]
        Tc = [190.564, 425.12]
        Pc = [45.99, 37.96]
        w = [0.0115478, 0.200164]

        model = PengRobinson76(Tc, Pc, w)

        P = 60
        t = 294
        k0 = (Pc/P)*exp(5.373*(1 + w)*(1 - Tc/T))
        
        flash_result = flash(model, n, t=t, p_spec=p, k0=k0, iters=iters)
        print *, "X:", flash_result%x
        print *, "Y:", flash_result%y
    end subroutine
end module