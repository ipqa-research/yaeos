module flashing
    use yaeos, only: pr, EquilibriaState, flash, PengRobinson76, ArModel, fugacity_tp
    implicit none

contains

    subroutine main()
        class(ArModel), allocatable :: model
        type(EquilibriaState) :: flash_result

        real(pr) :: tc(2), pc(2), w(2)

        real(pr) :: n(2), p, v, t, k0(2)
        integer :: iters, i, j

        print *, "FLASH EXAMPLE:"
        
        n = [0.4, 0.6]
        tc = [190.564, 425.12]
        pc = [45.99, 37.96]
        w = [0.0115478, 0.200164]
        model = PengRobinson76(tc, pc, w)

        P = 60
        t = 294
        k0 = (PC/P)*exp(5.373*(1 + w)*(1 - TC/T))
        print *, k0
        
        flash_result = flash(model, n, t=t, p_spec=p, k0=k0, iters=iters)
        print *, "X:", flash_result%x, sum(flash_result%x)
        print *, "Y:", flash_result%y, sum(flash_result%y)
    end subroutine
end module