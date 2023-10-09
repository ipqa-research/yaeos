! Dirty code to benchmark the implementation of PR76 compared to legacy code (faster) and a more dirty implementation.
module triangular_matrix
    use yaeos_constants, only: pr
    implicit none

contains

    pure function matrix_to_vector(m) result(v)
        real(pr), intent(in) :: m(:, :)
        real(pr) :: v(size(m, dim=1) * (size(m, dim=1)-1)/2 )
        integer :: i, j, n, nv

        n = size(m, dim=1)

        nv = 0
        do i=1, n
            do j=i+1,n
                nv = nv + 1
                v(nv) = m(i, j)
            end do
        end do
    end function
end module

module dirtypengrobinson76
    !-| Peng Robinson 76 Equation of State.
    use yaeos_constants, only: pr, R
    use yaeos_autodiff
    use yaeos_ar_models, only: set_ar_function
    use yaeos_interfaces, only: dual_property
    use triangular_matrix
    
    implicit none

    private
    public :: setup_dirty_pr76

    real(pr), allocatable :: kij(:, :), lij(:, :), kij_vect(:)
    real(pr), allocatable :: ac(:), b(:), k(:)
    real(pr), allocatable :: tc(:), pc(:), w(:)

    real(pr), parameter :: del1 = 1._pr + sqrt(2._pr)
    real(pr), parameter :: del2 = 1._pr - sqrt(2._pr)
contains
    subroutine setup_dirty_pr76(n, tc_in, pc_in, w_in, kij_in, lij_in)
        !-| Setup the enviroment to use the PengRobinson 76 Equation of State
        !   It uses the Cubic Van der Waals mixing rules
        integer :: n !! Number of components
        real(pr) :: tc_in(n)
        real(pr) :: pc_in(n)
        real(pr) :: w_in(n)
        real(pr) :: kij_in(n, n)
        real(pr) :: lij_in(n, n)

        tc = tc_in
        pc = pc_in
        w = w_in

        ac = 0.45723553_pr * R**2 * tc**2 / pc
        b = 0.07779607_pr * R * tc/pc
        k = 0.37464_pr + 1.54226_pr * w - 0.26993_pr * w**2

        kij = kij_in
        lij = lij_in

        kij_vect = matrix_to_vector(kij)

        call set_ar_function(arfun)
    end subroutine

    pure subroutine arfun(z, v, t, ar)
        type(hyperdual), intent(in) :: z(:), v, t
        type(hyperdual), intent(out) :: ar
    
        type(hyperdual) :: amix, a(size(z)), ai(size(z)), z2(size(z)), zij
        type(hyperdual) :: bmix
        type(hyperdual) :: b_v

        integer :: i, j

        a = 1.0_pr + k * (1.0_pr - sqrt(t/tc))
        ai = ac * a * a
        a = sqrt(ai)
        z2 = z * z

        amix = 0.0_pr
        bmix = 0.0_pr

        do i=1,size(z)-1
            do j=i+1,size(z)
                zij = z(i) * z(j)
                amix = amix + zij * (a(i) * a(j)) * (1 - kij(i, j))
                bmix = bmix + zij * (b(i) + b(j)) * (1 - lij(i, j))
             end do
        end do

        amix = 2 * amix + sum(z2 * ai)
        bmix = (bmix + sum(z2 * b)) / sum(z)
        
        b_v = bmix/v
        
        ar = (&
              - sum(z) * log(1.0_pr - b_v) &
              - amix / (R*t*bmix)*1.0_pr / (del1 - del2) &
              * log((1.0_pr + del1 * b_v) / (1.0_pr + del2 * b_v)) &
        )
    end subroutine
end module

module state
    use yaeos_constants, only: pr
    use yaeos_autodiff

    integer :: n
    real(pr), allocatable :: z(:), tc(:), pc(:), w(:), kij(:, :), lij(:, :), v, t
    type(hyperdual), allocatable :: z_d(:)
    
    type(hyperdual) :: v_d
    type(hyperdual) :: t_d
    type(hyperdual) :: ar_d
contains
    subroutine alloc(n)
        integer, intent(in) :: n
        integer :: i
        if (allocated(z)) deallocate(z)
        if (allocated(z_d)) deallocate(z_d)
        if (allocated(tc)) deallocate(tc)
        if (allocated(pc)) deallocate(pc)
        if (allocated(w)) deallocate(w)
        if (allocated(kij)) deallocate(kij)
        if (allocated(lij)) deallocate(lij)
        allocate(z(n), z_d(n), tc(n), pc(n), w(n), kij(n,n), lij(n,n))

        z = [(i, i=1,n)]
        z = z/sum(z)

        tc = z * 2
        pc = z * 100
        w  = z/10

        kij = 0.5
        lij = 0.5

        z_d = z
        v_d = 150.0_pr
        t_d = 500.0_pr
        v = 150
        t = 500
    end subroutine
end module

program test_pr76_two
    use state, only: n
    use yaeos_autodiff
    use yaeos_constants, only: pr
    use yaeos_thermo_properties, only: ln_phi
    use triangular_matrix, only: matrix_to_vector
    implicit none
    
    integer, parameter :: evals=1000
    integer :: i, ncomp

    real(pr) :: et, st
    real(pr) :: time_pr, time_dpr, time_lpr, desv_pr, desv_dpr, desv_lpr

    real(pr) :: kij(3,3)
    real(pr), allocatable :: k(:)

    print *, "n PR76 PR76err dirtyPR76 dirtyPR76err LPR76 LPR76err"

    do n=1,30
        call bench("PR76", time_pr, desv_pr)
        call bench("dirtyPR76", time_dpr, desv_dpr)
        call bench("leg_PR76", time_lpr, desv_lpr)
        print *, n , &
            time_pr * 1e6, desv_pr * 1e6, &
            time_dpr * 1e6, desv_dpr * 1e6, &
            time_lpr * 1e6, desv_lpr * 1e6
    end do
contains
    subroutine bench(model, time, desv)
        use state
        use pengrobinson76, only: setup_pr76
        use dirtypengrobinson76, only: setup_dirty_pr76
        use legacy_ar_models, only: PR76_factory, setup, lkij => kij, llij => lij
        use ar_interface, only: ar_fun
        character(len=*), intent(in) :: model
        real(pr), intent(out) :: time, desv

        call alloc(n)

        select case(model)
            case("PR76")
                call setup_pr76(n, tc, pc, w, kij, lij)
                call timeit(time, desv, sub)
            case("dirtyPR76")
                call setup_dirty_pr76(n, tc, pc, w, kij, lij)
                call timeit(time, desv, sub)
            case("leg_PR76")
                call setup(n, 2, 0, 1)
                lkij = 0
                llij = 0
                call PR76_factory(z, tc_in=tc, pc_in=pc, w_in=w)
                call timeit(time, desv, legsub)
        end select
    end subroutine

    subroutine sub()
        use state, only: z, v, t
        use yaeos_ar_models, only: residual_helmholtz, ar_fun
        real(pr) :: lnphi(size(z))
        
        call ln_phi(z, v, t, lnphi)
    end subroutine

    subroutine legsub()
        use state
        use legacy_ar_models, only: ArVnder, zTVTERMO

        real(pr) :: ar, arn(n), arnn(n, n)
        call zTVTermo(n, 5, 150.0_pr, z, 500.0_pr, ar, ar, arn, arn, arn, arnn)
    end subroutine

    subroutine timeit(time, desv, abs_sub)
        real(pr), intent(out) :: time, desv
        real(pr) :: times(evals)

        interface
            subroutine abs_sub()
            end subroutine
        end interface

        times = 0

        do i=1,evals
            call cpu_time(st)
                call abs_sub()
            call cpu_time(et)
            times(i) = (et - st)
        end do

        time = sum(times)/evals
        desv = sqrt(sum((times * times) - time**2)/evals)
    end subroutine
end program
