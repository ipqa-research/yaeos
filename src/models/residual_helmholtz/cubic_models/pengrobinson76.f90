module pengrobinson76
    use constants, only: pr, R
    use hyperdual_mod
    use cubic_eos, only: CubicEoS, a_classic, b_classic, c_classic, alloc
    use ar_models, only: size

    implicit none

    type, extends(CubicEOS) :: PR76
    end type

    interface setup_model
        module procedure :: setup_pr76
    end interface

    interface attractive_parameter
        module procedure :: a_pr
    end interface

    interface repulsive_parameter
        module procedure :: b_pr
    end interface
    
    interface del1_parameter
        module procedure :: del1_pr
    end interface
    
    interface del2_parameter
        module procedure :: del2_pr
    end interface

contains

    subroutine setup_pr76(model, n, tc, pc, w)
        use stdlib_optval, only: optval
        type(PR76) :: model
        integer :: n
        real(pr) :: tc(n), pc(n), w(n)

        call alloc(model%CubicEoS, n)

        model%del1 = 1._pr + sqrt(2.0_pr)
        model%del2 = 1._pr - sqrt(2.0_pr)

        model%tc = tc
        model%pc = pc
        model%w = w

        model%ac = 0.45723553_pr * R**2 * tc**2 / pc
        model%b = 0.07779607_pr * R * tc/pc
        model%k = 0.37464_pr + 1.54226_pr * w - 0.26993_pr * w**2
        model%c = 0
    end subroutine

    subroutine a_pr(model, z, v, t, a)
        type(PR76) :: model
        type(hyperdual) :: z(size(model))
        type(hyperdual) :: p, v, t
        type(hyperdual) :: a(size(model))

        call a_classic(model%CubicEoS, z, v, t, a)
    end subroutine

    subroutine b_pr(model, z, p, v, t, b)
        type(PR76) :: model
        type(hyperdual) :: z(size(model))
        type(hyperdual) :: p, v, t
        type(hyperdual) :: b(size(model))

        call b_classic(model%CubicEoS, z, v, t, b)
    end subroutine

    subroutine del1_pr(model, z, v, t, del1)
        type(PR76) :: model
        type(hyperdual) :: z(size(model))
        type(hyperdual) :: p, v, t
        type(hyperdual) :: del1(size(model))

        del1 = 1 + sqrt(2.0_pr)
    end subroutine

    subroutine del2_pr(model, z, v, t, del2)
        type(PR76) :: model
        type(hyperdual) :: z(size(model))
        type(hyperdual) :: v, t
        type(hyperdual) :: del2(size(model))

        del2 = 1 - sqrt(2.0_pr)
    end subroutine

end module
